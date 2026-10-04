import logging
import numpy as np
import random
import re
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# ---------------- annotation carried on a FASTA header ----------------
# The tag is matched at a word boundary so that an element whose NAME happens to contain
# "TSD=" cannot be misread as an annotation. The value is grabbed loosely, as a run of
# non-space, so that "TSD=" with nothing after it is caught as malformed rather than missed.
_TSD_TAG = re.compile(r"(?:\A|\s)TSD=(\S*)")
_TSD_VALUE = re.compile(r"\A(\d+)(?:-(\d+))?\Z")


def description_tail(record) -> str:
    """Whatever a FASTA header carries after the sequence ID, or "" if it carries nothing.

    Biopython sets description to the WHOLE title line, ID included, so a record parsed from
    ">copia#LTR/Copia TSD=5" has the id "copia#LTR/Copia" and the description
    "copia#LTR/Copia TSD=5". Stripping the ID back off leaves the annotation on its own,
    which is what can be copied onto a derived record without repeating the ID in its header.
    """
    if record.description.startswith(record.id):
        return record.description[len(record.id):].strip()
    return record.description.strip()


def parse_tsd_tag(description: str, record_id: str = ""):
    """The (min, max) TSD length a FASTA header asks for, or None if it asks for nothing.

    Grammar, both bounds inclusive::

        TSD=4       a fixed 4 bp duplication
        TSD=5-15    a length drawn uniformly from 5..15
        TSD=0       no duplication at all

    Both forms are needed because TSD length is a fixed, mechanistic property for some clades
    and a variable one for others. A cut-and-paste transposon's TSD is set by the stagger
    between its transposase's two cuts, so hAT and P duplicate 8 bp and Tc1/mariner 2 bp every
    time; a non-LTR element inserting by target-primed reverse transcription has no fixed
    stagger and produces a spread. Helitrons copy in by rolling-circle replication and
    duplicate nothing at all, hence TSD=0 -- a real value, not a missing one, which is why
    absence has to be spelled None rather than 0.

    A tag that is present but malformed raises rather than falling back: a header that meant
    to say something about TSD and failed is a typo to fix, not a default to accept.
    """
    match = _TSD_TAG.search(description)
    if match is None:
        return None
    where = f" in FASTA header {record_id!r}" if record_id else ""
    value = _TSD_VALUE.match(match.group(1))
    if value is None:
        raise ValueError(f"Malformed TSD tag 'TSD={match.group(1)}'{where}; "
                         f"expected TSD=<int> or TSD=<min>-<max>")
    low = int(value.group(1))
    high = int(value.group(2)) if value.group(2) is not None else low
    if high < low:
        raise ValueError(f"Inverted TSD tag 'TSD={match.group(1)}'{where}; "
                         f"the minimum must not exceed the maximum")
    return low, high


# ---------------- sampling TE insertions with min distance ----------------
def sample_TEins(regions, deletions, n: int, TEdistance: int, target_strands: list[str|None]):
    current_regions = np.array([[r[0], r[1], r[2], r[5]] for r in regions], dtype=object)
    
    def keep_TEdistance(regs, chrom, excl_start, excl_end):
        if len(regs) == 0:
            return regs
        
        reg_starts = regs[:, 1].astype(int)
        reg_ends = regs[:, 2].astype(int)
        
        mask = (regs[:, 0] != chrom) | (reg_ends <= excl_start) | (reg_starts >= excl_end)
        keep_regs = regs[mask]
        overlap_regs = regs[~mask]
        split_segments = []
        
        for r in overlap_regs:
            r_chrom, r_start, r_end, r_strand = r[0], int(r[1]), int(r[2]), r[3]
            if r_start < excl_start:
                split_segments.append([r_chrom, r_start, excl_start, r_strand])
            if r_end > excl_end:
                split_segments.append([r_chrom, excl_end, r_end, r_strand])
        
        if not split_segments:
            return keep_regs
        return np.vstack([keep_regs, np.array(split_segments, dtype=object)])

    # Subtract Deletions: exclude the deletion span itself plus a TEdistance buffer on each side.
    for d in (deletions or []):
        exclusion_start = d[1] - TEdistance
        exclusion_end = d[2] + TEdistance
        current_regions = keep_TEdistance(current_regions, d[0], exclusion_start, exclusion_end)

    # Sample n Positions
    positions = []
    for target_s in target_strands:
        # Filter by strand if a target is specified
        if target_s is not None:
            eligible_regs = current_regions[current_regions[:, 3] == target_s]
            if len(eligible_regs) == 0:
                logging.warning(
                    f"No regions remain on strand '{target_s}' (likely consumed by --TEdistance); "
                    f"falling back to any-strand placement for this insertion."
                )
                eligible_regs = current_regions
        else:
            eligible_regs = current_regions

        lengths = (eligible_regs[:, 2].astype(int) - eligible_regs[:, 1].astype(int))
        total_len = np.sum(lengths)

        if total_len <= 0:
            raise ValueError("TEdistance too large or genomic space exhausted.")

        # Select region
        r_val = np.random.randint(0, total_len)
        cum_lengths = np.cumsum(lengths)
        idx = np.searchsorted(cum_lengths, r_val, side='right')

        # Calculate coordinates
        prev_cum = cum_lengths[idx-1] if idx > 0 else 0
        offset = r_val - prev_cum
        chosen = eligible_regs[idx]
        actual_pos = int(chosen[1]) + offset
        
        positions.append((chosen[0], actual_pos, chosen[3]))

        # Update regions to enforce TEdistance for the next pos
        current_regions = keep_TEdistance(
            current_regions, chosen[0], 
            actual_pos - TEdistance, 
            actual_pos + TEdistance
        )

    return positions

# ---------------- adding background SVs ----------------
def bgSV(bedin:str, bedout:str, nSV:int, ins_ratio:float, fasta_in:str, fasta_out:str,
         margin:int = 30, seed=None):
    """
    Add background SVs (INS/DEL) to a BED file.
    Args:
        bedin: input BED file (only TE insertions/deletions)
        bedout: output BED file with added SVs
        nSV: number of SVs to add
        ins_ratio: ratio of insertions among the background SVs
        fasta_in / fasta_out: the pool, and the pool with the new insertion sequences added
        margin: background SVs stay this many bp clear of every TE event and of one another

    Insertions are points (start == end), as TErandom writes a TE insertion: an end one past
    the start made Simulate skip that reference base in every carrier. Positions are drawn per
    contig, from the gaps between that contig's events weighted by length -- the gaps were
    once pooled across contigs and every SV written to the last contig read -- and not only
    from the longest gaps, so a background SV can sit anywhere a TE event does not.
    """
    rnd = random.Random(seed) if seed is not None else random   # the global stream TErandom seeds
    TEs = []
    with open(bedin, "r") as fin:
        for line in fin:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            TEs.append(tuple(fields[:3]) + tuple(fields[3:]))
    events = [(f[0], int(f[1]), int(f[2])) for f in TEs]

    nINS = int(nSV * ins_ratio)
    nDEL = nSV - nINS
    SVmin, SVmax = 30, 300

    # Free gaps per contig, between that contig's events (with the margin), from its first event to its last.
    by_chrom = {}
    for chrom, s0, e0 in events:
        by_chrom.setdefault(chrom, []).append((s0 - margin, e0 + margin))
    gaps = []
    for chrom, spans in by_chrom.items():
        spans.sort()
        prev_end = spans[0][1]
        for s0, e0 in spans[1:]:
            if s0 > prev_end:
                gaps.append((chrom, prev_end, s0))
            prev_end = max(prev_end, e0)
    if not gaps:
        raise ValueError("No room between TE events for background SVs.")

    taken = {}
    def free(chrom, lo, hi):
        return all(hi + margin <= a or b + margin <= lo for a, b in taken.get(chrom, []))

    def draw(length):
        weights = [b - a - length for _, a, b in gaps]
        weights = [max(0, w) for w in weights]
        if not any(weights):
            raise ValueError(f"No gap between TE events holds a {length} bp background SV.")
        for _ in range(1000):
            chrom, a, b = rnd.choices(gaps, weights=weights, k=1)[0]
            lo = rnd.randint(a, b - length)
            if free(chrom, lo, lo + max(length, 1)):
                taken.setdefault(chrom, []).append((lo, lo + max(length, 1)))
                return chrom, lo
        raise ValueError("Too many background SVs for the room between TE events; ask for fewer.")

    new_bed, bgINS_seqs = [], []
    for idx in range(nDEL):
        del_len = rnd.randint(SVmin, SVmax)
        chrom, lo = draw(del_len)
        new_bed.append((chrom, lo, lo + del_len, f"bgDEL_{idx}_{del_len}"))
    for idx in range(nINS):
        ilen = rnd.randint(SVmin, SVmax - 1)
        name = f"bgINS_{idx}_{ilen}"
        bgINS_seqs.append((name, ''.join(rnd.choices('ATGC', k=ilen))))
        chrom, lo = draw(0)
        new_bed.append((chrom, lo, lo, name))

    records = list(SeqIO.parse(fasta_in, "fasta"))
    records.extend(SeqRecord(Seq(seq), id=name) for name, seq in bgINS_seqs)
    SeqIO.write(records, fasta_out, "fasta")

    rows = [tuple(t) for t in TEs] + [(c, str(a), str(b), n) for c, a, b, n in new_bed]
    rows.sort(key=lambda r: (r[0], int(r[1]), int(r[2])))
    with open(bedout, "w") as fout:
        for r in rows:
            fout.write("\t".join(map(str, r)) + "\n")

# ---------------- select TEs with the restrict of minimum number ----------------
def make_min_TE(TE_list: list, nMIN: int, nTE: int, TEtype: set, target_strands: list):
    # Organize TEs by family
    te_by_family = {}
    for entry in TE_list:
        family = entry[4]
        if family in TEtype:
            te_by_family.setdefault(family, []).append(entry)

    # Feasibility Check
    if len(te_by_family) * nMIN > nTE:
        raise ValueError(f"nMIN error: nMIN ({nMIN}) * num_families ({len(te_by_family)}) exceeds number of deletions requested ({nTE}).")

    selected_te = []
    
    # Satisfy nMIN for each family
    for family, members in te_by_family.items():
        if len(members) < nMIN:
            raise ValueError(f"Family {family} has fewer than {nMIN} members.")
        
        # Pick nMIN randomly to start
        selected_te.extend(random.sample(members, nMIN))
    
    nTE -= len(selected_te)
    if nTE == 0:
        return selected_te

    # nMIN forces us to pick at least nMIN members per family, regardless of strand. Consume the
    # matching slot from target_strands if possible, else burn a None slot, else pop arbitrarily.
    # Best-effort: the requested sense/antisense ratio may be skewed if many forced picks land on
    # the wrong strand.
    for te in selected_te:
        try:
            target_strands.remove(te[6])
        except ValueError:
            try:
                target_strands.remove(None)
            except ValueError:
                target_strands.pop()

    return selected_te + pick_stranded(
        [te for te in TE_list if te not in selected_te],
        nTE,
        target_strands
    )

def pick_stranded(TE_list: list, nTE: int, target_strands: list):
    selected_te = []
    n_either = target_strands.count(None)
    n_sense = target_strands.count("+")
    n_antisense = nTE - n_either - n_sense

    if n_sense:
        sense_pool = [x for x in TE_list if x[6] == "+"]
        if n_sense > len(sense_pool):
            selected_te += sense_pool
            n_either += n_sense - len(sense_pool)
        else:
            selected_te.extend(random.sample(sense_pool,n_sense))
    if n_antisense:
        antisense_pool = [x for x in TE_list if x[6] == "-"]
        if n_antisense > len(antisense_pool):
            selected_te += antisense_pool
            n_either += n_antisense - len(antisense_pool)
        else:
            selected_te.extend(random.sample(antisense_pool,n_antisense))
    if (n_sense or n_antisense) and n_either:
        TE_list = [x for x in TE_list if x not in selected_te]
    selected_te.extend(random.sample(TE_list,n_either))

    return selected_te