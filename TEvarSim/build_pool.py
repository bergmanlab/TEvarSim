import numpy as np
import random
import logging
import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from .utils import (sample_TEins, bgSV, make_min_TE, pick_stranded,
                    description_tail, parse_tsd_tag)

def is_dna(seq: Seq) -> bool:
    return set(str(seq)) <= set("ATGCN")

def check_TEid(seq: str) -> bool:
    return "#" in seq and "/" in seq

def check_output_file(output_path):
    if os.path.exists(output_path):
        raise FileExistsError(f"Output file '{output_path}' already exists. Please choose a different name.")

def CHRnorm(chr,chr_list):
    ### normalize chr according the style of chrs in fasta index
    if chr in chr_list:
        return chr
    chr_permutations = {
        chr.replace("chr",""),
        chr.replace("chr","CHR"),
        chr.replace("CHR","chr"),
        chr.replace("CHR",""),
        f"chr{chr}",
        f"CHR{chr}"
    }
    existing_chrs = set(iter(chr_list))
    matches = chr_permutations & existing_chrs
    if len(matches) == 1:
        return matches.pop()
    raise ValueError(f"Chromosome {chr} not found in chromosome list: {','.join(iter(chr_list))}")

class TEtype(list):
    def __init__(self,te_list):
        if te_list is None:
            te_list = []
        super().__init__(te_list)
    def __contains__(self, other):
        if not self:
            return True
        else:
            return super().__contains__(other)

# An element is only full-length if its INTERNAL region is. The LTR-I-LTR fragment pattern
# alone is not enough: a decayed relic keeps recognisable LTRs long after its internal region
# has gone, so RepeatMasker still lays down LTR / I / LTR and the pattern matches. One such
# relic in DG4005 (chrVII_RagTag:400163) retains 355bp of a 5249bp internal region -- 6.8% --
# and without this would be an eligible excision candidate, producing an "excision" of a
# 1.6kb pseudo-element leaving a 280bp "solo LTR".
MIN_INTERNAL_COVERAGE = 0.5


def ltr_element_core(te_info):
    """The LTR-I-LTR span of a full-length LTR element, and its 5' LTR length.

    RepeatMasker groups fragments by its ID column, and that group can hold more than the
    element: an adjacent solo LTR gets folded in whenever it resembles the element's own
    LTRs closely enough for RepeatMasker to join them. Whether it does depends on the
    library -- with Ty2 present the introduced Ty1c element's LTRs match TY2-LTR at 4%
    while the neighbouring relic matches TY1-LTR at 20%, so they stay apart; with a
    Ty1-only library both become TY1-LTR at 12.6% and 20.1% and RepeatMasker joins them.
    Taking min(start)/max(end) over the group then overshoots the element by that
    neighbour, and an excision built from it removes 309bp too much.

    An element is LTR-I-LTR: three fragments, not four. So bound it structurally instead --
    the last LTR before the internal region, and the first LTR after it. That is the same
    answer under either library, since it reads the element's own composition rather than
    how well its LTRs happen to score.

    Returns (start, end, ltr_len), or None if the group is not LTR-I-LTR shaped.
    """
    bounding = bounding_ltrs(te_info)
    if bounding is None:
        return None
    five_prime, three_prime = bounding
    return five_prime[5], three_prime[6], five_prime[6] - five_prime[5]


def bounding_ltrs(te_info):
    """The two LTR fragments an element's internal region sits between, or None.

    The last LTR before the internal region and the first one after it -- not the first and
    last fragments in RepeatMasker's ID group, which may include a neighbour it folded in.
    """
    ordered = sorted(te_info, key=lambda record: record[5])
    internal = [i for i, record in enumerate(ordered) if record[9][-2:] == "-I"]
    if not internal:
        return None
    before = [i for i in range(internal[0]) if ordered[i][9][-4:] == "-LTR"]
    after = [i for i in range(internal[-1] + 1, len(ordered)) if ordered[i][9][-4:] == "-LTR"]
    if not before or not after:
        return None
    return ordered[before[-1]], ordered[after[0]]


def ltrs_agree_on_strand(te_info):
    """Are this element's two bounding LTRs on the same strand?

    An LTR retrotransposon's LTRs are DIRECT repeats -- same sequence, same orientation --
    because the element is one transcriptional unit and both LTRs are copies of the same
    end. Two LTRs facing each other are not one element's ends; they are an inverted pair,
    or two unrelated LTRs RepeatMasker bridged across whatever lies between them, and the
    LTR-I-LTR shape is coincidental rather than structural.

    It also rules the span out of excision on mechanism: LTR-LTR recombination needs the
    two repeats aligned in the same direction. An inverted pair recombines to an inversion,
    not to a solo LTR.

    Returns True when the LTRs cannot be identified, leaving that judgement to
    ``bounding_ltrs``, which is the check that owns it.
    """
    bounding = bounding_ltrs(te_info)
    if bounding is None:
        return True
    five_prime, three_prime = bounding
    return five_prime[8] == three_prime[8]


# How far apart the two LTRs' divergences may be, in percentage points, for the element to
# be excisable.
MAX_LTR_DIVERGENCE_DIFFERENCE = 5.0


def ltrs_can_recombine(te_info, tolerance=MAX_LTR_DIVERGENCE_DIFFERENCE):
    """Are this element's two LTRs alike enough to recombine with each other?

    An excision is LTR-LTR recombination, so it needs homology between the two LTRs -- this
    is a requirement of the mechanism, not a data-quality filter. An element whose LTRs have
    drifted apart cannot recombine them, however well each one matches a consensus.

    Divergence from the consensus is the available proxy: two LTRs that are both close to it
    are close to each other, and a large gap between them means at least one has drifted.
    It is only a proxy -- two LTRs equally divergent from the consensus could still differ
    from each other -- so the tolerance is loose, meant to reject the clearly implausible
    rather than to measure recombination potential.

    Returns True when divergence cannot be read, so a malformed column never silently drops
    an element.
    """
    bounding = bounding_ltrs(te_info)
    if bounding is None:
        return False
    try:
        divergences = [float(record[1]) for record in bounding]
    except (TypeError, ValueError):
        return True
    return abs(divergences[0] - divergences[1]) <= tolerance


def internal_is_full_length(te_info, threshold=MIN_INTERNAL_COVERAGE):
    """Does this element's internal region cover enough of its consensus to be full-length?

    Read from RepeatMasker's own repeat-coordinate columns, so no consensus FASTA is needed:
    each fragment records where in the consensus it matched and how much was left over, which
    together give the consensus length. Fragments are summed, since one internal region is
    often reported in several pieces.

    Returns True when the coverage cannot be determined, so a missing or malformed coordinate
    never silently drops an element from the pool.
    """
    internal = [record for record in te_info if record[9][-2:] == "-I"]
    if not internal:
        return False
    covered = 0
    consensus = 0
    for record in internal:
        begin, end, left = record[11], record[12], record[13]
        if begin is None or end is None:
            return True
        covered += max(0, end - begin + 1)
        if left is not None:
            consensus = max(consensus, end + abs(left))
    if not consensus:
        return True
    return covered / consensus >= threshold


class RandomTE:
    def __init__(self, args):
        self.pool_fasta = args.outprefix + ".fa"
        self.out_fasta = args.outprefix + ".bgSV.fa"
        self.DELfile = args.existingTEs
        self.EXCfile = getattr(args, "excCandidates", None)
        self.sense_strand_ratio = args.sense_strand_ratio
        self.CHR = {}
        if not os.path.isfile(f"{args.ref}.fai"):
            raise FileNotFoundError(f"Reference fasta index required: {args.ref}.fai not found.")
        with open(f"{args.ref}.fai") as fin:
            for line in fin:
                line = line.strip().split("\t")
                if line:
                    self.CHR[line[0]] = int(line[1])
        if args.CHR:
            chosen_chrs = [CHRnorm(chr,self.CHR) for chr in args.CHR.split(",")]
            self.CHR = {chr:chr_len for chr,chr_len in self.CHR.items() if chr in chosen_chrs}
        if args.regions:
            try:
                with open(args.regions) as fin:
                    regions = fin.read().strip().split("\n")
            except FileNotFoundError:
                raise FileNotFoundError(f"Regions file {args.regions} not found.")
            self.regions = []
            for region in regions:
                region = region.split("\t")
                if region[0] in self.CHR:
                    region[1] = int(region[1])
                    region[2] = int(region[2])
                    if len(region) > 5:
                        self.regions.append(region)
                    else:
                        self.regions.append(region[:3]+[".",".","+"])
                        self.regions.append(region[:3]+[".",".","-"])
            if regions and not self.regions:
                raise ValueError("Chromosomes in regions file not found in reference.")
        else:
            # Skip position 0: simulate.py reads seq[start-1], so start=0 wraps to seq[-1].
            # End is half-open per BED convention.
            self.regions = [[chr, 1, chr_len, ".", ".", strand]
                for strand in ["+", "-"]
                for chr, chr_len in self.CHR.items()]
        if args.exclude:
            try:
                with open(args.exclude) as fin:
                    excluded_regions = fin.read().strip().split("\n")
            except FileNotFoundError:
                raise FileNotFoundError(f"Excluded regions file {args.exclude} not found.")
            for excluded_region in excluded_regions:
                excluded_region = excluded_region.split("\t")
                excluded_region[1] = int(excluded_region[1])
                excluded_region[2] = int(excluded_region[2])
                new_regions = []
                for region in self.regions:
                    if region[0] == excluded_region[0] and (len(excluded_region) < 6 or region[5] == excluded_region[5]):
                        if excluded_region[1] < region[2] and excluded_region[2] > region[1]:
                            if region[1] < excluded_region[1]:
                                new_regions.append((region[0],region[1],excluded_region[1],".",".",region[5]))
                            if region[2] > excluded_region[2]:
                                new_regions.append((region[0],excluded_region[2],region[2],".",".",region[5]))
                            continue
                    new_regions.append(region)
                self.regions = new_regions

        self.nINS = args.nINS
        self.nDEL = args.nDEL
        self.nEXC = args.nEXC  # number of LTR-LTR recombinations (excisions)
        self.TEtype = TEtype(args.TEtype)
        self.prefix = args.outprefix
        self.DELlen = args.DELlen
        self.TEdistance = args.TEdistance
        self.nSV = args.nSV # number of background SVs
        self.nMIN = args.nMIN # minimum number of TEs for each family
        self.random_seed = args.seed
        if self.random_seed is not None:
            np.random.seed(self.random_seed)
            random.seed(self.random_seed)

        # Strand targeting applies to placed insertions and selected deletions; excisions inherit
        # the orientation of the reference element, so they are not part of the strand budget.
        n_strand = self.nDEL + self.nINS
        if args.sense_strand_ratio is not None:
            n_sense = round(n_strand*args.sense_strand_ratio)
            self.target_strands = (["+"] * n_sense) + (["-"] * (n_strand - n_sense))
            np.random.shuffle(self.target_strands)
        else:
            self.target_strands = [None] * n_strand

    def _run(self):
        # 1. parse DEL file (populates self.DEL, and self.EXC candidates for RepeatMasker input)
        self.parse_DEL()
        # 2. select excisions (full-length LTR elements) first, then drop those loci from the
        #    deletion candidate pool so a locus is never used as both a DEL and an EXC.
        self.select_EXC()
        # 3. select deletions from the remaining candidates
        if self.nDEL > len(self.DEL):
            raise ValueError(
                f"Requested {self.nDEL} deletions but only {len(self.DEL)} available in --existingTEs "
                f"after filtering (and after removing excision loci). Lower --nDEL, or relax --TEtype/--DELlen."
            )
        logging.info(f"Generating {self.nINS} INS, {self.nDEL} DEL and {self.nEXC} EXC for chromosome(s) {','.join(self.CHR.keys())}")
        nplaced = self.nDEL + self.nINS
        if self.nMIN > 0 and self.nMIN >= nplaced:
            raise ValueError(f"minimum number of a TE family ({self.nMIN}) should be less than nDEL+nINS ({nplaced})")
        if self.nMIN > 0:
            self.DEL = make_min_TE(self.DEL, self.nMIN, self.nDEL, self.TEtype, self.target_strands[:self.nDEL])
        else:
            self.DEL = pick_stranded(self.DEL, self.nDEL, self.target_strands[:self.nDEL])
        # 4. place insertions (avoiding both selected DEL and EXC spans)
        self.parse_TEpool()
        # 5. Generate BED file
        self.build_bed()
        logging.info(f"Generated TE BED file: {self.prefix}.bed")
        # 6. Add background SVs if specified
        if self.nSV > 0:
            bedin = self.prefix + ".bed"
            bedout = self.prefix + ".bgSV.bed"
            bg_ins_ratio = self.nINS / (self.nINS + self.nDEL) if (self.nINS + self.nDEL) > 0 else 0.5
            bgSV(bedin, bedout, self.nSV, bg_ins_ratio, self.pool_fasta, self.out_fasta)
            logging.info(f"Generated BED file and fasta file with background SVs are: {bedout} and {self.out_fasta}")

    def select_EXC(self):
        """Pick nEXC full-length LTR elements to excise, and remove their loci from self.DEL."""
        if not hasattr(self, "EXC"):
            self.EXC = []
        if self.nEXC == 0:
            self.EXC = []
            return
        if self.nEXC > len(self.EXC):
            raise ValueError(
                f"Requested {self.nEXC} excisions but only {len(self.EXC)} full-length LTR elements "
                f"available in --existingTEs after filtering. Lower --nEXC, or relax --TEtype/--DELlen."
            )
        chosen = random.sample(self.EXC, self.nEXC)
        # A deletion candidate that overlaps a chosen element is that element (with --excCandidates the
        # deletion's span carries the TSD copy the excision's does not), so it cannot be deleted as well.
        by_chrom = {}
        for c in chosen:
            by_chrom.setdefault(c[0], []).append((c[1], c[2]))
        self.DEL = [d for d in self.DEL if not any(d[1] < e and s < d[2] for s, e in by_chrom.get(d[0], []))]
        self.EXC = chosen

    def build_bed(self):
        # merge INS, DEL and EXC
        merged = self.INS + self.DEL + self.EXC
        merged.sort()
        # bedfile output
        bed_name = self.prefix + ".bed"
        check_output_file(bed_name)
        with open(bed_name, "w") as f:
            for chrom, start, end, teID, _cls, _typ, strand, *rest in merged:
                cols = [chrom, str(start), str(end), teID, ".", strand]
                if rest:  # EXC records carry the LTR length as a 7th column
                    cols.append(str(rest[0]))
                f.write("\t".join(cols) + "\n")
        logging.info(f"Generated BED file: {bed_name}")
    
    def parse_DEL(self):
        self.DEL = []
        self.EXC = []  # full-length LTR elements eligible for excision (RepeatMasker input, or --excCandidates)
        self.parse_EXC_candidates()
        if not self.DELfile:
            return
        # process file
        ext = self.DELfile.split(".")[-1]
        if ext == "out":
            # Parse RepeatMasker input
            self.parse_DEL_repeatmasker()
            return
        elif ext == "txt":
            # Parse UCSC input
            def readline(line):
                if line.startswith("#"):
                    return
                fields = line.strip().split('\t')
                #       chrom,      repClass    start           end             strand      name        class_fam
                return  fields[5],  fields[12], int(fields[6]), int(fields[7]), fields[9],  fields[10], fields[12]
        elif ext == "bed":
            # Parse bed input
            def readline(line):
                if not line.strip():
                    return
                fields = line.strip().split("\t")
                name, class_fam = fields[3].split("#")
                #       chrom,      repClass                  start           end             strand     name    class_fam
                return fields[0],   class_fam.split("/")[0],  int(fields[1]), int(fields[2]), fields[5],  name,   class_fam
        else:
            raise ValueError(f"Input file format not recognized for --existingTEs: {ext}")
        with open(self.DELfile) as f:
            for line in f:
                try:
                    chrom, repClass, start, end, strand, name, class_fam = readline(line)
                except ValueError as e:
                    if "not enough values to unpack" in str(e):
                        continue
                    raise
                if end - start < self.DELlen:
                    continue
                if repClass in self.TEtype:
                    teID = f"DEL-{chrom}-{start}-{end}-{class_fam}-{name}"
                    self.DEL.append((chrom, start, end, teID, repClass,"DEL",strand))
    
    def parse_EXC_candidates(self):
        """Excision candidates from --excCandidates: full-length LTR elements whose span is known exactly,
        given as the element alone (no TSD) with its 5' LTR length in a 7th column. The excision leaves
        that LTR, between the element's two target-site copies, as RepeatMasker-derived candidates do."""
        if not self.EXCfile:
            return
        with open(self.EXCfile) as f:
            for line in f:
                if not line.strip() or line.startswith("#"):
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 7:
                    raise ValueError(f"--excCandidates needs 7 columns (the 7th the 5' LTR length): {line.strip()}")
                chrom, start, end = fields[0], int(fields[1]), int(fields[2])
                name, class_fam = fields[3].split("#")
                repClass = class_fam.split("/")[0]
                strand, ltr_len = fields[5], int(fields[6])
                if repClass not in self.TEtype or not 0 < ltr_len < end - start:
                    continue
                exc_id = f"EXC-{chrom}-{start}-{end}-{ltr_len}-{class_fam}-{name}"
                self.EXC.append((chrom, start, end, exc_id, repClass, "EXC", strand, ltr_len))

    def parse_DEL_repeatmasker(self):
        repeatmasker_records = {}
        with open(self.DELfile) as f:
            for line in f:
                if not line.strip() or not line.strip()[0].isdigit(): 
                    continue
                line = line.strip().split()
                if line[8] == "C":
                    line[8] = "-"
                for i in (5,6,14):
                    line[i] = int(line[i])
                # Repeat (consensus) coordinates, used to tell a full-length element from a
                # relic that has kept its LTRs. RepeatMasker writes these as
                # "begin end (left)" on +, and "(left) end begin" on C, so put them back in
                # order before use. Parentheses are stripped; a malformed field is left None.
                begin, end, left = line[11], line[12], line[13]
                if line[8] == "-":
                    begin, left = left, begin
                try:
                    line[11] = int(begin.strip("()"))
                    line[12] = int(end.strip("()"))
                    line[13] = int(left.strip("()"))
                except (ValueError, AttributeError):
                    line[11] = line[12] = line[13] = None
                if line[14] not in repeatmasker_records:
                    repeatmasker_records[line[14]] = []
                repeatmasker_records[line[14]].append(line)
        for te_info in repeatmasker_records.values():
            chrom = te_info[0][4]
            start = min(record[5] for record in te_info)
            end = max(record[6] for record in te_info)
            if end - start < self.DELlen:
                continue
            match_name = {}
            match_class_fam = {}
            match_strand = {"+":0,"-":0}
            for record in te_info:
                if record[9] not in match_name:
                    match_name[record[9]] = 0
                if record[10] not in match_class_fam:
                    match_class_fam[record[10]] = 0
                record_len = record[6] - record[5]
                match_name[record[9]] += record_len
                match_class_fam[record[10]] += record_len
                match_strand[record[8]] += record_len
            if match_strand["+"] > match_strand["-"]:
                strand = "+"
            else:
                strand = "-"
            is_full_ltr = False
            if len(match_name) == 1:
                name = next(iter(match_name))
                class_fam = next(iter(match_class_fam))
            else:
                ty_i = {n[:-2]:n_len for n,n_len in match_name.items() if n[-2:] == "-I"}
                if ty_i:
                    match_name = ty_i
                match_name = dict(sorted(match_name.items(),key=lambda item: item[1], reverse=True))
                match_class_fam = dict(sorted(match_class_fam.items(),key=lambda item: item[1], reverse=True))
                match_max = next(iter(match_name.values()))
                matches = [n for n in match_name if match_name[n] >= 0.2*match_max]
                name = "|".join(matches)
                # Full-length LTR element: an internal (-I) region flanked by two LTR (-LTR)
                # fragments. These are the elements eligible for LTR-LTR recombination (excision).
                is_full_ltr = (bool(ty_i) and te_info[0][9][-4:] == "-LTR"
                               and te_info[-1][9][-4:] == "-LTR"
                               and internal_is_full_length(te_info)
                               and ltrs_agree_on_strand(te_info))
                if is_full_ltr:
                    name += "-FULL"
                class_fam = next(iter(match_class_fam))
            repClass = class_fam.split("/")[0]
            # An LTR element is bounded by its own LTR-I-LTR structure, not by the extent of
            # RepeatMasker's ID group, which can have an adjacent solo LTR folded into it
            # depending on the library. Narrow the span here, before anything is built from
            # it, so the deletion and the excision describe the same element: deleting a
            # neighbouring relic along with the element is the same error as excising it.
            ltr_len = None
            if is_full_ltr:
                core = ltr_element_core(te_info)
                if core is not None:
                    start, end, ltr_len = core
                    # The narrowed element can fall under --DELlen where the ID group did not.
                    if end - start < self.DELlen:
                        continue
            if repClass in self.TEtype:
                teID = f"DEL-{chrom}-{start}-{end}-{class_fam}-{name}"
                self.DEL.append((chrom, start, end, teID, repClass,"DEL",strand))
                # The solo LTR left behind spans [start, start+ltr_len). Excision is
                # recombination between the two LTRs, so they have to be alike enough to
                # recombine -- an element whose LTRs have drifted apart is still deletable,
                # just not excisable.
                if (ltr_len is not None and 0 < ltr_len < end - start
                        and ltrs_can_recombine(te_info)):
                    exc_id = f"EXC-{chrom}-{start}-{end}-{ltr_len}-{class_fam}-{name}"
                    self.EXC.append((chrom, start, end, exc_id, repClass,
                                     "EXC", strand, ltr_len))
    
    def parse_TEpool(self):
        self.INS = []
        if self.nINS == 0:
            return
        records = list(SeqIO.parse(self.pool_fasta, "fasta"))
        # Guard against lst[-0:] == lst[0:] (whole list); explicitly take the last nINS entries.
        ins_strands = self.target_strands[-self.nINS:]
        # Insertions must avoid both selected deletion and excision spans.
        exclusion = self.DEL + self.EXC
        INSpos = sample_TEins(self.regions, exclusion, self.nINS, TEdistance=self.TEdistance, target_strands=ins_strands)
        for record, (chrom, pos, strand) in zip(records,INSpos):
            self.INS.append((chrom, pos, pos, record.id, record.id.split("/")[1], "INS", strand))

class TEPoolBuilder:
    def __init__(self, args):
        # base
        self.consensus = args.consensus
        self.nINS = args.nINS
        self.prefix = args.outprefix
        # SNP and INDEL
        self.snp_rate = args.snp_rate
        self.indel_rate = args.indel_rate
        self.indel_ins = args.indel_ins
        self.indel_geom_p = args.indel_geom_p
        # truncate
        self.truncated_ratio = args.truncated_ratio
        self.truncated_max_length = args.truncated_max_length
        # polyA
        self.polyA_ratio = args.polyA_ratio
        self.polyA_min = args.polyA_min
        self.polyA_max = args.polyA_max
        # other
        self.random_seed = args.seed
        if self.random_seed is not None:
            np.random.seed(self.random_seed)
        logging.basicConfig(level=logging.INFO, format="[%(levelname)s] %(message)s")
        
        # data stored
        self.BASES = ['A', 'C', 'G', 'T']
        self.CHOICES_DICT = {b: [c for c in self.BASES if c != b] for b in self.BASES}
        self.current_seq = []
        self.seqID_suffix = ""

    def _run(self):
        records = list(SeqIO.parse(self.consensus, "fasta"))
        if not records:
            raise ValueError(f"No sequences found in consensus FASTA: {self.consensus}")
        logging.info(f"Found {len(records)} sequences in consensus FASTA.")
        # check if all sequence are legal
        parsed_records = []
        potential_LTRs = {}
        for record in records:
            record.seq = record.seq.upper()
            if not is_dna(record.seq):
                raise ValueError(f"Invalid characters found in {record.id}")
            if not check_TEid(record.id):
                raise ValueError(f"Invalid format in TE ID: {record.id}. Please follow the format '>TEname#class/superfamily', e.g., '>AluY#SINE/Alu'")
            # The library is the only file that carries an element's annotation, and this is
            # the only command that reads it, so a malformed tag has to be caught here or it
            # surfaces two commands later as a silent fallback to the global TSD flags.
            parse_tsd_tag(description_tail(record), record.id)
            if record.id.split("#")[0][-2:] == "-I" or record.id.split("#")[0][-4:] == "-LTR":
                id, class_fam = record.id.split("#")
                ltr_name, ltr_part = id.rsplit("-",1)
                if ltr_name not in potential_LTRs:
                    potential_LTRs[ltr_name] = {"class_fam":class_fam}
                potential_LTRs[ltr_name][ltr_part] = record
            else:
                parsed_records.append(record)
        if potential_LTRs:
            for ltr_name, ltr_parts in potential_LTRs.items():
                if "I" in ltr_parts and "LTR" in ltr_parts:
                    ltr_rec = ltr_parts["LTR"]
                    internal_rec = ltr_parts["I"]
                    
                    new_id = f"{ltr_name}-FULL#{ltr_parts['class_fam']}"
                    # The synthesised full-length element inherits the LTR part's annotation
                    # rather than describing itself. A target site duplication is made by the
                    # integration machinery cutting the host at the element's ends, so it is a
                    # property of the LTR, and the -LTR record is where a library states it.
                    # (The description this replaces said "Merged LTR-I-LTR sequence for X",
                    # which never reached disk anyway: the pool writer below discarded it.)
                    new_record = SeqRecord(
                        ltr_rec.seq + internal_rec.seq + ltr_rec.seq,
                        id=new_id,
                        description=description_tail(ltr_rec) or description_tail(internal_rec)
                    )
                    
                    parsed_records.append(new_record)
                else:
                    # An -I with no -LTR, or the reverse: emit the parts as themselves. The
                    # dict also holds the class/family under "class_fam", which is a string
                    # and not a record, so it has to be skipped -- appending it put a str
                    # into the record list and crashed on .seq one loop later.
                    for part, part_rec in ltr_parts.items():
                        if part != "class_fam":
                            parsed_records.append(part_rec)
            records = parsed_records
                    

        # generate random masks for truncation and polyA
        truncated_num = np.random.random(self.nINS)
        polyA_num = np.random.random(self.nINS)
        trunc_mask = truncated_num < self.truncated_ratio
        polyA_mask = polyA_num < self.polyA_ratio
        mutate = np.select(
            [trunc_mask & polyA_mask, trunc_mask, polyA_mask],
            [3, 2, 1],
            default=0
        )

        # map numbers to corresponding functions
        func_map = {
            0: lambda: self.INDEL_mutate().SNP_mutate(),
            1: lambda: self.INDEL_mutate().SNP_mutate().apply_polyA(),
            2: lambda: self.INDEL_mutate().SNP_mutate().apply_truncate(),
            3: lambda: self.INDEL_mutate().SNP_mutate().apply_truncate().apply_polyA()
        }

        out_records = []
        for j in mutate:
            idx = np.random.randint(0, len(records))
            record = records[idx]
            self.current_seq = list(str(record.seq))
            self.seqID_suffix = ""  # reset suffix each loop
            func_map[j]()
            new_id = f"{record.id}_{self.seqID_suffix}"
            # The pool is the only file Simulate reads -- it never opens --consensus -- so
            # anything the library said about an element has to survive into the pool header
            # or it is lost. The ID picks up the modification suffix and so is no longer a
            # prefix of the library description, hence carrying the tail separately.
            new_record = SeqRecord(Seq("".join(self.current_seq)), id=new_id,
                                   description=description_tail(record))
            out_records.append(new_record)

        # output
        output_fasta = self.prefix + ".fa"
        check_output_file(output_fasta)
        SeqIO.write(out_records, output_fasta, "fasta")
        logging.info(f"Generated {len(out_records)} sequences -> {output_fasta}")
    

    def apply_truncate(self):
        max_trunc = int(len(self.current_seq) * self.truncated_max_length)
        if max_trunc < 1:
            return self
        trunc_len = np.random.randint(1, max_trunc + 1)  # inclusive upper bound
        self.current_seq = self.current_seq[trunc_len:]
        self.seqID_suffix += f"{trunc_len}truncate"
        return self

    def apply_polyA(self):
        polyA_len = np.random.randint(self.polyA_min, self.polyA_max + 1)
        self.current_seq.extend(["A"] * polyA_len)
        self.seqID_suffix += f"{polyA_len}polyA"
        return self
    
    def SNP_mutate(self):
        L = len(self.current_seq)
        n_snp = np.random.poisson(self.snp_rate * L)
        if n_snp < 1:
            return self
        self.seqID_suffix += f"{n_snp}SNP"
        # SNP positions, ensure not exceeding sequence length
        snp_positions = np.random.choice(L, size=min(n_snp, L), replace=False)
        for pos in snp_positions:
            current_base = self.current_seq[pos]
            if current_base in self.CHOICES_DICT:
                self.current_seq[pos] = np.random.choice(self.CHOICES_DICT[current_base])
        return self
    
    def INDEL_mutate(self):
        # total INDEL number
        L = len(self.current_seq)
        n_indel = np.random.poisson(self.indel_rate * L)
        if n_indel < 1:
            return self
        self.seqID_suffix += f"{n_indel}INDEL"
        if n_indel == 0:
            return self
        # generate INDEL lengths
        indel_len_list = np.random.geometric(self.indel_geom_p, n_indel)
        indel_len_list = np.clip(indel_len_list, a_min=None, a_max=30)
        # select positions for INDELs
        indel_positions = np.random.choice(L, size=min(n_indel, L), replace=False)
        for pos, indel_len in zip(sorted(indel_positions, reverse=True), indel_len_list):
            if np.random.random() < self.indel_ins:
                ins_seq = [np.random.choice(self.BASES) for _ in range(indel_len)]
                self.current_seq[pos:pos] = ins_seq
            else:
                del_end = min(pos + indel_len, len(self.current_seq))
                del self.current_seq[pos:del_end]
        return self

def run(args):
    # generate TE pool
    TEPoolBuilder(args)._run()
    # generate bed file with random TE insertions
    RandomTE(args)._run()


