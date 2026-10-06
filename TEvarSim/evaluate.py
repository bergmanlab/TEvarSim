'''
``tevarsim Evaluate`` -- score a prediction against every simulated locus at once.

``Compare`` answers "how well was *this* genome called?". It takes one truth sample and one
prediction sample and reports site-level precision/recall plus a genotype accuracy for that
pair, so benchmarking a whole simulation means looping it over every genome and stitching
the per-genome outputs back together by hand.

``Evaluate`` turns the question around. The unit of analysis is the simulated *locus*, not
the genome: each record in the truth VCF is one locus -- one mobile element and everything
that happened to it -- and each is scored once, against all of the simulated genomes at the
same time. Truth and prediction samples are paired by name (the pairing ``Compare`` is
normally driven with, ``-I $sample -J $sample``), so no genome has to be named on the
command line. Each locus then carries

  * whether the prediction recovered the locus at all, and how far off its breakpoint and
    allele length were,
  * which genomes were called as carriers versus which genomes actually carry it, and
  * whether the genotype of each paired genome agrees by allele content,

and those per-locus records are aggregated into overall metrics plus breakdowns by event
type, TE family, element size and allele frequency.

A locus is not the same thing as an event: an element inserted in one generation and
excised in a later one is one locus carrying two events. Every count here is of loci
except the event-type breakdown, where a locus is counted under each event class in its
history -- so those rows total the events, and total more than the loci.

Allele comparison, family parsing and the confusion-count arithmetic are shared with
``Compare`` so the two subcommands agree on what "the same allele" means.
'''
import json
import os
import re
import sys
from collections import OrderedDict

import pysam

from .compare_vcf import (
    calculate_metrics,
    convert_to_ploidy,
    genotype_agrees,
    parse_te_ids,
    resolve_alleles,
)
from .report import STRATUM_TITLES, write_report

# Size strata, in bp of net length change. The defaults are chosen around the shapes an LTR
# simulation produces: a solo LTR (~340 bp), a full-length element (~6 kb) and a stacked
# pair (~12 kb) each land in their own bin.
DEFAULT_SIZE_BINS = (100, 500, 1000, 5000, 10000)
# Allele-frequency strata. Rare events are the hard ones, so the low end is finer.
DEFAULT_AF_BINS = (0.1, 0.25, 0.5, 0.75)
# Above this many simulated genomes an exact-carrier-count table is more noise than signal.
MAX_GENOMES_FOR_CARRIER_TABLE = 20


class Locus:
    '''
    One simulated locus: a single record of the truth VCF, with all of its genomes.

    That is one mobile element and its whole history -- an insertion, and any excision that
    later reduced it to a solo LTR -- so a locus may carry more than one event.
    '''

    def __init__(self, chrom, pos, varID, ref, alts, info, family, superfamily,
                 gts, descs, alt_descs):
        self.chrom = chrom
        self.pos = pos
        self.id = varID
        self.ref = ref
        self.alts = alts
        self.info = info
        self.family = family
        self.superfamily = superfamily
        self.gts = gts              # truth sample -> GT tuple
        self.descs = descs          # truth sample -> resolved Allele tuple (or None)
        self.alt_descs = alt_descs  # every non-reference allele of the record
        self.match = None           # the prediction record paired with this locus
        self.displaced = None       # a call that found it, but one solo LTR downstream
        self.index = None           # position in the truth file


class VcfRow:
    '''
    One record of a file as it stands in that file, whether or not it was scored.

    The report shows truth and prediction side by side, whole, so a record that never
    became a Locus still needs a row: dropping the ones a ``--TEtype`` filter removed would
    make the file look like a smaller simulation than it is, and dropping the prediction's
    all-reference records would hide that the caller emitted them. ``scored`` is the Locus
    or PredRecord this row became, or None with ``skipped`` naming why it became nothing.
    '''

    def __init__(self, index, chrom, pos, text, scored=None, skipped=None):
        self.index = index      # position in the file, which is what the report keys on
        self.chrom = chrom
        self.pos = pos
        self.text = text
        self.scored = scored
        self.skipped = skipped


class PredRecord:
    '''One record of the prediction, kept with every sample it genotypes.'''

    def __init__(self, chrom, pos, varID, ref, alts, gts, descs, alt_descs, nonvariant):
        self.chrom = chrom
        self.pos = pos
        self.id = varID
        self.ref = ref
        self.alts = alts
        self.gts = gts
        self.descs = descs
        self.alt_descs = alt_descs
        self.nonvariant = nonvariant  # allele indices that mean "no variant here"
        self.matched = False
        self.matched_for = None       # the locus this call was paired with
        self.displaced_for = None     # the locus this call found at the wrong anchor
        self.index = None             # position in the prediction file
        self.carrier_alts = None      # --carrier_info: the ALT indices that are carrier alleles (None: not given)
        self.span_end = pos           # last position the record stands for (see tsd_span_end)


# ---- loading -----------------------------------------------------------------


def alt_descriptors(ref, alts, skip=()):
    '''
    Resolve every ALT of a record into an Allele, independent of any sample's genotype.

    ``skip`` names ALT strings that do not describe a variant, such as ``<*>``, so they
    neither pad the allele list nor offer themselves as a match.
    '''
    if not alts:
        return ()
    descs = resolve_alleles(ref, alts, tuple(range(1, len(alts) + 1)))
    if descs is None:
        return ()
    return tuple(d for d, a in zip(descs, alts) if a not in skip)


def event_type(record, info):
    '''
    The event class of a truth record: ``INS``, ``DEL``, ``EXC`` or ``OTHER``.

    ``INFO/TYPE`` is authoritative when present (every VCF ``Simulate`` writes carries it);
    a hand-built or third-party truth file without it falls back to the REF/ALT shapes.
    '''
    declared = info.get("TYPE")
    if declared:
        return str(declared)
    alt = record.alts[0] if record.alts else ""
    if len(record.ref) == 1 and len(alt) > 1:
        return "INS"
    if len(record.ref) > 1 and len(alt) == 1:
        return "DEL"
    return "OTHER"


def alt_event_types(record, info):
    '''
    The event behind each ALT allele, one value per ALT.

    ``INFO/EVENTTYPE`` is Number=A and carries exactly this. ``INFO/TYPE`` is Number=1 and
    names only the event that created the record, so an element inserted and later excised
    reads TYPE=INS and says nothing about the excision -- which is why the per-allele field
    is what the event classes are counted from.
    '''
    n_alts = len(record.alts or ())
    raw = info.get("EVENTTYPE")
    if not raw:
        return [event_type(record, info)] * n_alts
    types = [str(t) for t in _as_tuple(raw)]
    if len(types) < n_alts:
        types += [types[-1] if types else "?"] * (n_alts - len(types))
    return types[:n_alts]


# What each event class does to the length of the locus. An insertion adds sequence, an
# excision and a deletion take it away.
EVENT_SHAPE = {"INS": "gain", "EXC": "loss", "DEL": "loss"}


def allele_deltas(alt_descs):
    '''
    How much sequence each ALT adds or removes relative to the allele before it.

    Read against the previous allele rather than against REF, because that is the allele it
    was built from: on a record where an element was inserted and later excised, the
    excision allele is still far longer than REF -- it holds the anchor and the solo LTR --
    and only reads as a loss beside the full-length allele that preceded it.
    '''
    deltas = []
    previous = 0                                  # REF
    for desc in alt_descs:
        if desc.symbolic or desc.size is None:
            deltas.append(None)
            continue
        deltas.append(desc.size - previous)
        previous = desc.size
    return deltas


def event_classes(alt_descs, event_types, tol):
    '''
    Which event classes a record's alleles bear out, and which its EVENTTYPE declares
    without an allele of the matching shape.

    An event is only counted where the allele it names actually does what that event does:
    an INS allele has to add sequence and an EXC allele has to take it away. Otherwise a
    mislabelled or half-merged record would be scored under a class it never demonstrates.
    Changes within ``tol`` count as neither, since that is the slack at which two alleles
    are treated as the same length everywhere else.

    Returns (supported classes, descriptions of the declared-but-unsupported ones).
    '''
    deltas = allele_deltas(alt_descs)
    supported, unsupported = set(), []
    for index, event in enumerate(event_types):
        wanted = EVENT_SHAPE.get(event)
        delta = deltas[index] if index < len(deltas) else None
        if wanted is None or delta is None:
            # An event class with no length signature, or a symbolic allele with no length
            # to read: nothing to contradict, so take the record at its word.
            supported.add(event)
            continue
        found = "gain" if delta > tol else "loss" if delta < -tol else "flat"
        if found == wanted:
            supported.add(event)
        else:
            unsupported.append(f"{event} allele {index + 1} {found}s {abs(delta)}bp")
    return supported, unsupported


def event_history(record, info):
    '''
    What happened to this record's element over the whole simulation, e.g. ``INS`` for a
    plain insertion, ``INS,EXC`` for one that was later excised to a solo LTR, and
    ``nested INS`` for the inner element of a stacked pair.

    ``INFO/EVENTTYPE`` carries one value per ALT allele and ``INFO/MEPRESENT`` says which
    of those alleles actually hold this record's element, so only the alleles that hold it
    say anything about its history. A truth VCF straight out of ``Simulate`` writes one
    event per record, so the history is just the type; a backtracked lineage VCF merges an
    element's later history onto the same record and the history is where that shows up.
    '''
    n_alts = len(record.alts or ())
    raw_types = info.get("EVENTTYPE")
    if not raw_types:
        return event_type(record, info)
    types = [str(t) for t in _as_tuple(raw_types)]
    types += [types[-1] if types else "?"] * (n_alts - len(types))
    raw_present = info.get("MEPRESENT")
    if raw_present:
        present = [int(v) == 1 for v in _as_tuple(raw_present)]
        present += [True] * (n_alts - len(present))
    else:
        present = [True] * n_alts
    history = [t for t, p in zip(types, present) if p]
    collapsed = [t for i, t in enumerate(history) if i == 0 or t != history[i - 1]]
    label = ",".join(collapsed) or "-"
    # The element of a stacked pair that is absent from the locus's first allele is the one
    # that landed inside the other.
    if present and not present[0]:
        label = f"nested {label}"
    return label


def _as_tuple(value):
    return value if isinstance(value, (tuple, list)) else (value,)


def event_size(alt_descs):
    '''
    The size of an event, as the net length change of its largest allele, keeping the sign
    so that an excision reads negative. Returns None when no allele carries a sequence.
    '''
    sizes = [d.size for d in alt_descs if not d.symbolic and d.size is not None]
    if not sizes:
        return None
    return max(sizes, key=abs)


def load_truth_events(vcf_file, INSonly, TEtype):
    '''
    Read every record of the truth VCF as one event, keeping all of its genomes.

    Mirrors ``Compare``'s ``--TEtype``/``--INSonly`` filtering, including its refusal to run
    against a family that is not in the file: filtering to a family that does not exist
    yields an empty benchmark whose 0% detection rate says nothing about the prediction.

    Returns (events, samples, rows), where ``rows`` is every record of the file in order --
    the ones the filters dropped included -- so the report can show the file whole.
    '''
    vcf = pysam.VariantFile(vcf_file)
    samples = list(vcf.header.samples)
    if not samples:
        raise ValueError(f"Truth VCF {vcf_file} has no sample columns to evaluate")
    wanted = TEtype.casefold() if TEtype is not None else None
    seen_families = set()
    events = []
    rows = []
    for index, record in enumerate(vcf):
        varID = record.id or f"{record.chrom}:{record.pos}"
        info = dict(record.info)
        etype = event_type(record, info)
        text = str(record).rstrip("\n")
        if INSonly and etype != "INS":
            rows.append(VcfRow(index, record.chrom, record.pos, text,
                               skipped=f"not an insertion ({etype}); dropped by --INSonly"))
            continue
        family, superfamily = parse_te_ids(varID)
        seen_families.add(family)
        if wanted is not None and family.casefold() != wanted:
            rows.append(VcfRow(index, record.chrom, record.pos, text,
                               skipped=f"family {family}; dropped by --TEtype {TEtype}"))
            continue
        gts, descs = {}, {}
        for sample in samples:
            gt = record.samples[sample].get("GT")
            gt = tuple(gt) if gt is not None else (None,)
            gts[sample] = gt
            descs[sample] = resolve_alleles(record.ref, record.alts, gt)
        alt_descs = alt_descriptors(record.ref, record.alts)
        event = Locus(record.chrom, record.pos, varID, record.ref, tuple(record.alts or ()),
                      info, family, superfamily, gts, descs, alt_descs)
        event.type = etype
        event.history = event_history(record, info)
        event.alt_events = alt_event_types(record, info)
        event.size = event_size(alt_descs)
        event.index = index
        event.text = text
        events.append(event)
        rows.append(VcfRow(index, record.chrom, record.pos, text, scored=event))
    if wanted is not None and not any(f.casefold() == wanted for f in seen_families):
        raise ValueError(
            f"--TEtype '{TEtype}' matches no TE family in {vcf_file}. "
            f"Families present: {', '.join(sorted(seen_families))}. "
            "Note --TEtype matches the family (e.g. TY1-FULL), not the superfamily (Copia).")
    return events, samples, rows


# A carrier-info value marking less than this fraction of the REF value's annotated bp is a remnant of the REF element
# (an excision's solo LTR), not the element: its allele carries the element's absence.
REMNANT_FRAC = 0.5


def _annotated_bp(value):
    '''The bp a carrier-info value annotates, when it is written as miniME's ME_INFO is -- "|"-joined
    contig/start/end/family pieces -- or None.'''
    total = 0
    for piece in value.split("|"):
        parts = piece.split("/")
        try:
            total += int(parts[2]) - int(parts[1])
        except (IndexError, ValueError):
            return None
    return total


def carrier_alleles(info_text, key, n_alts):
    '''
    The ALTs an INFO field marks as carrier alleles, for --carrier_info: a set of 1-based ALT indices.

    A field with one value per allele -- REF first (Number=R, like miniME's ME_INFO) or ALTs only (Number=A) --
    marks each ALT by its own value, and "." (or empty) means that ALT is not a carrier allele -- unless the REF
    value is itself marked: then the reference holds the element, the event is its absence (a deletion, whose
    carriers are the genomes LACKING the element), and the carrier alleles are the ALTs marked "." or marking a remnant of
    it (under REMNANT_FRAC of its annotated bp: an excision's solo LTR). Any other field on a
    single-ALT record (GraffiTE's repeat_ids) marks that one ALT by being present and not ".". A record without the
    field has no carrier allele. None where a multi-ALT record's field cannot be read per allele, so the caller
    falls back to inferring it.
    '''
    value = None
    for item in info_text.split(";"):
        name, _, val = item.partition("=")
        if name == key:
            value = val if _ else ""
            break
    if value is None:
        return set()
    values = value.split(",")
    if len(values) == n_alts + 1 and n_alts > 0:
        if values[0] not in ("", "."):
            # An excision's ALT keeps a solo LTR, so it is marked too: what it holds is a remnant of the REF element.
            ref_bp = _annotated_bp(values[0])
            return {i + 1 for i, v in enumerate(values[1:])
                    if v in ("", ".") or (ref_bp and (_annotated_bp(v) or ref_bp) < REMNANT_FRAC * ref_bp)}
        per = values[1:]
    elif len(values) == n_alts:
        per = values
    elif n_alts == 1:
        per = ["" if value in ("", ".") else value]
    else:
        return None
    return {i + 1 for i, v in enumerate(per) if v not in ("", ".")}


def pred_is_carrier(record, gt):
    '''Whether a genotype of a prediction record carries a variant: by --carrier_info's alleles when given.'''
    if record.carrier_alts is not None:
        return any(a in record.carrier_alts for a in gt if a is not None)
    return is_carrier(gt, record.nonvariant)


def tsd_span_end(pos, ref, info_text):
    '''
    The last position a prediction stands for: POS, unless its REF allele starts with the TSD its INFO names.

    A caller that writes the TSD into REF (miniME: REF is the single copy of the target site, a carrier's ALT is
    TSD + element + TSD) anchors POS at the start of that copy, while a simulated insertion sits at the point of
    insertion past it -- the same event, TSD length - 1 bp apart. Such a record matches anywhere along
    POS .. POS + len(TSD) - 1, so neighbouring insertions a few bp apart are not paired across.
    '''
    for item in info_text.split(";"):
        name, _, val = item.partition("=")
        if name == "TSD" and val and val != "." and ref.upper().startswith(val.upper()):
            return pos + len(val) - 1
    return pos


def pred_offset(record, pos):
    '''Signed distance of a prediction from a position: 0 inside the positions it stands for, else from the nearer
    end, positive where the call lies to the right.'''
    if pos < record.pos:
        return record.pos - pos
    end = getattr(record, "span_end", record.pos)
    return end - pos if pos > end else 0


def load_pred_vcf(vcf_file, carrier_info=None):
    '''
    Read the prediction VCF, keeping every record that any sample calls as a variant.

    A record all of whose samples are reference (or ``<*>``, the "no variant in this
    genome" allele) describes nothing that could match a simulated event, and counting it
    among the predictions would deflate precision for free.

    Returns (records, samples, rows), where ``rows`` is every record of the file in order,
    including the all-reference ones no prediction was read from.
    '''
    vcf = pysam.VariantFile(vcf_file)
    samples = list(vcf.header.samples)
    records = []
    rows = []
    for index, record in enumerate(vcf):
        nonvariant = {0}
        if record.alts and "<*>" in record.alts:
            nonvariant.add(record.alts.index("<*>") + 1)
        gts, descs = {}, {}
        called = False
        for sample in samples:
            gt = record.samples[sample].get("GT")
            gt = tuple(gt) if gt is not None else (None,)
            gts[sample] = gt
            descs[sample] = resolve_alleles(record.ref, record.alts, gt)
            if any(a is not None and a not in nonvariant for a in gt):
                called = True
        text = str(record).rstrip("\n")
        # A sample-less VCF (a sites-only call set) still describes predictions; keep it.
        if samples and not called:
            rows.append(VcfRow(index, record.chrom, record.pos, text,
                               skipped="no sample calls this a variant"))
            continue
        alt_descs = alt_descriptors(record.ref, record.alts, skip=("<*>",))
        pred = PredRecord(record.chrom, record.pos, record.id or ".", record.ref,
                          tuple(record.alts or ()), gts, descs, alt_descs, nonvariant)
        pred.index = index
        pred.text = text
        pred.span_end = tsd_span_end(record.pos, record.ref, text.split("\t")[7])
        if carrier_info:
            pred.carrier_alts = carrier_alleles(text.split("\t")[7], carrier_info, len(record.alts or ()))
        records.append(pred)
        rows.append(VcfRow(index, record.chrom, record.pos, text, scored=pred))
    return records, samples, rows


def load_pred_bed(bed_file):
    '''
    Read a BED prediction. A BED carries no allele sequences and no sample columns, so it
    can only be scored at the locus level: detection and breakpoint offset, no genotypes.
    '''
    records = []
    rows = []
    with open(bed_file) as fin:
        for index, line in enumerate(fin):
            if line.startswith("#") or not line.strip():
                continue
            fields = line.split()
            chrom, start = fields[0], int(fields[1])
            varID = fields[3] if len(fields) > 3 else "."
            pred = PredRecord(chrom, start, varID, None, (), {}, {}, (), {0})
            pred.index = index
            pred.text = line.rstrip("\n")
            records.append(pred)
            rows.append(VcfRow(index, chrom, start, pred.text, scored=pred))
    return records, [], rows


# ---- sample pairing ----------------------------------------------------------


def read_sample_map(path):
    '''Read a two-column ``truth_sample<TAB>pred_sample`` mapping file.'''
    pairs = []
    with open(path) as fin:
        for lineno, line in enumerate(fin, 1):
            if line.startswith("#") or not line.strip():
                continue
            fields = line.split()
            if len(fields) < 2:
                raise ValueError(
                    f"{path}:{lineno}: expected 'truth_sample<TAB>pred_sample', got {line!r}")
            pairs.append((fields[0], fields[1]))
    return pairs


def pair_samples(truth_samples, pred_samples, sample_map=None):
    '''
    Decide which prediction sample stands for which simulated genome.

    An explicit ``--sample_map`` wins. Otherwise the samples are paired by name, which is
    the pairing a per-genome ``Compare`` loop is normally driven with; as a convenience a
    lone sample on each side is paired whatever it is called. Nothing is paired by column
    order: a silently wrong pairing would report confident, meaningless genotype accuracy.
    '''
    truth_set, pred_set = set(truth_samples), set(pred_samples)
    if sample_map:
        pairs = []
        for tname, pname in sample_map:
            if tname not in truth_set:
                raise ValueError(f"--sample_map names '{tname}', which is not a truth sample")
            if pname not in pred_set:
                raise ValueError(f"--sample_map names '{pname}', which is not a prediction sample")
            pairs.append((tname, pname))
        return pairs, "sample map"
    shared = [s for s in truth_samples if s in pred_set]
    if shared:
        return [(s, s) for s in shared], "matching names"
    if len(truth_samples) == 1 and len(pred_samples) == 1:
        return [(truth_samples[0], pred_samples[0])], "the only sample on each side"
    return [], "nothing"


# ---- matching ----------------------------------------------------------------


def alleles_overlap(truth_descs, pred_descs, tol):
    '''Does any non-reference allele of the event also appear in the prediction record?'''
    for t in truth_descs:
        for p in pred_descs:
            if t.symbolic or p.symbolic:
                if t.symbolic and p.symbolic and t.seq == p.seq:
                    return True
                continue
            if t.seq == p.seq or abs(t.size - p.size) <= tol:
                return True
    return False


def nearby_events(events, max_dist):
    '''Each event's simulated neighbours: the other events on its contig within ``max_dist`` of it.'''
    by_chrom = {}
    for e in events:
        by_chrom.setdefault(e.chrom, []).append(e)
    out = {}
    for chrom, group in by_chrom.items():
        group.sort(key=lambda e: e.pos)
        for i, e in enumerate(group):
            near = []
            for j in range(i - 1, -1, -1):
                if e.pos - group[j].pos > max_dist:
                    break
                near.append(group[j])
            for j in range(i + 1, len(group)):
                if group[j].pos - e.pos > max_dist:
                    break
                near.append(group[j])
            out[id(e)] = near
    return out


def _closest(size, descs):
    sizes = [d.size for d in descs if not d.symbolic and d.size is not None]
    return min((abs(size - s) for s in sizes), default=None)


# An ALT of a multi-allelic call is the event's own allele only if it is at least element-sized: within this fraction
# of the event's allele length, the same direction (insertion or deletion). Wide on purpose -- a caller writing a
# site's alleles whole folds flank variation into them (a TY4 beside a 218 bp deletion; a call 560 bp short) --
# and narrow enough that a SNP or small indel written as its own ALT is not taken for a 6 kb element.
EVENT_ALT_LENGTH_FRACTION = 0.5


def carries_event(event, neighbours, tsample, p_descs):
    '''
    Whether a genome's allele in a multi-allelic call is this event's.

    A record with several ALTs holds several alleles at one site, and a genome with any of them used to count as
    carrying the matched event. Two ways that is wrong: the genome's ALT is a different, much smaller variant at
    the site (a SNP or small indel written as its own allele), or it is a neighbouring simulated insertion the
    genome carries -- where two target sites overlap, the neighbour interrupts this event's site, and the caller
    rightly writes that as another ALT of this record. So an ALT is the event's when it is element-sized
    (EVENT_ALT_LENGTH_FRACTION) and no event within --max_dist that the genome really carries is closer to it in
    length.
    '''
    if not p_descs:
        return True
    sizes = [d.size for d in p_descs if d.idx != 0 and d.size is not None and not d.symbolic]
    if not sizes or any(d.symbolic for d in p_descs if d.idx != 0) or \
            any(t.symbolic or t.size is None for t in event.alt_descs) or not event.alt_descs:
        return True                 # nothing comparable by length: every non-reference allele counts, as before
    for size in sizes:
        sized = [t.size for t in event.alt_descs
                 if t.size and (size > 0) == (t.size > 0)
                 and abs(size - t.size) <= EVENT_ALT_LENGTH_FRACTION * abs(t.size)]
        if not sized:
            continue
        mine = min(abs(size - s) for s in sized)
        if not any(is_carrier(n.gts.get(tsample, (None,)), {0}) and
                   (_closest(size, n.alt_descs) is not None and _closest(size, n.alt_descs) < mine)
                   for n in neighbours):
            return True
    return False


def best_length_error(truth_descs, pred_descs):
    '''
    The smallest length disagreement between any pair of truth and predicted alleles, i.e.
    how far the prediction is from the simulated element at its closest. None when either
    side carries no comparable sequence.
    '''
    errors = [p.size - t.size
              for t in truth_descs if not t.symbolic and t.size is not None
              for p in pred_descs if not p.symbolic and p.size is not None]
    if not errors:
        return None
    return min(errors, key=abs)


def _haplotype(fasta, chrom, start, end, pos, ref, allele):
    '''
    The sequence a genome carrying ``allele`` has over reference [start, end): the record's
    REF span replaced by the allele. Upper-cased, because a simulated TSD copies soft-masked
    (lower-case) reference bases while callers write them in upper case -- the same bases.
    '''
    return (fasta.fetch(chrom, start, pos - 1) + allele
            + fasta.fetch(chrom, pos - 1 + len(ref), end)).upper()


def haplotype_identical(event, match, pairs, fasta):
    '''
    Is the called allele, written into the reference, base-for-base the simulated one?

    Evaluate otherwise compares alleles by net length only (``--gt_len_tol``), so a call of
    the right size with the wrong sequence -- a TSD copy dropped, a junction misplaced, a
    neighbour's allele -- passes. Comparing the haplotypes rather than the ALT strings is
    what makes this fair across representations: an insertion with a TSD can be written at
    any anchor along the duplicated bases (the simulation writes it after the TSD, a
    left-aligning caller before it), and every such record yields the same haplotype.

    Where both files call a genome a carrier, each such genome's non-reference alleles must
    give the same haplotypes. Where no genome is a carrier in both, the locus passes if any
    simulated allele and any called allele give the same haplotype. None when nothing is
    comparable: a symbolic allele, or a contig the reference does not have.
    '''
    if match is None or event.chrom != match.chrom or event.chrom not in fasta.references:
        return None
    length = fasta.get_reference_length(event.chrom)
    start = max(0, min(event.pos, match.pos) - 2)
    end = min(length, max(event.pos - 1 + len(event.ref), match.pos - 1 + len(match.ref)) + 1)
    cache = {}

    def hap(record, allele):
        key = (id(record), allele)
        if key not in cache:
            cache[key] = _haplotype(fasta, record.chrom, start, end, record.pos, record.ref, allele)
        return cache[key]

    def non_ref(descs, nonvariant=()):
        if not descs:
            return []
        return [d for d in descs if d.idx != 0 and d.idx not in nonvariant]

    compared = False
    for tsample, psample in pairs:
        t_descs = non_ref(event.descs.get(tsample))
        p_descs = non_ref(match.descs.get(psample), match.nonvariant)
        if not t_descs or not p_descs:
            continue
        if any(d.symbolic for d in t_descs + p_descs):
            return None
        compared = True
        if sorted(hap(event, d.seq) for d in t_descs) != sorted(hap(match, d.seq) for d in p_descs):
            return False
    if compared:
        return True
    t_alts = [d for d in event.alt_descs if not d.symbolic]
    p_alts = [d for d in match.alt_descs if not d.symbolic and d.seq != "<*>"]
    if not t_alts or not p_alts:
        return None
    return any(hap(event, t.seq) == hap(match, p.seq) for t in t_alts for p in p_alts)


def match_events(events, records, max_dist, tol):
    '''
    Pair each simulated event with at most one prediction record, and vice versa.

    Every candidate pair within ``max_dist`` is scored, and the pairs are then taken in
    order of preference: allele agreement first, breakpoint proximity second. The one-to-one
    constraint is what keeps a stacked pair of elements honest -- two events at one locus
    need two prediction records to both count as recovered, rather than both claiming the
    single record that happens to sit there.
    '''
    by_chrom = {}
    for record in records:
        by_chrom.setdefault(record.chrom, []).append(record)
    candidates = []
    for i, event in enumerate(events):
        for j, record in enumerate(by_chrom.get(event.chrom, [])):
            distance = abs(pred_offset(record, event.pos))
            if distance > max_dist:
                continue
            agree = alleles_overlap(event.alt_descs, record.alt_descs, tol)
            candidates.append((not agree, distance, i, event.chrom, j))
    candidates.sort()
    taken_events, taken_records = set(), set()
    for _disagree, _distance, i, chrom, j in candidates:
        if i in taken_events or (chrom, j) in taken_records:
            continue
        taken_events.add(i)
        taken_records.add((chrom, j))
        events[i].match = by_chrom[chrom][j]
        by_chrom[chrom][j].matched = True
        # The pairing is read from both ends by the side-by-side viewer, which has to get
        # from a prediction back to the locus that claimed it as well as the other way.
        by_chrom[chrom][j].matched_for = events[i]
    match_displaced(events, by_chrom, taken_events, taken_records, max_dist, tol)


def ltr_shift(event):
    '''
    How far downstream of an excision's POS a call would sit if it were anchored past the
    solo LTR instead of at its start, or None if this locus has no excision allele.

    ``INFO/LTRLEN`` gives the length directly; failing that the excision allele is the
    anchor base plus the solo LTR, so its own length says the same thing.
    '''
    if "EXC" not in getattr(event, "alt_events", ()):
        return None
    declared = event.info.get("LTRLEN")
    if declared:
        try:
            return int(_as_tuple(declared)[0])
        except (TypeError, ValueError):
            pass
    index = list(event.alt_events).index("EXC")
    if index < len(event.alts):
        return max(0, len(event.alts[index]) - 1)
    return None


def match_displaced(events, by_chrom, taken_events, taken_records, max_dist, tol):
    '''
    Pair a leftover excision with a call sitting one solo LTR downstream of it.

    A caller that anchors an SV where the two alleles start to differ puts an excision at
    the far end of the surviving LTR, since the LTR itself is identical either side and
    collapses. The simulation writes POS at the start of the element, which is where the
    solo LTR starts, so the two conventions differ by exactly the LTR length -- 338bp on
    Ty1, well past any sane --max_dist. Left alone that scores one excision twice over, as
    a miss and as a false positive.

    Such a pair is recorded as ``displaced``, never as ``match``: the locus was found, so it
    is not a miss, but the call is in the wrong place, so it is not a correct prediction
    either. It earns detection credit and nothing else -- no allele or breakpoint
    statistics, no carrier or genotype credit, and the record stays among the unmatched
    predictions, so precision is unaffected.
    '''
    for i, event in enumerate(events):
        if i in taken_events:
            continue
        shift = ltr_shift(event)
        if not shift:
            continue
        expected = event.pos + shift
        best = None
        for j, record in enumerate(by_chrom.get(event.chrom, [])):
            if (event.chrom, j) in taken_records:
                continue
            distance = abs(pred_offset(record, expected))
            if distance > max_dist:
                continue
            agree = alleles_overlap(event.alt_descs, record.alt_descs, tol)
            key = (not agree, distance, j)
            if best is None or key < best[0]:
                best = (key, j)
        if best is None:
            continue
        j = best[1]
        taken_events.add(i)
        taken_records.add((event.chrom, j))
        event.displaced = by_chrom[event.chrom][j]
        by_chrom[event.chrom][j].displaced_for = event


# ---- per-event scoring -------------------------------------------------------


def is_carrier(gt, nonvariant):
    '''Does this genotype hold the variant at all? Missing calls count as non-carriers.'''
    return any(a is not None and a not in nonvariant for a in gt)


def match_ploidy(t_gt, t_descs, p_gt, p_descs):
    '''
    Replicate the shorter genotype so a haploid truth can be compared against a diploid
    call (``1`` versus ``1/1``), the same widening ``Compare`` applies. Genotypes whose
    lengths are not multiples of one another are left alone and simply fail to agree.
    '''
    lt, lp = len(t_gt), len(p_gt)
    if lt == lp or lt == 0 or lp == 0:
        return t_gt, t_descs, p_gt, p_descs
    if lp > lt and lp % lt == 0:
        k = lp // lt
        return t_gt * k, (t_descs * k if t_descs else t_descs), p_gt, p_descs
    if lt > lp and lt % lp == 0:
        k = lt // lp
        return t_gt, t_descs, p_gt * k, (p_descs * k if p_descs else p_descs)
    return t_gt, t_descs, p_gt, p_descs


def score_event(event, pairs, tol, neighbours=()):
    '''Turn one simulated event and its matched prediction into a scored record.'''
    match = event.match
    truth_carriers, pred_carriers = [], []
    tp = fp = fn = 0
    compared = concordant = 0
    by_class = OrderedDict((cls, _blank_counts()) for cls in ALLELE_CLASSES)
    multi = match is not None and \
        len([a for a in (match.alts or ()) if a != "<*>"]) > 1
    untyped = []
    for tsample, psample in pairs:
        t_gt = event.gts.get(tsample, (None,))
        t_carrier = is_carrier(t_gt, {0})
        if t_carrier:
            truth_carriers.append(tsample)
        if match is None:
            p_carrier = False
            typed = False
        else:
            p_gt = match.gts.get(psample, (None,))
            p_carrier = pred_is_carrier(match, p_gt)
            if p_carrier and multi and match.carrier_alts is None and not carries_event(event, neighbours, tsample,
                                                         match.descs.get(psample)):
                p_carrier = False
            typed = any(a is not None for a in p_gt)
            if p_carrier:
                pred_carriers.append(psample)
            # A genotype is only comparable where a prediction record describes the locus.
            wt_gt, wt_descs, wp_gt, wp_descs = match_ploidy(
                t_gt, event.descs.get(tsample), p_gt, match.descs.get(psample))
            same = genotype_agrees(wt_descs, wp_descs, tol)
            if same is None:
                same = sorted(a for a in wt_gt if a is not None) == \
                       sorted(a for a in wp_gt if a is not None)
            compared += 1
            concordant += bool(same)
        tp += t_carrier and p_carrier
        fp += (not t_carrier) and p_carrier
        fn += t_carrier and (not p_carrier)
        # The same genome scored by allele class. A genome no correctly placed call
        # genotypes -- the locus was missed or displaced, or the call says "." -- has no
        # allele here at all: an FN of its own class and an FP of neither.
        t_cls = "carrier" if t_carrier else "noncarrier"
        counts = by_class[t_cls]
        counts["total"] += 1
        if not typed:
            counts["fn"] += 1
            untyped.append(tsample)
            continue
        counts["genotyped"] += 1
        p_cls = "carrier" if p_carrier else "noncarrier"
        if p_cls == t_cls:
            counts["tp"] += 1
        else:
            counts["fn"] += 1
            by_class[p_cls]["fp"] += 1

    # Allele frequency is a property of the simulation, so it is read off every genome the
    # truth VCF holds -- not only the ones that could be paired with a prediction.
    alleles = [a for gt in event.gts.values() for a in gt if a is not None]
    an = len(alleles)
    ac = sum(1 for a in alleles if a != 0)
    n_truth_carriers = sum(1 for gt in event.gts.values() if is_carrier(gt, {0}))

    scored = OrderedDict()
    scored["chrom"] = event.chrom
    scored["pos"] = event.pos
    scored["id"] = event.id
    scored["type"] = event.type
    scored["history"] = event.history
    supported, unsupported = event_classes(event.alt_descs, event.alt_events, tol)
    scored["event_classes"] = sorted(supported)
    scored["unsupported_events"] = unsupported
    scored["family"] = event.family
    scored["superfamily"] = event.superfamily
    scored["size_bp"] = event.size
    scored["allele_frequency"] = round(ac / an, 4) if an else None
    scored["allele_count"] = ac
    scored["allele_number"] = an
    scored["n_carrier_genomes"] = n_truth_carriers
    # "detected" stays the strict sense -- a correctly placed call. A displaced one
    # found the locus but at the wrong anchor, so it counts toward recovery and nothing
    # else: no allele, breakpoint, carrier or genotype credit, and its record stays among
    # the unmatched predictions.
    scored["detected"] = match is not None
    scored["displaced"] = None if event.displaced is None else OrderedDict([
        ("id", event.displaced.id),
        ("pos", event.displaced.pos),
        ("pos_offset", pred_offset(event.displaced, event.pos)),
        ("allele_bp", event_size(event.displaced.alt_descs)),
        ("expected_at", event.pos),
    ])
    scored["recovered"] = match is not None or event.displaced is not None
    if match is None:
        scored["match"] = None
    else:
        length_error = best_length_error(event.alt_descs, match.alt_descs)
        scored["match"] = OrderedDict([
            ("id", match.id),
            ("pos", match.pos),
            ("pos_offset", pred_offset(match, event.pos)),
            ("allele_bp", event_size(match.alt_descs)),
            ("length_error", length_error),
            ("allele_match", alleles_overlap(event.alt_descs, match.alt_descs, tol)
                             if match.alt_descs else None),
        ])
    scored["carriers"] = OrderedDict([
        ("truth", truth_carriers),
        ("predicted", pred_carriers),
        ("tp", tp), ("fp", fp), ("fn", fn),
    ])
    scored["genotypes"] = OrderedDict([("compared", compared), ("concordant", concordant)])
    scored["alleles"] = OrderedDict(list(by_class.items()) + [("untyped", untyped)])
    return scored


# ---- aggregation -------------------------------------------------------------


def _spread(values):
    '''
    Mean and standard deviation of a signed error.

    Signed rather than absolute: a call that lands short of the simulated value reads
    negative and one that overshoots positive, so a systematic bias shows up as a mean away
    from zero instead of being folded in with random jitter -- a caller consistently 3bp
    downstream and one scattering +-3bp at random are the same number in absolute terms and
    quite different problems. The SD is the population SD: these are all the events that
    were recovered, not a sample drawn from more of them.
    '''
    if not values:
        return OrderedDict([("n", 0), ("mean", None), ("sd", None)])
    mean = sum(values) / len(values)
    variance = sum((value - mean) ** 2 for value in values) / len(values)
    return OrderedDict([
        ("n", len(values)),
        ("mean", round(mean, 2)),
        ("sd", round(variance ** 0.5, 2)),
    ])


# The two classes each paired genome's allele at a locus is scored as, by what the truth says
# it holds there. For a deletion the carrier is the genome that LACKS the reference's element.
ALLELE_CLASSES = ("carrier", "noncarrier")


def _blank_counts():
    return OrderedDict([("total", 0), ("genotyped", 0), ("tp", 0), ("fp", 0), ("fn", 0)])


def with_rates(counts):
    '''
    Add recall, precision and F1 to a block of confusion counts. Precision and F1 are left
    out (None) where FP is: a breakdown row has no locus FPs, since a call on no simulated
    locus belongs to no stratum.
    '''
    out = OrderedDict(counts)
    tp, fp, fn = out["tp"], out["fp"], out["fn"]
    out["recall"] = round(tp / (tp + fn), 4) if tp + fn else None
    if fp is None:
        out["precision"] = out["f1"] = None
    else:
        f1, precision, _ = calculate_metrics(tp, fp, fn)
        out["precision"] = round(precision, 4) if tp + fp else None
        out["f1"] = round(f1, 4) if tp + fp and tp + fn else None
    return out


def counts_block(n_loci, n_recalled, allele_counts, locus_fp=None):
    '''
    TP, FP, FN, recall, precision and F1 by locus and by allele.

    By locus: a simulated locus a call was placed on (displaced or not) is a TP, one with no
    call an FN; an FP is a call on no simulated locus, known only for the whole run.

    By allele, over every paired genome at every simulated locus: each class in turn is the
    positive one. TP is an allele genotyped as its own class; FN is one that was not --
    genotyped as the other class, or never genotyped; FP is an allele of the other class
    genotyped as this one. ``genotyped`` counts the alleles that got a genotype at all, so
    recall is at most genotyped / total. "all" sums the two classes.
    '''
    loci = OrderedDict([("total", n_loci), ("genotyped", n_recalled), ("tp", n_recalled),
                        ("fp", locus_fp), ("fn", n_loci - n_recalled)])
    block = OrderedDict([("loci", with_rates(loci))])
    total = _blank_counts()
    for cls in ALLELE_CLASSES:
        for key in total:
            total[key] += allele_counts[cls][key]
        block[cls] = with_rates(allele_counts[cls])
    block["all"] = with_rates(total)
    return block


def aggregate(scored_events):
    '''Roll a list of scored events up into one block of metrics.'''
    n = len(scored_events)
    detected = [e for e in scored_events if e["detected"]]
    displaced = [e for e in scored_events if e.get("displaced")]
    recovered = len(detected) + len(displaced)
    allele_ok = sum(1 for e in detected if e["match"]["allele_match"])
    allele_comparable = sum(1 for e in detected if e["match"]["allele_match"] is not None)
    seq_checked = [e for e in detected if "haplotype_identical" in e["match"]]
    seq_exact = sum(1 for e in seq_checked if e["match"]["haplotype_identical"])
    seq_comparable = sum(1 for e in seq_checked if e["match"]["haplotype_identical"] is not None)
    tp = sum(e["carriers"]["tp"] for e in scored_events)
    fp = sum(e["carriers"]["fp"] for e in scored_events)
    fn = sum(e["carriers"]["fn"] for e in scored_events)
    compared = sum(e["genotypes"]["compared"] for e in scored_events)
    concordant = sum(e["genotypes"]["concordant"] for e in scored_events)
    f1, precision, recall = calculate_metrics(tp, fp, fn)
    # Signed, both of them: pos_offset is the call's position minus the simulated one and
    # length_error the called allele's length minus the simulated one, so short reads
    # negative and long reads positive.
    allele_counts = {cls: _blank_counts() for cls in ALLELE_CLASSES}
    for e in scored_events:
        for cls in ALLELE_CLASSES:
            for key, value in e["alleles"][cls].items():
                allele_counts[cls][key] += value
    offsets = [e["match"]["pos_offset"] for e in detected]
    errors = [e["match"]["length_error"] for e in detected
              if e["match"]["length_error"] is not None]
    return OrderedDict([
        ("n_loci", n),
        ("n_recovered", recovered),
        ("recall", round(recovered / n, 4) if n else None),
        ("recovery_rate", round(recovered / n, 4) if n else None),   # = recall; the old name
        ("n_detected", len(detected)),
        ("n_displaced", len(displaced)),
        ("n_allele_concordant", allele_ok),
        ("allele_concordance_rate",
         round(allele_ok / allele_comparable, 4) if allele_comparable else None),
    ] + ([
        ("n_haplotype_identical", seq_exact),
        ("n_haplotype_comparable", seq_comparable),
        ("haplotype_identity_rate", round(seq_exact / seq_comparable, 4) if seq_comparable else None),
    ] if seq_checked else []) + [
        ("breakpoint_offset_bp", _spread(offsets)),
        ("allele_length_error_bp", _spread(errors)),
        ("carriers", OrderedDict([
            ("tp", tp), ("fp", fp), ("fn", fn),
            ("recall", round(recall, 4)), ("precision", round(precision, 4)),
            ("f1", round(f1, 4)),
        ])),
        ("genotypes", OrderedDict([
            ("compared", compared),
            ("concordant", concordant),
            ("concordance", round(concordant / compared, 4) if compared else None),
        ])),
        ("counts", counts_block(n, recovered, allele_counts)),
    ])


def add_unmatched(counts, unmatched_records):
    '''
    Charge the calls that matched no simulated locus to the run's counts: each is a locus
    FP, and each genome it calls a carrier a carrier FP. A displaced call is left out of
    both -- its locus is already counted as recalled, and its genomes as not genotyped.
    '''
    spurious = [r for r in unmatched_records if r.displaced_for is None]
    carriers = sum(1 for r in spurious for gt in r.gts.values() if pred_is_carrier(r, gt))
    counts["loci"]["fp"] = len(spurious)
    counts["loci"].update(with_rates(counts["loci"]))
    counts["carrier"]["fp"] += carriers
    counts["carrier"].update(with_rates(counts["carrier"]))
    counts["all"]["fp"] += carriers
    counts["all"].update(with_rates(counts["all"]))


def load_nonmobilizing(vcf_file):
    '''The non-mobilizing SVs Simulate wrote (TErandom --nNM), one dict per record: its span, kind
    and carrier count.'''
    sites = []
    vcf = pysam.VariantFile(vcf_file)
    for record in vcf:
        carriers = sum(1 for s in record.samples.values()
                       if any(a not in (None, 0) for a in (s.get("GT") or ())))
        sites.append(OrderedDict([
            ("chrom", record.chrom), ("pos", record.pos), ("end", record.pos + len(record.ref) - 1),
            ("id", record.id), ("kind", record.info.get("NMKIND", ".")),
            ("host", record.info.get("NMHOST", ".")), ("n_carriers", carriers)]))
    return sites


BG_SV_SIZE_BINS = ((500, "50-499 bp"), (2000, "500-1,999 bp"), (None, "2-10 kb"))


def load_background_sv(vcf_file):
    '''The background SVs Simulate wrote (--bg-sv-rate), in the shape load_nonmobilizing gives: a
    duplication spans the stretch it copies (INFO/DUPSTART to POS), a deletion what it removes. The
    kind is the type and a size class.'''
    sites = []
    vcf = pysam.VariantFile(vcf_file)
    for record in vcf:
        carriers = sum(1 for s in record.samples.values()
                       if any(a not in (None, 0) for a in (s.get("GT") or ())))
        typ = record.info.get("TYPE", ".")
        size = abs(int(record.info.get("SVLEN", 0)))
        label = next(lab for edge, lab in BG_SV_SIZE_BINS if edge is None or size < edge)
        start = int(record.info["DUPSTART"]) if typ == "DUP" and "DUPSTART" in record.info else record.pos
        end = record.pos if typ == "DUP" else record.pos + len(record.ref) - 1
        sites.append(OrderedDict([
            ("chrom", record.chrom), ("pos", start), ("end", end), ("id", record.id),
            ("kind", f"{typ} {label}"), ("size", size), ("n_carriers", carriers)]))
    return sites


def score_nonmobilizing(sites, unmatched_records, max_dist):
    '''
    The calls a caller made at non-mobilizing SVs. A call that matched no simulated locus and whose
    span comes within max_dist of one is charged to it. Those calls are already locus FPs (they
    are unmatched predictions); this says how many there are and at which kinds of SV, since a
    TE caller that reports one is mistaking a variant no TE made for a TE event. A displaced call
    is left out, as it is from the FPs.
    '''
    by_chrom = {}
    for k, site in enumerate(sites):
        by_chrom.setdefault(site["chrom"], []).append((site["pos"] - max_dist, site["end"] + max_dist, k))
    calls_at = {}
    n_calls = 0
    for r in unmatched_records:
        if r.displaced_for is not None:
            continue
        end = max(r.span_end or r.pos, r.pos + len(r.ref or "N") - 1)
        hits = [k for lo, hi, k in by_chrom.get(r.chrom, ()) if r.pos <= hi and end >= lo]
        if hits:
            n_calls += 1
        for k in hits:
            calls_at.setdefault(k, []).append(f"{r.chrom}:{r.pos}:{r.id}")
    kinds = OrderedDict()
    for k, site in enumerate(sites):
        row = kinds.setdefault(site["kind"], OrderedDict([("sites", 0), ("sites_called", 0), ("calls", 0)]))
        row["sites"] += 1
        row["sites_called"] += k in calls_at
        row["calls"] += len(calls_at.get(k, ()))
    return OrderedDict([
        ("sites", len(sites)),
        ("sites_called", len(calls_at)),
        ("calls", n_calls),
        ("by_kind", kinds),
        ("called", [OrderedDict(list(sites[k].items()) + [("calls", calls_at[k])]) for k in sorted(calls_at)]),
    ])


def stratify(scored_events, key):
    '''Group events by ``key(event)`` and aggregate each group. Sorted by group size.'''
    return stratify_multi(scored_events, lambda event: (key(event),))


def stratify_multi(scored_events, keys):
    '''
    Group events by the set of labels ``keys(event)`` returns, so an event counts once in
    each label it bears. The rows then do not sum to the number of events -- an element
    that was inserted and later excised is one event under both INS and EXC.
    '''
    groups = OrderedDict()
    for event in scored_events:
        for label in keys(event) or ("?",):
            groups.setdefault(label, []).append(event)
    ordered = sorted(groups.items(), key=lambda kv: (-len(kv[1]), str(kv[0])))
    return OrderedDict((str(label), aggregate(members)) for label, members in ordered)


def size_bin(size, edges):
    '''Label the size stratum an event falls in, by the magnitude of its length change.'''
    if size is None:
        return "unknown"
    magnitude = abs(size)
    labels = _size_bin_labels(edges)
    for edge, label in zip(edges, labels):
        if magnitude < edge:
            return label
    return labels[-1]


def _size_bin_labels(edges):
    labels = [f"<{_bp(edges[0])}"]
    for low, high in zip(edges, edges[1:]):
        labels.append(f"{_bp(low)}-{_bp(high)}")
    labels.append(f">={_bp(edges[-1])}")
    return labels


def _bp(value):
    if value >= 1000 and value % 1000 == 0:
        return f"{value // 1000}kb"
    return f"{value}bp"


def af_bin(af, edges):
    '''Label the allele-frequency stratum an event falls in.'''
    if af is None:
        return "unknown"
    previous = 0.0
    for edge in edges:
        if af <= edge:
            return f"{previous:g}-{edge:g}"
        previous = edge
    return f"{previous:g}-1"


# ---- reporting ---------------------------------------------------------------


def _pct(value):
    return "n/a" if value is None else f"{value * 100:.2f}%"


def _num(value):
    return "-" if value is None else f"{value:g}" if isinstance(value, float) else str(value)


def _signed(value):
    '''
    Render a signed error with its sign always shown, so the direction of a bias reads at a
    glance. "n/a" rather than "-" for a missing value, which would look like a minus sign.
    '''
    return "n/a" if value is None else f"{value:+g}"


def _table(headers, rows):
    '''Render a fixed-width table. Numbers right-aligned, the first column left-aligned.'''
    if not rows:
        return ["  (no events)"]
    widths = [max(len(str(h)), max(len(str(r[i])) for r in rows))
              for i, h in enumerate(headers)]
    def render(cells):
        out = [str(cells[0]).ljust(widths[0])]
        out += [str(c).rjust(w) for c, w in zip(cells[1:], widths[1:])]
        return "  " + "  ".join(out)
    lines = [render(headers), "  " + "  ".join("-" * w for w in widths)]
    lines += [render(r) for r in rows]
    return lines


def _stratum_rows(strata):
    rows = []
    for label, stats in strata.items():
        rows.append([
            label,
            stats["n_loci"],
            stats["n_recovered"],
            stats["n_displaced"],
            _pct(stats["recall"]),
            _pct(stats["allele_concordance_rate"]),
            _signed(stats["breakpoint_offset_bp"]["mean"]),
            _signed(stats["allele_length_error_bp"]["mean"]),
            _num(stats["counts"]["carrier"]["f1"]),
            _num(stats["counts"]["noncarrier"]["f1"]),
            _num(stats["counts"]["all"]["f1"]),
            _pct(stats["genotypes"]["concordance"]),
        ])
    return rows


COUNT_HEADERS = ["level", "total", "called", "TP", "FP", "FN", "recall", "precision", "F1"]


def count_rows(counts, unit="haplotype"):
    '''Rows of the locus and allele table: loci, then carrier, noncarrier and all alleles.'''
    rows = []
    for key, label in (("loci", "loci"), ("carrier", f"carrier {unit}s"),
                       ("noncarrier", f"noncarrier {unit}s"), ("all", f"all {unit}s")):
        c = counts[key]
        rows.append([label, c["total"], c["genotyped"], c["tp"], _num(c["fp"]), c["fn"],
                     _num(c["recall"]), _num(c["precision"]), _num(c["f1"])])
    return rows


def stratum_headers(unit="loci"):
    '''
    Column headings for a breakdown table. The first count is labelled by what the rows
    actually hold: every table counts loci except the event-type one, where a locus is
    counted under each event class in its history and the column therefore counts events.
    '''
    return ["stratum", unit, "found", "displ", "recall",
            "allele ok", "mean off", "mean len err", "carrier F1", "noncarr F1", "all F1",
            "genotype"]


def format_report(summary, meta, title="tevarsim Evaluate", preamble=()):
    '''The printed summary: the headline metrics, then one table per breakdown.'''
    lines = []
    lines.append("")
    lines.append(title)
    lines.extend(preamble)
    lines.append(f"  truth : {meta['truth']}")
    lines.append(f"          {meta['n_loci']} simulated loci over "
                 f"{meta['n_truth_samples']} genomes")
    lines.append(f"  pred  : {meta['pred']}")
    lines.append(f"          {meta['n_pred_records']} records over "
                 f"{meta['n_pred_samples']} samples")
    lines.append(f"  paired: {len(meta['sample_pairs'])} genome(s) by {meta['pairing']}")
    if meta["sample_pairs"]:
        shown = ", ".join(f"{t}={p}" if t != p else t for t, p in meta["sample_pairs"][:6])
        if len(meta["sample_pairs"]) > 6:
            shown += f", ... (+{len(meta['sample_pairs']) - 6})"
        lines.append(f"          {shown}")
    else:
        lines.append("          no genomes could be paired -- carrier and genotype "
                     "statistics are unavailable.")
        lines.append("          Pass --sample_map truth_sample<TAB>pred_sample to pair them.")
    if meta.get("carrier_info"):
        lines.append(f"  carriers: the ALTs INFO/{meta['carrier_info']} marks as carrier alleles")
    if meta["filters"]:
        lines.append(f"  filter: {meta['filters']}")
        lines.append("          the truth is filtered but the prediction is not, so "
                     "unmatched predictions include events the filter removed.")

    overall = summary["overall"]
    lines.append("")
    lines.append(f"Locus detection ({overall['n_loci']} simulated loci)")
    lines.append(f"  recall               : {overall['n_recovered']} / {overall['n_loci']}"
                 f"  ({_pct(overall['recall'])})")
    if overall["n_displaced"]:
        lines.append(f"    correctly placed   : {overall['n_detected']}")
        lines.append(f"    displaced          : {overall['n_displaced']}"
                     "   (call sits past the solo LTR, not at its start -- found, but not "
                     "a correct call)")
    lines.append(f"  allele concordant    : {overall['n_allele_concordant']} of the "
                 f"{overall['n_detected']} correctly placed  "
                 f"({_pct(overall['allele_concordance_rate'])})")
    if "n_haplotype_identical" in overall:
        lines.append(f"  haplotype identical  : {overall['n_haplotype_identical']} of the "
                     f"{overall['n_haplotype_comparable']} comparable  "
                     f"({_pct(overall['haplotype_identity_rate'])})   (the called allele, written "
                     "into the reference, is base-for-base the simulated one)")
    displaced_note = (f"   ({summary['predictions']['displaced']} of them displaced calls)"
                      if summary["predictions"].get("displaced") else "")
    lines.append(f"  unmatched predictions: {summary['predictions']['unmatched']} of "
                 f"{summary['predictions']['total']}{displaced_note}")
    lines.append(f"  locus precision      : {_pct(summary['predictions']['precision'])}"
                 "   (prediction records correctly placed on a simulated locus)")

    offsets = overall["breakpoint_offset_bp"]
    errors = overall["allele_length_error_bp"]
    lines.append("")
    lines.append("Accuracy of the recalled loci")
    lines.append(f"  breakpoint offset   : mean {_signed(offsets['mean'])} bp, "
                 f"SD {_num(offsets['sd'])} bp   (n={offsets['n']})")
    lines.append(f"  allele length error : mean {_signed(errors['mean'])} bp, "
                 f"SD {_num(errors['sd'])} bp   (n={errors['n']})")
    lines.append("  (- is short of the simulated value, + is past it)")

    genotypes = overall["genotypes"]
    lines.append("")
    if meta["sample_pairs"]:
        unit = meta["genotype_label"]
        lines.append(f"By locus and by {unit}, over {len(meta['sample_pairs'])} paired genomes")
        lines.extend(_table(COUNT_HEADERS, count_rows(overall["counts"], unit)))
        lines.append(f"  a {unit} with no genotype (missed or displaced locus, or \".\") is an FN of "
                     "its own class and an FP of neither;")
        lines.append("  the carriers of a call on no simulated locus are carrier FPs")
        lines.append(f"  {unit} concordance: "
                     f"{genotypes['concordant']} / {genotypes['compared']} "
                     f"({_pct(genotypes['concordance'])})")
    else:
        lines.append("Allele accuracy: not available (no genomes paired)")

    unsupported = summary.get("unsupported_events") or []
    if unsupported:
        lines.append("")
        lines.append(f"Loci declaring an event no allele bears out ({len(unsupported)})")
        lines.append("  not counted under that type; the allele does not change length the "
                     "way the event would")
        for item in unsupported[:5]:
            lines.append(f"  {item['chrom']}:{item['pos']}  {'; '.join(item['declared'])}")
        if len(unsupported) > 5:
            lines.append(f"  ... and {len(unsupported) - 5} more")

    bg = summary.get("background_sv")
    if bg is not None:
        lines.append("")
        lines.append(f"Background SVs ({bg['sites']} in {meta.get('background_sv')})")
        lines.append("  deletions and tandem duplications that are no TE event, placed at random; a TE caller should")
        lines.append("  call none. A call on none of the simulated loci that comes within --max_dist of one is charged")
        lines.append("  to it, and is among the unmatched predictions above")
        lines.append(f"  sites called         : {bg['sites_called']} / {bg['sites']}"
                     f"  ({_pct(bg['sites_called'] / bg['sites'] if bg['sites'] else None)})")
        lines.append(f"  calls at them        : {bg['calls']} of the {summary['predictions']['unmatched']} unmatched predictions")
        lines.extend(_table(["type, size", "sites", "called", "calls"],
                            [[kind, row["sites"], row["sites_called"], row["calls"]]
                             for kind, row in bg["by_kind"].items()]))

    nm = summary.get("nonmobilizing")
    if nm is not None:
        lines.append("")
        lines.append(f"Non-mobilizing SVs ({nm['sites']} in {meta.get('nonmobilizing')})")
        lines.append("  deletions that cut into or swallow a reference element by no TE mechanism; a TE caller should")
        lines.append("  call none. A call on none of the simulated loci that comes within --max_dist of one is charged")
        lines.append("  to it, and is among the unmatched predictions above")
        lines.append(f"  sites called         : {nm['sites_called']} / {nm['sites']}"
                     f"  ({_pct(nm['sites_called'] / nm['sites'] if nm['sites'] else None)})")
        lines.append(f"  calls at them        : {nm['calls']} of the {summary['predictions']['unmatched']} unmatched predictions")
        lines.extend(_table(["kind", "sites", "called", "calls"],
                            [[kind, row["sites"], row["sites_called"], row["calls"]]
                             for kind, row in nm["by_kind"].items()]))

    for title, key in STRATUM_TITLES:
        strata = summary.get(key)
        if not strata:
            continue
        lines.append("")
        lines.append(title)
        if key == "by_event_type":
            lines.append("  a locus counts under each event class in its history, so "
                         "these rows count events, not loci")
        unit = "events" if key == "by_event_type" else "loci"
        lines.extend(_table(stratum_headers(unit), _stratum_rows(strata)))
    lines.append("")
    return "\n".join(lines)


# ---- driver ------------------------------------------------------------------


class Evaluator:
    def __init__(self, args):
        self.truth_file = args.truth
        self.pred_file = args.pred
        self.predType = getattr(args, "predType", "VCF")
        self.outprefix = args.outprefix
        self.TEtype = getattr(args, "TEtype", None)
        self.INSonly = getattr(args, "INSonly", False)
        self.max_dist = getattr(args, "max_dist", 100)
        self.gt_len_tol = getattr(args, "gt_len_tol", 50)
        self.nHap = getattr(args, "nHap", 1)
        self.sample_map_file = getattr(args, "sample_map", None)
        self.carrier_info = getattr(args, "carrier_info", None)
        self.size_bins = tuple(getattr(args, "size_bins", None) or DEFAULT_SIZE_BINS)
        self.af_bins = tuple(getattr(args, "af_bins", None) or DEFAULT_AF_BINS)
        self.no_html = getattr(args, "no_html", False)
        self.reference = getattr(args, "reference", None)
        self.nonmobilizing_file = getattr(args, "nonmobilizing", None)
        self.nonmobilizing = load_nonmobilizing(self.nonmobilizing_file) if self.nonmobilizing_file else None
        self.background_sv_file = getattr(args, "background_sv", None)
        self.background_sv = load_background_sv(self.background_sv_file) if self.background_sv_file else None

    def _summarize(self, scored, records, truth_samples):
        '''Every summary table for one scoring of the loci; also returns the unmatched records.'''
        unmatched_records = [r for r in records if not r.matched]
        matched = len(records) - len(unmatched_records)
        unmatched = len(unmatched_records)

        summary = OrderedDict()
        summary["overall"] = aggregate(scored)
        add_unmatched(summary["overall"]["counts"], unmatched_records)
        # A displaced call stays unmatched: it found a locus but put the breakpoint a
        # whole solo LTR away, so it is not a correct prediction and precision must not
        # credit it. Counted separately so it is not read as a spurious call either.
        n_displaced_calls = sum(1 for r in unmatched_records if r.displaced_for is not None)
        summary["predictions"] = OrderedDict([
            ("total", len(records)),
            ("matched", matched),
            ("unmatched", unmatched),
            ("displaced", n_displaced_calls),
            ("precision", round(matched / len(records), 4) if records else None),
        ])
        # Counted per event class, not per record: an element inserted and later excised is
        # one event under both INS and EXC, so its excision is scored as an excision rather
        # than disappearing into the INS row because INFO/TYPE named the event that created
        # the record. Rows therefore do not sum to the event count.
        summary["by_event_type"] = stratify_multi(scored, lambda e: e["event_classes"])
        summary["unsupported_events"] = [
            OrderedDict([("chrom", e["chrom"]), ("pos", e["pos"]), ("id", e["id"]),
                         ("declared", e["unsupported_events"])])
            for e in scored if e["unsupported_events"]]
        if any(e["history"] != e["type"] for e in scored):
            summary["by_event_history"] = stratify(scored, lambda e: e["history"])
        summary["by_family"] = stratify(scored, lambda e: e["family"])
        # The family is parsed exactly as Compare's --TEtype parses it, so the labels here
        # are the values that flag accepts. For a pool whose IDs name the source locus
        # (the MEI callset's chr1-683234-AluSp#SINE/Alu) that makes almost every event its
        # own family, and the superfamily is the level that actually groups them -- so add
        # that table, but only where it groups anything the family table did not.
        superfamilies = {e["superfamily"] for e in scored if e["superfamily"]}
        if superfamilies and len(superfamilies) < len(summary["by_family"]):
            summary["by_superfamily"] = stratify(scored, lambda e: e["superfamily"] or "-")
        summary["by_size"] = _order_bins(
            stratify(scored, lambda e: size_bin(e["size_bp"], self.size_bins)),
            _size_bin_labels(self.size_bins))
        summary["by_allele_frequency"] = _order_bins(
            stratify(scored, lambda e: af_bin(e["allele_frequency"], self.af_bins)),
            _af_bin_labels(self.af_bins))
        if self.nonmobilizing is not None:
            summary["nonmobilizing"] = score_nonmobilizing(self.nonmobilizing, unmatched_records, self.max_dist)
        if self.background_sv is not None:
            # The same charging as for non-mobilizing SVs; the kinds sort by type, then size.
            scored = score_nonmobilizing(self.background_sv, unmatched_records, self.max_dist)
            order = [f"{t} {lab}" for t in ("DEL", "DUP") for _, lab in BG_SV_SIZE_BINS]
            scored["by_kind"] = OrderedDict((k, scored["by_kind"][k]) for k in order if k in scored["by_kind"])
            summary["background_sv"] = scored
        if len(truth_samples) <= MAX_GENOMES_FOR_CARRIER_TABLE:
            summary["by_carrier_count"] = OrderedDict(sorted(
                stratify(scored, lambda e: e["n_carrier_genomes"]).items(),
                key=lambda kv: int(kv[0])))
        return summary, unmatched_records

    def _open_reference(self):
        if not self.reference:
            return None
        try:
            return pysam.FastaFile(self.reference)
        except (OSError, ValueError) as err:
            raise SystemExit(f"[ERROR] --reference {self.reference}: cannot open as an indexed FASTA "
                             f"({err}); index it with samtools faidx, or copy it somewhere writable") from err

    def _run(self):
        # Ground truth. --nHap merges consecutive haplotype columns into one individual,
        # exactly as Compare does, before anything is paired with the prediction.
        truth_file = self.truth_file
        if self.nHap > 1:
            truth_file = f"{self.outprefix}.polyhap.vcf"
            convert_to_ploidy(self.truth_file, self.nHap, truth_file)
        events, truth_samples, truth_rows = load_truth_events(
            truth_file, self.INSonly, self.TEtype)

        if self.predType == "VCF":
            records, pred_samples, pred_rows = load_pred_vcf(self.pred_file, self.carrier_info)
        else:
            records, pred_samples, pred_rows = load_pred_bed(self.pred_file)

        sample_map = read_sample_map(self.sample_map_file) if self.sample_map_file else None
        pairs, pairing = pair_samples(truth_samples, pred_samples, sample_map)

        match_events(events, records, self.max_dist, self.gt_len_tol)
        ordered = sorted(events, key=lambda e: (e.chrom, e.pos))
        near = nearby_events(events, self.max_dist)
        scored = [score_event(event, pairs, self.gt_len_tol, near.get(id(event), ()))
                  for event in ordered]
        # Which record of the truth file each scored locus came from, in the same order.
        # Kept beside the scored loci rather than inside them: the report links a locus to
        # its record with it, and it is a fact about this file rather than about the locus,
        # so it has no business in the JSON a locus is written to.
        locus_rows = [event.index for event in ordered]

        filters = ", ".join(
            part for part in (
                f"--TEtype {self.TEtype}" if self.TEtype else "",
                "--INSonly" if self.INSonly else "") if part)

        # --reference: is each correctly placed call's allele, written into the reference, the
        # simulated haplotype base for base? Recorded on every matched locus, and the basis of
        # the second, sequence-exact evaluation written after the usual one.
        fasta = self._open_reference()
        if fasta is not None:
            for event, entry in zip(ordered, scored):
                if entry["match"] is not None:
                    entry["match"]["haplotype_identical"] = haplotype_identical(event, event.match, pairs, fasta)
        summary, unmatched_records = self._summarize(scored, records, truth_samples)

        # With one allele per genome there is no genotype to get right beyond which
        # haplotype is carried, so the concordance that is reported says "haplotype".
        haploid = all(len(gt) == 1 for e in events for gt in e.gts.values()) and \
                  all(len(gt) == 1 for r in records for gt in r.gts.values())
        meta = OrderedDict([
            ("truth", self.truth_file),
            ("pred", self.pred_file),
            ("predType", self.predType),
            ("n_loci", len(events)),
            ("n_truth_samples", len(truth_samples)),
            ("n_pred_records", len(records)),
            ("n_pred_samples", len(pred_samples)),
            ("sample_pairs", pairs),
            ("pairing", pairing),
            ("genotype_label", "haplotype" if haploid else "genotype"),
            ("filters", filters),
            ("max_dist", self.max_dist),
            ("gt_len_tol", self.gt_len_tol),
            ("nHap", self.nHap),
            ("carrier_info", self.carrier_info),
            ("reference", self.reference),
            ("nonmobilizing", self.nonmobilizing_file),
            ("background_sv", self.background_sv_file),
            ("size_bins", list(self.size_bins)),
            ("af_bins", list(self.af_bins)),
        ])

        report = format_report(summary, meta)
        print(report)

        # One file per locus, alongside the combined summary. The name each event was
        # written under goes into the summary too, so the two outputs can be joined
        # without reconstructing the naming rules.
        out_path = f"{self.outprefix}.json"
        locus_dir = f"{self.outprefix}_loci"
        names = locus_filenames(scored, unmatched_records)
        locus_ref = os.path.basename(locus_dir)
        for event, name in zip(scored, names):
            event["locus_file"] = f"{locus_ref}/{name}"
        run = OrderedDict([
            ("truth", self.truth_file),
            ("pred", self.pred_file),
            ("predType", self.predType),
            ("filters", filters),
            ("max_dist", self.max_dist),
            ("gt_len_tol", self.gt_len_tol),
            ("nHap", self.nHap),
            ("n_paired_genomes", len(pairs)),
            ("summary_file", out_path),
        ])
        removed = write_locus_files(locus_dir, run, scored, unmatched_records, names)

        # The HTML report shows the shape the summary can only average: the per-locus
        # breakpoint and length errors as distributions, and which sizes were missed.
        html_path = None
        if not self.no_html:
            html_path = write_report(f"{self.outprefix}.html", meta, summary, scored,
                                     locus_dir=locus_dir, summary_file=out_path,
                                     truth_rows=truth_rows, pred_rows=pred_rows,
                                     locus_rows=locus_rows)

        # The same evaluation again, counting a call only where its haplotype is the simulated
        # one. A correctly placed call whose sequence differs -- or cannot be compared -- is
        # unpaired from its locus, which becomes a miss, and the record an unmatched
        # prediction; a displaced call loses its detection credit, its allele never having
        # been compared. Done after every output of the usual evaluation is written, because
        # it unpairs the loci those outputs describe.
        exact_summary = exact_scored = exact_path = exact_html = None
        if fasta is not None:
            for event, entry in zip(ordered, scored):
                if event.displaced is not None:
                    event.displaced.displaced_for = None
                    event.displaced = None
                if event.match is not None and not entry["match"].get("haplotype_identical"):
                    event.match.matched = False
                    event.match.matched_for = None
                    event.match = None
            exact_scored = [score_event(event, pairs, self.gt_len_tol, near.get(id(event), ()))
                            for event in ordered]
            for entry, original in zip(exact_scored, scored):
                entry["locus_file"] = original.get("locus_file")
                if entry["match"] is not None:
                    entry["match"]["haplotype_identical"] = True
            exact_summary, exact_unmatched = self._summarize(exact_scored, records, truth_samples)
            print(format_report(exact_summary, meta,
                                title="tevarsim Evaluate -- sequence-exact (--reference)",
                                preamble=[
                                    "  as above, but a call counts only where its allele, written into the reference, is",
                                    "  base-for-base the simulated haplotype; any other call is a miss at its locus and an",
                                    f"  unmatched prediction. reference: {self.reference}"]))
            exact_path = f"{self.outprefix}.exact.json"
            with open(exact_path, "w") as fo:
                json.dump(OrderedDict([("meta", meta), ("summary", exact_summary),
                                       ("loci", exact_scored)]), fo, indent=2)
                fo.write("\n")
            if not self.no_html:
                exact_html = write_report(f"{self.outprefix}.exact.html", meta, exact_summary, exact_scored,
                                          locus_dir=locus_dir, summary_file=exact_path,
                                          truth_rows=truth_rows, pred_rows=pred_rows,
                                          locus_rows=locus_rows)

        with open(out_path, "w") as fo:
            json.dump(OrderedDict([("meta", meta), ("summary", summary)]
                                  + ([("summary_sequence_exact", exact_summary)] if exact_summary else [])
                                  + [("loci", scored)]), fo, indent=2)
            fo.write("\n")

        # Flush first: the report goes to a block-buffered stdout when it is piped, and the
        # unbuffered stderr note would otherwise land above a report it comes after.
        sys.stdout.flush()
        print(f"[INFO] per-event results written to {out_path}", file=sys.stderr)
        print(f"[INFO] {len(names)} per-locus files written to {locus_dir}/ "
              f"({len(scored)} simulated events, {len(unmatched_records)} unmatched "
              f"predictions"
              + (f"; {removed} stale file(s) replaced)" if removed else ")"),
              file=sys.stderr)
        if html_path:
            print(f"[INFO] HTML report written to {html_path}", file=sys.stderr)
        if exact_path:
            print(f"[INFO] sequence-exact evaluation written to {exact_path}"
                  + (f" and {exact_html}" if exact_html else ""), file=sys.stderr)

        self.meta = meta
        self.summary = summary
        self.loci = scored
        self.locus_dir = locus_dir
        self.locus_files = names
        self.html_file = html_path
        self.truth_rows = truth_rows
        self.pred_rows = pred_rows
        self.locus_rows = locus_rows
        self.exact_summary = exact_summary
        self.exact_loci = exact_scored
        return self


def _order_bins(strata, labels):
    '''Put binned strata back in bin order; ``stratify`` sorts by size, which scrambles it.'''
    ordered = OrderedDict((label, strata[label]) for label in labels if label in strata)
    for label, stats in strata.items():
        if label not in ordered:
            ordered[label] = stats
    return ordered


# ---- per-locus files ---------------------------------------------------------

# Contig names end up in a filename, so anything unsafe there is replaced.
_UNSAFE_IN_FILENAME = re.compile(r"[^A-Za-z0-9._-]")


def locus_filenames(scored, unmatched):
    '''
    Name one file per locus, keyed by where the event was simulated: ``<chrom>_<pos>.json``.

    Two elements stacked at one position -- and an unmatched call sitting on a simulated
    event that was paired with a different record -- share a key, so repeats take a
    ``-2``, ``-3`` suffix in the order they are written: simulated events first, then the
    predictions that matched none of them.
    '''
    seen = {}
    names = []
    keyed = [(e["chrom"], e["pos"]) for e in scored] + [(r.chrom, r.pos) for r in unmatched]
    for chrom, pos in keyed:
        key = f"{_UNSAFE_IN_FILENAME.sub('_', str(chrom))}_{pos}"
        seen[key] = seen.get(key, 0) + 1
        names.append(f"{key}.json" if seen[key] == 1 else f"{key}-{seen[key]}.json")
    return names


def nearest_event(record, scored):
    '''
    The closest simulated locus on this contig to a prediction that matched none of them,
    which is what separates a call landing just outside --max_dist from one with nothing
    simulated anywhere near it.
    '''
    best = None
    for event in scored:
        if event["chrom"] != record.chrom:
            continue
        distance = abs(pred_offset(record, event["pos"]))
        if best is None or distance < best[0]:
            best = (distance, event)
    if best is None:
        return None
    distance, event = best
    return OrderedDict([("chrom", event["chrom"]), ("pos", event["pos"]),
                        ("id", event["id"]), ("distance_bp", distance),
                        ("detected", event["detected"])])


def unmatched_payload(record, scored):
    '''Describe a prediction that no simulated event claimed.'''
    carriers = [sample for sample, gt in record.gts.items()
                if pred_is_carrier(record, gt)]
    return OrderedDict([
        ("chrom", record.chrom),
        ("pos", record.pos),
        ("id", record.id),
        ("allele_bp", event_size(record.alt_descs)),
        ("carriers", OrderedDict([("predicted", carriers)])),
        ("nearest_simulated_locus", nearest_event(record, scored)),
        # Set when this call did find a simulated locus, one solo LTR downstream of where
        # it should have anchored. Recorded so it is not mistaken for a spurious call.
        ("displaced_match_for", None if record.displaced_for is None else OrderedDict([
            ("chrom", record.displaced_for.chrom),
            ("pos", record.displaced_for.pos),
            ("id", record.displaced_for.id),
            ("pos_offset", pred_offset(record, record.displaced_for.pos)),
        ])),
    ])


def write_locus_files(directory, run, scored, unmatched, names):
    '''
    Write one self-contained JSON per locus: every simulated event, and every prediction
    that matched none of them. Each file repeats enough of the run to be read on its own
    and names the combined summary it came from; strip ``run`` and ``kind`` and what is
    left is exactly that locus's entry in the summary's ``events`` list.

    Files left by an earlier run of the same prefix would read as this run's results, so
    the directory's own .json files are cleared first. Returns how many were removed.
    '''
    os.makedirs(directory, exist_ok=True)
    removed = 0
    for stale in sorted(os.listdir(directory)):
        path = os.path.join(directory, stale)
        if stale.endswith(".json") and os.path.isfile(path):
            os.remove(path)
            removed += 1

    payloads = [(name, "simulated_locus", event)
                for name, event in zip(names, scored)]
    payloads += [(name, "unmatched_prediction", unmatched_payload(record, scored))
                 for name, record in zip(names[len(scored):], unmatched)]
    for name, kind, body in payloads:
        payload = OrderedDict([("run", run), ("kind", kind)])
        payload.update(body)
        with open(os.path.join(directory, name), "w") as fo:
            json.dump(payload, fo, indent=2)
            fo.write("\n")
    return removed


def _af_bin_labels(edges):
    labels, previous = [], 0.0
    for edge in edges:
        labels.append(f"{previous:g}-{edge:g}")
        previous = edge
    labels.append(f"{previous:g}-1")
    return labels


def run(args):
    Evaluator(args)._run()
