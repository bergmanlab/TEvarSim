import re
import hashlib
import numpy as np
from Bio import SeqIO
from .utils import description_tail, parse_tsd_tag
from contextlib import ExitStack
from traceback import format_exc
import sys
import os

def SeqDiverse(seq: str,
               snp_rate: float = 0.02,
               indel_rate: float = 0.005,
               ins_ratio: float = 0.5,
               indel_geom_p: float = 0.7,
               max_indel_len: int = 20) -> str:
    seq_array = list(seq)
    L = len(seq_array)
    BASES = ['A', 'C', 'G', 'T']
    CHOICES_DICT = {b: [c for c in BASES if c != b] for b in BASES}
    # INDEL
    n_indel = np.random.poisson(indel_rate * L)
    #print(n_indel)
    n_indel = min(n_indel, int(L * 0.1))
    #print(n_indel)
    if n_indel > 0:
        # generate INDEL lengths
        indel_len_list = np.random.geometric(indel_geom_p, n_indel)
        indel_len_list = np.clip(indel_len_list, a_min=None, a_max=int(max_indel_len))
    # select positions for INDELs
        indel_positions = np.random.choice(L, size=min(n_indel, L), replace=False)
        for pos, indel_len in zip(sorted(indel_positions, reverse=True), indel_len_list):
            if np.random.random() < ins_ratio:
                ins_seq = [np.random.choice(BASES) for _ in range(indel_len)]
                seq_array[pos:pos] = ins_seq
            else:
                del_end = min(pos + indel_len, len(seq_array))
                del seq_array[pos:del_end]
    # SNP
    L = len(seq_array)
    n_snp = np.random.poisson(snp_rate * L)
    if n_snp > 0:
        snp_positions = np.random.choice(L, size=min(n_snp, L), replace=False)
        for pos in snp_positions:
            current_base = seq_array[pos]
            if current_base in CHOICES_DICT:
                seq_array[pos] = np.random.choice(CHOICES_DICT[current_base])
    return ''.join(seq_array)

def Get_config(config_file):
    divConfig = {}
    with open(config_file, 'r') as f:
        for line in f:
            if line.strip() and not line.startswith('#'):
                key, value = line.strip().split('=')
                divConfig[key.strip()] = float(value.strip())
    return divConfig

class Simulator:
    def __init__(self, args):
        self.reference = args.ref
        self.pool_fasta = args.pool
        self.bed_file = args.bed
        self.output_prefix = args.outprefix
        self.num_genomes = args.num
        self.af_min = args.af_min
        self.af_max = args.af_max
        # --tsd-min/--tsd-max name the range an insertion's TSD is drawn from when nothing else
        # says, 5-20 unless given. Given explicitly, they also clamp what a TSD= tag says
        # (tsd_clamp), so one library serves a series of runs that differ only in TSD: --tsd-max 0
        # is the same insertions with no TSD at all.
        given_min = getattr(args, "tsd_min", None)
        given_max = getattr(args, "tsd_max", None)
        self.tsd_clamp = (given_min if given_min is not None else 0, given_max)
        self.tsd_min = given_min if given_min is not None else min(5, given_max if given_max is not None else 5)
        self.tsd_max = given_max if given_max is not None else max(20, self.tsd_min)
        self.tsd_from_header = args.tsd_from_header
        self.af_dist = args.af_dist
        self.af_mean = args.af_mean
        # Deletions get their own allele frequency spec, each part falling back to the
        # insertion's. They need one because the two event types measure opposite things: for
        # an insertion GT=1 means the genome HAS the element, but for a deletion it means the
        # genome LACKS the reference's, so a deletion's allele frequency is the frequency of
        # an ABSENCE. A reference element that most genomes lack -- a reference-private
        # insertion, seen from the other side -- is therefore a deletion at HIGH frequency,
        # which is what --del-af-dist reflected_exponential expresses.
        #
        # Excisions are left on the insertion spec. An excision is a derived event that arises
        # and spreads from rare, exactly like an insertion; only DEL carries the inverted
        # sense.
        self.del_af_dist = args.del_af_dist if args.del_af_dist is not None else args.af_dist
        self.del_af_mean = args.del_af_mean if args.del_af_mean is not None else args.af_mean
        self.del_af_min = args.del_af_min if args.del_af_min is not None else args.af_min
        self.del_af_max = args.del_af_max if args.del_af_max is not None else args.af_max
        self._af_by_type = any(getattr(args, f"del_af_{k}") is not None
                               for k in ("dist", "mean", "min", "max"))
        self.allow_zero_carriers = args.allow_zero_carriers
        self.sense_strand_ratio = args.sense_strand_ratio
        self.diverse = args.diverse
        self.diverse_config = args.diverse_config
        self.random_seed = args.seed
        # Background variation (see _add_background). Read with defaults so callers that build
        # their own args object, the tests among them, need not know about it.
        self.bg_pi = getattr(args, "bg_pi", 0.0) or 0.0
        self.bg_indel_frac = getattr(args, "bg_indel_frac", 0.1)
        self.bg_indel_max = getattr(args, "bg_indel_max", 50)
        self.bg_margin = getattr(args, "bg_margin", 30)
        self.bg_block = getattr(args, "bg_block", 50000)
        # Background structural variants (see _add_background): deletions and tandem duplications.
        self.bg_sv_rate = getattr(args, "bg_sv_rate", 0.0) or 0.0
        self.bg_sv_dup_frac = getattr(args, "bg_sv_dup_frac", 0.0)
        self.bg_sv_min = getattr(args, "bg_sv_min", 50)
        self.bg_sv_max = getattr(args, "bg_sv_max", 10000)

        self.TEevents = []
        if self.random_seed is not None:
            np.random.seed(self.random_seed)
        # Every event draws its allele frequency, genotypes, TSD and rescued carrier from streams
        # of its own, seeded by the run's seed and the event itself (_event_rng). So an event comes
        # out the same in every run of the seed that has it, whatever else the BED holds: adding
        # deletions to a set of insertions, or excisions to both, leaves each insertion's carriers
        # and TSD as they were, and a run differs from another in exactly the events that differ.
        self._run_entropy = (self.random_seed if self.random_seed is not None
                             else np.random.SeedSequence().entropy)
    
    # What an event's stream is for, mixed into its seed alongside the event.
    _AF, _GT, _TSD, _RESCUE = range(4)

    def _event_rng(self, event, purpose):
        """The random stream for one purpose of one event: the run's seed, the event's contig,
        span and ID, and the purpose, hashed together. The purpose keeps one draw from moving
        another, so --num changes no TSD and --tsd-max no genotype."""
        key = f'{event["chrom"]}\t{event["start"]}\t{event["end"]}\t{event["te_id"]}'.encode()
        words = np.frombuffer(hashlib.sha256(key).digest()[:16], dtype=np.uint32).tolist()
        return np.random.default_rng([self._run_entropy, purpose, *words])

    def _clamp_tsd(self, bounds):
        """A TSD= tag's (min, max) inside an explicit --tsd-min/--tsd-max. Clamping both ends keeps
        a tag that falls wholly outside on the nearer bound: TSD=5 under --tsd-max 4 is 4."""
        lo, hi = self.tsd_clamp
        def clamp(x):
            x = max(x, lo)
            return x if hi is None else min(x, hi)
        return clamp(bounds[0]), clamp(bounds[1])
    
    def _run(self):
        self._parse_bed()
        self._check_bed()
        self._random_sample_genotypes()
        self.get_TE_tag()
        self._add_background()
        self.generate_vcf()
        self.generate_background_vcf()
        self.generate_nonmobilizing_vcf()
        self.generate_background_sv_vcf()
        self.generate_genome()

    def _check_bed(self):
        self.CHR = {}
        for ref_record in SeqIO.parse(self.reference, "fasta"):
            self.CHR[ref_record.id] = {
                "len":len(ref_record.seq),
                "seq":str(ref_record.seq),
                "events":[],
                "chunks":[],
                "cols_to_replace":[],
                "gt_row":[],
                "flipped":[],
                "col_index":0,
                "start":0
            }
        # Ordering and overlap are properties of ONE CONTIG's event list, not of the file as a
        # whole: the events are bucketed per contig just below, and generate_genome walks each
        # bucket with its own cursor into that contig's sequence. So the previous event is
        # tracked per contig -- held globally, these checks would reject a bed sorted by contig
        # and then by position, which is exactly what TErandom writes, at every contig boundary.
        # A bed whose contigs are interleaved is still accepted, so long as each contig's own
        # events ascend and none of them overlap.
        #
        # One overlap is allowed: an insertion strictly inside a deletion, excision or
        # non-mobilizing SV (a reference element, or the span a non-mobilizing SV takes) is NESTED
        # in it. A genome without the host sequence cannot carry what was inserted into it;
        # _mask_nested sorts that out once genotypes are drawn.
        previous = {}
        for event in self.TEevents:
            chrom = event["chrom"]
            start = event["start"]
            end = event["end"]
            if start > end:
                raise ValueError(f"Invalid TE event: start >= end. Position: {chrom}\t{start}")
            if chrom not in self.CHR:
                raise ValueError(f"Contig not found in reference: {chrom}")
            if start < 0 or end > self.CHR[chrom]["len"]:
                raise ValueError(f"TE event out of genome bounds. Position: {chrom}\t{start}\t{end}")
            prev_start, prev_end, host = previous.get(chrom, (-1, -1, None))
            # Unsorted is diagnosed before overlapping because it is the more likely cause and
            # the more useful thing to be told: an out-of-order event always overlaps the one
            # before it as well, so the overlap message alone would point at a symptom.
            if start < prev_start:
                raise ValueError(
                    f"bed file is not sorted by start position within a contig. "
                    f"Position: {chrom}\t{start}, which follows {prev_start}")
            # An insertion is a point event (start == end), so two of them at the same position
            # are adjacent rather than overlapping and the comparison stays strict.
            if start < prev_end:
                if not (event["type"] == "INS" and not event["background"] and host is not None
                        and host["start"] < start < host["end"]):
                    raise ValueError(
                        f"Overlapping TE events detected. Position: {chrom}\t{start}, which "
                        f"starts inside the event ending at {prev_end}")
                event["host"] = host
            elif event["type"] != "INS" and not event["background"]:
                host = event
            else:
                host = None
            previous[chrom] = (start, max(end, prev_end), host)
            self.CHR[chrom]["events"].append(event)

    def _parse_bed(self):
        """
        """
        self.TEevents = []
        with open(self.bed_file) as f:
            for line in f:
                if line.startswith("#") or not line.strip():
                    continue
                line = line.strip().split("\t")
                chrom, start, end, te_id, *__ = line
                if len(line) > 5:
                    strand=line[5]
                else:
                    strand = "+" if np.random.rand() < self.sense_strand_ratio else "-"
                te_id = "DEL" if te_id == "-" else te_id # If TE deletion has no ID, assign "DEL"
                start = int(start)
                end = int(end)
                ltr_len = None
                if te_id.startswith("EXC"):
                    # LTR-LTR recombination (excision): the reference span [start,end) is a
                    # full-length LTR element that collapses to a solo LTR. The LTR length is
                    # carried in an optional 7th BED column, falling back to the teID encoding
                    # EXC-{chrom}-{start}-{end}-{ltrlen}-{class/fam}-{name}.
                    event_type = "EXC"
                    if len(line) > 6 and line[6]:
                        ltr_len = int(line[6])
                    else:
                        ltr_len = int(te_id.split("-")[4])
                else:
                    # A non-mobilizing SV (TErandom --nNM) is a deletion of whatever span it names.
                    event_type = "INS" if end-start <= 1 and not te_id.startswith("NM-") else "DEL"
                self.TEevents.append({
                    "chrom":chrom,
                    "start": start,
                    "end": end,
                    "te_id": te_id,
                    "type": event_type,
                    "strand":strand,
                    "ltr_len": ltr_len,
                    # TErandom --nSV's synthetic INS/DEL: in the genomes, but not TE truth.
                    "background": te_id.startswith("bg"),
                    # TErandom --nNM: in the genomes, in <prefix>.nonmobilizing.vcf, never TE truth.
                    "nonmobilizing": te_id.startswith("NM-"),
                })
        print(f"[INFO] Parsed {len(self.TEevents)} TE events from BED.",file=sys.stderr)
        print(f"[INFO] Example event: {self.TEevents[0] if self.TEevents else 'No events'}",file=sys.stderr)
    
    def _draw_afs(self, u, dist=None, mean=None, lo=None, hi=None):
        """Allele frequencies from --af-dist, one per uniform in u (one per event, each from the
        event's own stream).

        One uniform per event, transformed, so that an event's allele frequency is a function of
        its own draw and the spec alone: switching --af-dist reshapes every frequency without
        reshuffling which events are rare.

        The exponential is sampled by inverse CDF for the same reason. For rate lam truncated
        to [a, b],

            F(x) = (1 - exp(-lam*(x-a))) / (1 - exp(-lam*(b-a)))

        inverts to x = a - ln(1 - u*(1 - exp(-lam*(b-a)))) / lam, written below with expm1
        and log1p so it stays accurate when lam*(b-a) is small.

        Note what truncating from below does: an exponential conditioned on X > a is, by
        memorylessness, a plus an exponential. So af_mean is the mean ABOVE af_min, and the
        realised mean allele frequency is about af_min + af_mean, less whatever af_max cuts
        off the tail.
        """
        dist = self.af_dist if dist is None else dist
        mean = self.af_mean if mean is None else mean
        a = self.af_min if lo is None else lo
        b = self.af_max if hi is None else hi
        u = np.asarray(u, dtype=float)
        if dist == "uniform":
            return a + (b - a) * u
        if b <= a:
            # Degenerate truncation: one allowed value.
            return np.full(u.shape, a, dtype=float)
        lam = 1.0 / mean
        afs = a - np.log1p(u * np.expm1(-lam * (b - a))) / lam
        if dist == "reflected_exponential":
            # The same exponential, mirrored end for end within [a, b], so the mass piles up
            # against the TOP of the range instead of the bottom and the mean sits at
            # b - mean rather than a + mean. This is the shape a deletion needs: a reference
            # element that most genomes lack is one whose ABSENCE is at high frequency, and
            # the absence frequency is 1 minus an insertion frequency, which is rare. An
            # ordinary exponential cannot express it -- asking for a high mean just flattens
            # the distribution rather than concentrating it near the top.
            afs = a + b - afs
        # Bounded by [a, b] analytically, but rounding as u approaches 1 can step a few ulp
        # past b, and a frequency above 1 is meaningless wherever it is read afterwards.
        return np.clip(afs, a, b)

    def _event_afs(self, events):
        """One allele frequency per event, deletions drawn from the deletion spec (which is the
        insertion's wherever no --del-af-* option was given)."""
        n = len(events)
        u = np.array([self._event_rng(e, self._AF).random() for e in events], dtype=float)
        afs = np.empty(n, dtype=float)
        # A non-mobilizing SV is a deletion, but a new mutation rather than the absence of a
        # reference element, so it is rare like an insertion: GT=1 is the derived state.
        is_del = np.array([e["type"] == "DEL" and not e.get("nonmobilizing") for e in events], dtype=bool)
        if (~is_del).any():
            afs[~is_del] = self._draw_afs(u[~is_del])
        if is_del.any():
            afs[is_del] = self._draw_afs(u[is_del], dist=self.del_af_dist,
                                         mean=self.del_af_mean, lo=self.del_af_min,
                                         hi=self.del_af_max)
        return afs

    def _mask_nested(self, events, genotypes):
        """Take each nested insertion away from the genomes that lack the sequence it landed in.

        Every event draws its carriers alone, so a nested insertion's drawn genotypes are the same
        whatever its host does, and a run with the host's deletion differs from one without it only
        where it has to. A genome lacks the insertion's site when it carries the host's deletion or
        non-mobilizing SV, or its excision -- unless the insertion lies in the 5' LTR the excision
        leaves. Those genomes are recorded
        as having no site at all (event["applicable"]), which the truth VCF writes as ".", not as
        "does not carry". Events come in file order, so a host's row is final before its nested
        insertions are read against it.
        """
        row_of = {id(e): i for i, e in enumerate(events)}
        for i, e in enumerate(events):
            if e.get("host") is not None:
                h = e["host"]
                carried = genotypes[row_of[id(h)]].astype(bool)
                site = ~carried
                if h["type"] == "EXC" and e["start"] < h["start"] + h["ltr_len"]:
                    site |= carried
            else:
                continue
            e["applicable"] = site
            genotypes[i] &= site

    def _rescue_zero_carriers(self, events, genotypes, nested=None):
        """Give one uniformly chosen genome to each event no genome drew, and count them.

        An event whose binomial roll comes up empty is skipped by generate_vcf and never
        enters a genome, so the number of events requested and the number simulated diverge
        with no error. At the low allele frequencies that make insertions private that is most
        of them: the regime that should produce a population of private insertions instead
        produces an almost empty one.

        Modifies genotypes in place. A no-op at allele frequency 1, since a row of ones is
        never empty, so lineage simulations are unaffected. nested=False takes only events that
        are not nested in another, nested=True only those that are, None every event.
        """
        if self.allow_zero_carriers or self.num_genomes == 0 or genotypes.size == 0:
            return 0
        empty = np.flatnonzero(genotypes.sum(axis=1) == 0)
        rescued = 0
        for row in empty:
            is_nested = events[row].get("host") is not None
            if nested is not None and is_nested != nested:
                continue
            # Only a genome with the site can carry it; a nested insertion none has stays empty.
            site = events[row].get("applicable")
            allowed = np.arange(self.num_genomes) if site is None else np.flatnonzero(site)
            if allowed.size:
                genotypes[row, allowed[self._event_rng(events[row], self._RESCUE).integers(allowed.size)]] = 1
                rescued += 1
        return rescued

    def _random_sample_genotypes(self):
        """
        Randomly generate an allele frequency for each TE event, and then generate the genotype for each sample based on that frequency.
        """
        # 1. allele frequency
        if self.af_min > self.af_max:
            raise ValueError(
        "Parameter error: af_min must be less than or equal to af_max.")

        nTE_total = 0
        rescued = 0
        afs_10 = []

        # 2. sample genotypes
        for chrom in self.CHR:
            events = self.CHR[chrom]["events"]
            nTE = len(events)
            nTE_total += nTE
            afs = self._event_afs(events)
            self.CHR[chrom]["genotypes"] = np.zeros((nTE, self.num_genomes), dtype=int)
            for i, af in enumerate(afs):
                # A Bernoulli(af) per genome, as a uniform below af.
                self.CHR[chrom]["genotypes"][i] = (
                    self._event_rng(events[i], self._GT).random(self.num_genomes) < af)
                if len(afs_10) < 10:
                    afs_10.append(str(round(af,4)))
            # Hosts first, so a rescued carrier of a host is read as one by the masking; then the
            # nested insertions, among the genomes that still have their site.
            gts = self.CHR[chrom]["genotypes"]
            rescued += self._rescue_zero_carriers(events, gts, nested=False)
            self._mask_nested(events, gts)
            rescued += self._rescue_zero_carriers(events, gts, nested=True)

        print(f"[INFO] Generated genotypes for {nTE_total} events across {self.num_genomes} genomes.",file=sys.stderr)
        if self.af_dist != "uniform":
            print(f"[INFO] Allele frequencies drawn from {self.af_dist} of mean {self.af_mean}, "
                  f"truncated to [{self.af_min}, {self.af_max}].",file=sys.stderr)
        if self._af_by_type:
            print(f"[INFO] Deletions drawn separately from {self.del_af_dist} of mean "
                  f"{self.del_af_mean}, truncated to [{self.del_af_min}, {self.del_af_max}]; "
                  f"a deletion's allele frequency is the frequency of the ABSENCE of the "
                  f"reference's element.",file=sys.stderr)
        print(f"[INFO] Allele frequencies (first 10): {' '.join(afs_10)}",file=sys.stderr)
        if rescued:
            print(f"[INFO] {rescued} of {nTE_total} event(s) drew no carrier and were given one "
                  f"at random (allele frequency 1/{self.num_genomes}); pass "
                  f"--allow-zero-carriers to drop them instead.",file=sys.stderr)

    def get_TE_tag(self):
        """
        Generate additional metadata for each TE event:
        strand: the insertion orientation ("+" or "-")
        tsd_len: TSD length (from tsd_min to tsd_max)
        ref: reference sequence
        alt: variant sequence (including the result after TE insertion or deletion)"
        """
        # TE pool. The TSD tags are read into a dict of their own rather than folded into
        # te_pool, so the KeyError handling below keeps working on the shape it expects.
        te_pool = {}
        te_tsd = {}
        for rec in SeqIO.parse(self.pool_fasta, "fasta"):
            te_pool[rec.id] = rec.seq
            if self.tsd_from_header:
                bounds = parse_tsd_tag(description_tail(rec), rec.id)
                if bounds is not None:
                    te_tsd[rec.id] = bounds
        untagged = []
        for event in self.TEevents:
            chrom, start, end, te_id, strand = event["chrom"], event["start"], event["end"], event["te_id"], event["strand"]
            # tsd length, from --tsd-min/--tsd-max, or from the element's TSD= tag clamped to an
            # explicit --tsd-min/--tsd-max. Deletions and excisions duplicate nothing; theirs is
            # drawn from the flags all the same and only ever reaches INFO/TSD.
            lo, hi = self.tsd_min, self.tsd_max
            if te_id.startswith("bg") or event["nonmobilizing"]:
                # Background SVs are synthetic sequence with no element behind them, so they
                # duplicate nothing. Checked before the header lookup so they neither pick up
                # a tag nor turn up in the untagged warning: bgSV writes them into the pool
                # with no description at all.
                lo = hi = 0
            elif self.tsd_from_header and event["type"] == "INS":
                bounds = te_tsd.get(te_id)
                if bounds is None:
                    untagged.append(te_id)
                else:
                    lo, hi = self._clamp_tsd(bounds)
            tsd_len = lo if lo == hi else int(self._event_rng(event, self._TSD).integers(lo, hi + 1))
            if event["type"] == "INS":
                # bed file is 0-based
                ref_allele = self.CHR[chrom]["seq"][start - 1]
                try:
                    te_seq = te_pool[te_id]
                except KeyError:
                    print(format_exc(),file=sys.stderr)
                    raise ValueError(f"TE sequence {te_seq} not found in {self.pool_fasta}. This could be because the pool file doesn't match the input bed.")
                if strand == "-":
                    te_seq = str(te_seq.reverse_complement())
                # A target site duplication is the stretch of host sequence the integration
                # machinery cut twice, so the copy carried into the ALT allele is the tsd_len
                # bases IMMEDIATELY 5' of the insertion point -- the tail of the left flank.
                # The +1-shifted form this replaces took one base too far: it duplicated only
                # tsd_len-1 of them and then re-copied seq[start], the first base of the RIGHT
                # flank, leaving a spurious 1 bp tandem repeat at the junction. The net inserted
                # length was still tsd_len, so nothing that measured length could see it; it
                # showed up only as an off-by-one between INFO/TSD and the genome.
                # max() clamps an insertion closer to the contig start than tsd_len. Without it
                # the slice start goes negative, which Python reads from the far END of the
                # contig and so yields an empty string -- a silently absent TSD on a record
                # whose INFO/TSD claims one.
                tsd_seq = self.CHR[chrom]["seq"][max(0, start - tsd_len) : start]
                alt_allele = self.CHR[chrom]["seq"][start - 1] + str(te_seq) + tsd_seq
            elif event["type"] == "EXC":
                # LTR-LTR recombination: REF is the full LTR-I-LTR element, ALT is the anchor
                # base plus the solo LTR left behind (the element's 5' LTR).
                ltr_len = event["ltr_len"]
                ref_allele = self.CHR[chrom]["seq"][start-1:end]
                solo_ltr = self.CHR[chrom]["seq"][start:start+ltr_len]
                alt_allele = self.CHR[chrom]["seq"][start-1] + solo_ltr
            else:
                # REF sequence
                ref_allele = self.CHR[chrom]["seq"][start-1:end]
                alt_allele = self.CHR[chrom]["seq"][start-1]
            # update TE tags
            event.update({
                "tsd_len": tsd_len,
                "ref": ref_allele,
                "alt": alt_allele
            })
        if untagged:
            names = sorted(set(untagged))
            shown = ", ".join(names[:3])
            more = f" and {len(names) - 3} more" if len(names) > 3 else ""
            print(f"[WARN] --tsd-from-header: {len(untagged)} insertion(s) from {len(names)} "
                  f"pool record(s) carry no TSD= tag and fell back to --tsd-min/--tsd-max "
                  f"({self.tsd_min}-{self.tsd_max}). Examples: {shown}{more}",file=sys.stderr)
    
    _MODIFICATION_SUFFIX = re.compile(r"(?:\d+(?:SNP|INDEL|polyA|truncate))*")

    def _parse_te_modification(self, te_id: str):
        """
        Parse TE ID into family name and modifications.
        Example: "LINE_1SNP0INDEL5polyA" → ("LINE", {"nSNP":1, "npolyA":5})

        Only a suffix that IS a modification string is one. An insertion's ID always ends in
        "_" plus its modifications (possibly none), but a deletion's or an excision's ends in
        the family name, which may hold an underscore of its own: "DEL-...-LTR/Gypsy-TY3_1p-FULL"
        was cut to "...TY3", and an SVA_A deletion to "...SVA", in INFO/TEFAMILY.
        """
        te_family = te_id
        mods = {}
        if "_" in te_id and self._MODIFICATION_SUFFIX.fullmatch(te_id.rsplit("_", 1)[1]):
            te_family, mod_str = te_id.rsplit("_", 1)
            for key in ["SNP", "INDEL", "polyA", "truncate"]:
                match = re.search(rf"(\d+){key}", mod_str)
                if match:
                    val = int(match.group(1))
                    if val > 0:
                        mods["n" + key] = val
        return te_family, mods
    

    @staticmethod
    def _genealogy(n: int, rng) -> list:
        """A random genealogy of n genomes (Kingman coalescent, time in units of the mean pairwise
        coalescence time), as its branches: (length, carriers below it as a bool array). The root,
        above every genome, is not a branch -- a mutation there would be carried by all of them."""
        lineages = [np.zeros(n, dtype=bool) for _ in range(n)]
        for i in range(n):
            lineages[i][i] = True
        lengths = [0.0] * n
        branches = []
        k = n
        while k > 1:
            t = rng.exponential(1.0 / (k * (k - 1) / 2))
            lengths = [x + t for x in lengths]
            i, j = sorted(rng.choice(k, 2, replace=False))
            branches.append((lengths[i], lineages[i]))
            branches.append((lengths[j], lineages[j]))
            merged = lineages[i] | lineages[j]
            for idx in (j, i):
                del lineages[idx]
                del lengths[idx]
            lineages.append(merged)
            lengths.append(0.0)
            k -= 1
        return branches

    def _add_background(self):
        """Background variation: SNPs and short indels shared among the genomes by descent.

        Real genomes at one locus differ by far more than their elements, and a caller has to
        tell the two apart; genomes that differ only by the events simulated never ask it to.
        So, with --bg-pi, each block of --bg-block bp gets a genealogy of the genomes of its own
        (a stand-in for recombination), and mutations fall on its branches at a rate that gives
        pairwise diversity --bg-pi per bp: a mutation is carried by every genome below its
        branch. --bg-indel-frac of them are indels, half insertions of random sequence and half
        deletions, of geometric length up to --bg-indel-max. The reference is the genealogy's
        root, so it carries none of them.

        With --bg-sv-rate, structural variants fall on the same genealogies: deletions and tandem
        duplications (--bg-sv-dup-frac of them), log-uniform in length over --bg-sv-min to
        --bg-sv-max bp, placed anywhere -- over the reference's own elements as readily as anywhere
        else -- but never within --bg-margin of a TE event or of another one, nor over an N. The rate
        is what one genome carries against the reference, per Mb: a genome's mutations are those on
        its path to the root, whose expected length is 2(1 - 1/n) mean pairwise times. They are placed
        before the SNPs and indels, which then keep clear of them as of TE events, so a duplication
        copies the sequence its carriers hold. They are written to <prefix>.background_sv.vcf.

        Background stays --bg-margin bp clear of every TE event, so a TSD and the sequence an
        allele is scored on are the reference's, and is written to <prefix>.background.vcf, not
        the TE truth. It draws from a stream of its own, so the TE events, their genotypes and
        their TSDs are what they would have been without it.
        """
        if not (self.bg_pi or self.bg_sv_rate):
            return
        rng = np.random.default_rng(None if self.random_seed is None else self.random_seed + 7919)
        mu = self.bg_pi / 2          # pairwise diversity = 2 * mutation rate * mean pairwise time (1)
        n_genomes = self.num_genomes
        mu_sv = self.bg_sv_rate / 1e6 / (2 * (1 - 1 / n_genomes)) if n_genomes > 1 else self.bg_sv_rate / 1e6
        log_lo, log_hi = np.log(self.bg_sv_min), np.log(self.bg_sv_max + 1)
        n_sv = 0
        BASES = "ACGT"
        n_bg = 0
        for chrom, info in self.CHR.items():
            seq, L = info["seq"], info["len"]
            blocked = sorted((max(0, e["start"] - self.bg_margin), min(L, e["end"] + self.bg_margin))
                             for e in info["events"])
            candidates, sv_candidates = [], []
            for b0 in range(1, L, self.bg_block):
                b1 = min(L - 1, b0 + self.bg_block)
                if b1 <= b0:
                    continue
                branches = self._genealogy(self.num_genomes, rng)
                lengths = np.array([b[0] for b in branches])
                m = rng.poisson(mu * (b1 - b0) * lengths.sum())
                if m > 0:
                    picks = rng.choice(len(branches), size=m, p=lengths / lengths.sum())
                    for pos, br in zip(rng.integers(b0, b1, size=m), picks):
                        candidates.append((int(pos), branches[br][1]))
                if mu_sv:
                    # Drawn only when asked for, so a run without SVs consumes the stream as before.
                    m_sv = rng.poisson(mu_sv * (b1 - b0) * lengths.sum())
                    if m_sv:
                        picks = rng.choice(len(branches), size=m_sv, p=lengths / lengths.sum())
                        sizes = np.exp(rng.uniform(log_lo, log_hi, size=m_sv)).astype(int)
                        dups = rng.random(m_sv) < self.bg_sv_dup_frac
                        for pos, br, size, dup in zip(rng.integers(b0, b1, size=m_sv), picks, sizes, dups):
                            sv_candidates.append((int(pos), int(size), bool(dup), branches[br][1]))
            candidates.sort(key=lambda c: c[0])
            added, gts = [], []
            # Structural variants first: each keeps --bg-margin from TE events and from the others,
            # and its span joins the blocked set the SNPs and indels then keep clear of.
            sv_spans = []
            for pos, size, dup, carriers in sorted(sv_candidates, key=lambda c: c[0]):
                lo, hi = pos, pos + size
                if lo < 1 or hi > L - 1 or "N" in seq[lo - 1:hi]:
                    continue
                span = (lo - self.bg_margin, hi + self.bg_margin)
                if any(a < span[1] and span[0] < b for a, b in blocked) or \
                        any(a < span[1] and span[0] < b for a, b in sv_spans):
                    continue
                sv_spans.append((lo - self.bg_margin, hi + self.bg_margin))
                if dup:
                    # A tandem duplication: a copy of [lo, hi) inserted right after it.
                    event = {"start": hi, "end": hi, "type": "INS", "ref": seq[hi - 1],
                             "alt": seq[hi - 1] + seq[lo:hi], "sv_kind": "DUP"}
                else:
                    event = {"start": lo, "end": hi, "type": "DEL", "ref": seq[lo - 1:hi],
                             "alt": seq[lo - 1], "sv_kind": "DEL"}
                event.update({"chrom": chrom, "te_id": f"bgSV_{event['sv_kind']}", "strand": "+",
                              "ltr_len": None, "tsd_len": 0, "background": True, "bg_sv": True,
                              "sv_lo": lo, "sv_hi": hi})
                added.append(event)
                gts.append(carriers.astype(int))
                n_sv += 1
            if sv_spans:
                blocked = sorted(blocked + sv_spans)
            prev_end, bi = -1, 0
            for pos, carriers in candidates:
                r = rng.random()
                if r < self.bg_indel_frac / 2:
                    kind, n = "INS", int(min(self.bg_indel_max, rng.geometric(0.3)))
                elif r < self.bg_indel_frac:
                    kind, n = "DEL", int(min(self.bg_indel_max, rng.geometric(0.3)))
                else:
                    kind, n = "SNP", 1
                start, end = (pos, pos) if kind == "INS" else (pos, pos + n)
                # Every event, TE or background, must start after the previous one ends.
                if end > L - 1 or start <= prev_end:
                    continue
                while bi < len(blocked) and blocked[bi][1] <= start:
                    bi += 1
                if bi < len(blocked) and blocked[bi][0] <= max(end, start + 1) - 1:
                    continue
                if kind == "SNP":
                    if seq[pos] not in BASES:
                        continue
                    ref, alt = seq[pos], str(rng.choice([b for b in BASES if b != seq[pos]]))
                elif kind == "INS":
                    ref = seq[pos - 1]
                    alt = ref + "".join(rng.choice(list(BASES), size=n))
                else:
                    if "N" in seq[pos - 1:end]:
                        continue
                    ref, alt = seq[pos - 1:end], seq[pos - 1]
                added.append({"chrom": chrom, "start": start, "end": end, "te_id": f"bg{kind}",
                              "type": kind, "strand": "+", "ltr_len": None, "tsd_len": 0,
                              "ref": ref, "alt": alt, "background": True})
                gts.append(carriers.astype(int))
                prev_end = start if kind == "INS" else end
            if not added:
                continue
            rows = list(zip(info["events"], info["genotypes"])) + list(zip(added, gts))
            rows.sort(key=lambda r: (r[0]["start"], r[0]["end"]))
            info["events"] = [r[0] for r in rows]
            info["genotypes"] = np.array([r[1] for r in rows], dtype=int).reshape(len(rows), self.num_genomes)
            n_bg += len(added)
        print(f"[INFO] Added {n_bg - n_sv} background variants (pairwise diversity {self.bg_pi}/bp, "
              f"{self.bg_indel_frac:.0%} indels), written to {self.output_prefix}.background.vcf",
              file=sys.stderr)
        if self.bg_sv_rate:
            print(f"[INFO] Added {n_sv} background SVs ({self.bg_sv_rate}/genome/Mb against the reference, "
                  f"{self.bg_sv_min}-{self.bg_sv_max} bp, {self.bg_sv_dup_frac:.0%} duplications), written to "
                  f"{self.output_prefix}.background_sv.vcf", file=sys.stderr)

    def generate_background_vcf(self):
        """The background variants (_add_background, and TErandom --nSV's), apart from the TE truth."""
        events = [(chrom, i, e) for chrom, info in self.CHR.items()
                  for i, e in enumerate(info["events"]) if e.get("background") and not e.get("bg_sv")]
        if not events:
            return
        outprefix = os.path.basename(self.output_prefix)
        with open(f"{self.output_prefix}.background.vcf", "w") as vcf:
            vcf.write("##fileformat=VCFv4.2\n")
            for chrom in self.CHR:
                vcf.write(f"##contig=<ID={chrom},length={self.CHR[chrom]['len']}>\n")
            vcf.write('##INFO=<ID=TYPE,Number=1,Type=String,Description="Background variant type: SNP, INS or DEL; not a TE event">\n')
            vcf.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
            samples = [f"{outprefix}_{i}" for i in range(self.num_genomes)]
            vcf.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(samples) + "\n")
            for chrom, idx, e in events:
                genotypes = list(map(str, self.CHR[chrom]["genotypes"][idx]))
                if all(g == "0" for g in genotypes):
                    continue
                # A SNP's REF is its own base (0-based start); INS and DEL carry the base before.
                pos = e["start"] + 1 if e["type"] == "SNP" else e["start"]
                vcf.write(f"{chrom}\t{pos}\t{e['te_id']}\t{e['ref']}\t{e['alt']}\t.\tPASS\tTYPE={e['type']}\tGT\t"
                          + "\t".join(genotypes) + "\n")

    def generate_background_sv_vcf(self):
        """The background structural variants (--bg-sv-rate), apart from the TE truth and the small background:
        what Evaluate --background_sv counts a caller's records against. A duplication's record sits at the
        end of the stretch it copies, its ALT the anchor base plus the copy; INFO/DUPSTART is where the copied
        stretch starts."""
        events = [(chrom, i, e) for chrom, info in self.CHR.items()
                  for i, e in enumerate(info["events"]) if e.get("bg_sv")]
        if not events:
            return
        outprefix = os.path.basename(self.output_prefix)
        path = f"{self.output_prefix}.background_sv.vcf"
        with open(path, "w") as vcf:
            vcf.write("##fileformat=VCFv4.2\n")
            for chrom in self.CHR:
                vcf.write(f"##contig=<ID={chrom},length={self.CHR[chrom]['len']}>\n")
            vcf.write('##INFO=<ID=TYPE,Number=1,Type=String,Description="Background structural variant: DEL, or DUP '
                      '(a tandem duplication); not a TE event">\n')
            vcf.write('##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Net length change">\n')
            vcf.write('##INFO=<ID=DUPSTART,Number=1,Type=Integer,Description="DUP only: 1-based first base of the '
                      'copied stretch, which ends at POS">\n')
            vcf.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
            samples = [f"{outprefix}_{i}" for i in range(self.num_genomes)]
            vcf.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(samples) + "\n")
            n = 0
            for chrom, idx, e in events:
                genotypes = list(map(str, self.CHR[chrom]["genotypes"][idx]))
                if all(g == "0" for g in genotypes):
                    continue
                n += 1
                info = f"TYPE={e['sv_kind']};SVLEN={len(e['alt']) - len(e['ref'])}"
                if e["sv_kind"] == "DUP":
                    info += f";DUPSTART={e['sv_lo'] + 1}"
                vcf.write(f"{chrom}\t{e['start']}\tbgSV_{e['sv_kind']}_{n}\t{e['ref']}\t{e['alt']}\t.\tPASS\t{info}\tGT\t"
                          + "\t".join(genotypes) + "\n")
        print(f"[INFO] {n} segregating background SV(s) written to {path}", file=sys.stderr)

    def generate_nonmobilizing_vcf(self):
        """The non-mobilizing SVs (TErandom --nNM), apart from the TE truth: what Evaluate --nonmobilizing
        counts a caller's records against. Written only when the BED holds any."""
        events = [(chrom, i, e) for chrom, info in self.CHR.items()
                  for i, e in enumerate(info["events"]) if e.get("nonmobilizing")]
        if not events:
            return
        outprefix = os.path.basename(self.output_prefix)
        path = f"{self.output_prefix}.nonmobilizing.vcf"
        with open(path, "w") as vcf:
            vcf.write("##fileformat=VCFv4.2\n")
            for chrom in self.CHR:
                vcf.write(f"##contig=<ID={chrom},length={self.CHR[chrom]['len']}>\n")
            vcf.write('##INFO=<ID=TYPE,Number=1,Type=String,Description="Variant type: DEL, a deletion that cuts into or '
                      'swallows a reference element by no TE mechanism (TErandom --nNM); a TE caller should call none">\n')
            vcf.write('##INFO=<ID=NMKIND,Number=1,Type=String,Description="internal (inside the element), left (left flank '
                      'into the element), right (element into the right flank), tight or wide (the whole element and a '
                      'little or a lot of flank)">\n')
            vcf.write('##INFO=<ID=NMHOST,Number=1,Type=String,Description="The reference element it cuts into or swallows: '
                      'chrom-start-end-class/family-name, BED coordinates">\n')
            vcf.write('##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Net length change">\n')
            vcf.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
            samples = [f"{outprefix}_{i}" for i in range(self.num_genomes)]
            vcf.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(samples) + "\n")
            n = 0
            for chrom, idx, e in events:
                genotypes = list(map(str, self.CHR[chrom]["genotypes"][idx]))
                if all(g == "0" for g in genotypes):
                    continue
                _, kind, host = e["te_id"].split("-", 2)
                vcf.write(f"{chrom}\t{e['start']}\t{e['te_id']}\t{e['ref']}\t{e['alt']}\t.\tPASS\t"
                          f"TYPE=DEL;NMKIND={kind};NMHOST={host};SVLEN={len(e['alt']) - len(e['ref'])}\tGT\t"
                          + "\t".join(genotypes) + "\n")
                n += 1
        print(f"[INFO] {n} segregating non-mobilizing SV(s) written to {path}", file=sys.stderr)

    def generate_vcf(self):
        """
        Generate a VCF file from parsed TE events with detailed INFO.
        INFO includes:
            - TYPE: INS or DEL
            - STRAND: +/-
            - TSD: Target site duplication length
            - TEFAMILY: TE family ID (from consensus)
            - SNP/INDEL/polyA/truncate: modification counts (only for INS)
        """
        vcf_path = f"{self.output_prefix}.vcf"
        outprefix = os.path.basename(self.output_prefix)
        with open(vcf_path, "w") as vcf:
            # VCF Header
            vcf.write("##fileformat=VCFv4.2\n")
            for chr in self.CHR:
                vcf.write(f"##contig=<ID={chr},length={self.CHR[chr]['len']}>\n")
            vcf.write('##INFO=<ID=TYPE,Number=1,Type=String,Description="Variant type (INS, DEL, or EXC for LTR-LTR recombination)">\n')
            vcf.write('##INFO=<ID=EVENTTYPE,Number=A,Type=String,Description="Type of the event that produced each ALT allele (INS, DEL, or EXC), one per ALT allele. A single simulated generation only ever produces one event per record, but an allele may accumulate further events across generations, e.g. EVENTTYPE=INS,EXC for an element that was inserted and later excised to a solo LTR">\n')
            vcf.write('##INFO=<ID=STRAND,Number=1,Type=String,Description="Insertion strand (+ or -)">\n')
            vcf.write('##INFO=<ID=LTRLEN,Number=1,Type=Integer,Description="Length of the solo LTR left behind by LTR-LTR recombination (EXC only)">\n')
            vcf.write('##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Net length change of the variant, negative for excision (EXC only)">\n')
            vcf.write('##INFO=<ID=TSD,Number=1,Type=Integer,Description="Target site duplication length: for an insertion the ALT allele ends with the TSD bases immediately 5\' of POS. Meaningless for DEL and EXC, which duplicate nothing">\n')
            vcf.write('##INFO=<ID=TEFAMILY,Number=1,Type=String,Description="TE family ID from consensus">\n')
            vcf.write('##INFO=<ID=HOST,Number=1,Type=String,Description="The event this insertion landed in, as CHROM:POS:ID: a reference event (DEL, EXC, or a non-mobilizing SV in the .nonmobilizing.vcf) whose span holds POS. GT is . in a genome without the host sequence (no site to carry), 0 in one with the host sequence but not this insertion">\n')
            vcf.write('##INFO=<ID=SNP,Number=1,Type=Integer,Description="Number of SNP modifications in TE sequence">\n')
            vcf.write('##INFO=<ID=INDEL,Number=1,Type=Integer,Description="Number of INDEL modifications in TE sequence">\n')
            vcf.write('##INFO=<ID=polyA,Number=1,Type=Integer,Description="Length of added polyA tail in TE sequence">\n')
            vcf.write('##INFO=<ID=truncate,Number=1,Type=Integer,Description="Number of truncated bases at 5\' end of TE sequence">\n')
            vcf.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
            samples = [f"{outprefix}_{i}" for i in range(self.num_genomes)]
            vcf.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(samples) + "\n")

            for chrom,chr_info in self.CHR.items():
                for idx, event in enumerate(chr_info["events"]):
                    if event.get("background") or event.get("nonmobilizing"):
                        continue    # written to the background or non-mobilizing VCF instead: not TE truth
                    row = chr_info["genotypes"][idx]
                    if not row.any():
                        continue
                    site = event.get("applicable")
                    genotypes = [str(g) if site is None or site[i] else "." for i, g in enumerate(row)]
                    # Although BED files are 0-based and VCF files are 1-based, 
                    # VCF files require retrieving one base upstream, so the overall position remains unchanged.
                    pos = event["start"]
                    var_id = event["te_id"]
                    ref = event["ref"]
                    alt = str(event["alt"])
                    event_type = event["type"]
                    strand = event["strand"]
                    tsd = event["tsd_len"]

                    # parse TEID: TE family + modified information
                    te_family, mods = self._parse_te_modification(var_id)

                    # INFO
                    # EVENTTYPE is Number=A, so it carries one value per ALT allele. Every
                    # record written here is single-ALT (one event per site per generation),
                    # so it gets exactly one value; downstream tools that merge an element's
                    # later history onto the same record extend the list allele by allele.
                    info_parts = [
                        f"TYPE={event_type}",
                        f"EVENTTYPE={event_type}",
                        f"STRAND={strand}",
                        f"TSD={tsd}",
                        f"TEFAMILY={te_family}"
                    ]
                    if event_type == "INS":  # only insert needs to be changed
                        for key in ["SNP", "INDEL", "polyA", "truncate"]:
                            val = mods.get("n" + key)
                            if val is not None:
                                info_parts.append(f"{key}={val}")
                    elif event_type == "EXC":
                        info_parts.append(f"LTRLEN={event['ltr_len']}")
                        info_parts.append(f"SVLEN={len(alt) - len(ref)}")
                    host = event.get("host")
                    if host is not None:
                        info_parts.append(f"HOST={host['chrom']}:{host['start']}:{host['te_id']}")

                    info_str = ";".join(info_parts)
                    vcf.write(f"{chrom}\t{pos}\t{var_id}\t{ref}\t{alt}\t.\tPASS\t{info_str}\tGT\t" + "\t".join(genotypes) + "\n")

        print(f"[INFO] VCF file written to {vcf_path}",file=sys.stderr)


    def generate_genome(self):
        """
        Generate the genome sequence for each sample based on the genotype matrix and variant information.

        Memory, as measured before nested insertions: 80,000 events into 100 D. melanogaster
        genomes (137Mb reference) peaked at 1.45GB in 98 seconds, nearly flat in population size
        (1.30GB at 5 genomes), tracking the reference plus the total inserted sequence. Chunks
        are now built and written one contig at a time (_contig_chunks), the largest of the
        savings that measurement identified; each chunk holds a bool presence vector, or None
        for a chunk every genome holds.
        """
        if self.diverse:
            print("introduce sequence diversity for each TE-events",file=sys.stderr)
        divConfig = Get_config(self.diverse_config) if self.diverse and self.diverse_config else {}
        with ExitStack() as stack:
            files = [stack.enter_context(open(f"{self.output_prefix}_{i}.fa","w")) for i in range(self.num_genomes)]
            for chrom, chr_info in self.CHR.items():
                chunks = self._contig_chunks(chr_info)
                for idx in range(self.num_genomes):
                    present = [(c, te) for c, p, te in chunks if p is None or p[idx]]
                    if self.diverse:
                        # sequence diversity for each TE event's own sequence (not background's)
                        assembly = ''.join(SeqDiverse(c, **divConfig) if te else c for c, te in present)
                    else:
                        assembly = ''.join(c for c, _ in present)
                    files[idx].write(f">{chrom}_{idx}\n")
                    for i in range(0, len(assembly), 60):
                        files[idx].write(assembly[i:i+60] + "\n")

    def _contig_chunks(self, chr_info):
        """One contig as (sequence, presence, is_te) chunks: presence is the bool vector of the
        genomes that hold the chunk, or None for every genome; is_te marks a TE event's own
        sequence, which --diverse mutates.

        Walking the events in order, each top-level event contributes its alleles' chunks: an
        insertion its sequence, a
        deletion or non-mobilizing SV its span (held by the genomes that do not carry it), an
        excision the full element and the solo LTR, a background SNP its two bases. An insertion
        nested in a reference event splits that event's span at its point, in each allele that
        still has the point -- the solo LTR as well as the full element when it lies in the 5'
        LTR -- and is held by its carriers among the genomes holding that allele.
        """
        seq, events = chr_info["seq"], chr_info["events"]
        genotypes = chr_info["genotypes"]
        rows = {id(e): genotypes[i].astype(bool) for i, e in enumerate(events)}
        nested_in, tops = {}, []
        for e in events:
            if e.get("host") is not None:
                nested_in.setdefault(id(e["host"]), []).append(e)
            else:
                tops.append(e)
        chunks = []

        def both(a, b):
            return a if b is None else a & b

        def emit_ins(e, cond):
            chunks.append((e["alt"][1:], both(rows[id(e)], cond), not e.get("background")))

        def emit_span(lo, hi, held, inner, te):
            cur = lo
            for x in inner:
                if lo < x["start"] < hi:
                    chunks.append((seq[cur:x["start"]], held, te))
                    emit_ins(x, held)
                    cur = x["start"]
            chunks.append((seq[cur:hi], held, te))

        cursor = 0
        for e in tops:
            lo, hi = e["start"], e["end"]
            if lo > cursor:
                chunks.append((seq[cursor:lo], None, False))
            te = not (e.get("background") or e.get("nonmobilizing"))
            carried = rows[id(e)]
            inner = nested_in.get(id(e), [])
            if e["type"] == "INS":
                emit_ins(e, None)
            elif e["type"] == "SNP":
                # Background substitution: the reference base without it, the other base with it.
                chunks.append((e["ref"], ~carried, te))
                chunks.append((e["alt"], carried, te))
            elif e["type"] == "EXC":
                emit_span(lo, hi, ~carried, inner, te)                         # the full element
                emit_span(lo, lo + e["ltr_len"], carried, inner, te)           # the solo LTR left
            else:
                emit_span(lo, hi, ~carried, inner, te)
            cursor = max(cursor, hi)
        if cursor < chr_info["len"]:
            chunks.append((seq[cursor:chr_info["len"]], None, False))
        return chunks

def run(args):
    Simulator(args)._run()
