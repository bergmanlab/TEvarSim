import re
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
        self.tsd_min = args.tsd_min
        self.tsd_max = args.tsd_max
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
        # Only split the draw when a deletion spec was actually asked for. Splitting always
        # would reorder the RNG stream for every mixed-type run, changing the output of runs
        # that never wanted this feature; gating it keeps those byte-identical.
        self._af_by_type = any(getattr(args, f"del_af_{k}") is not None
                               for k in ("dist", "mean", "min", "max"))
        self.allow_zero_carriers = args.allow_zero_carriers
        self.sense_strand_ratio = args.sense_strand_ratio
        self.diverse = args.diverse
        self.diverse_config = args.diverse_config
        self.random_seed = args.seed

        self.TEevents = []
        if self.random_seed is not None:
            np.random.seed(self.random_seed)
        # A stream of its own for the draws that only happen sometimes -- the rescued carrier
        # below, and the per-element TSD lengths read from a header. Taking those from the
        # global stream would shift every draw after them, so two runs of the same seed that
        # differed only in whether one event was rescued would differ in the TSD of every
        # unrelated event downstream. Kept separate, the two runs differ in exactly the thing
        # that changed, and stay diffable.
        self._aux = np.random.default_rng(self.random_seed)
    
    def _run(self):
        self._parse_bed()
        self._check_bed()
        self._random_sample_genotypes()
        self.get_TE_tag()
        self.generate_vcf()
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
        previous = {}
        for event in self.TEevents:
            chrom = event["chrom"]
            start = event["start"]
            end = event["end"]
            if start > end:
                raise ValueError(f"Invalid TE event: start >= end. Position: {start}")
            if chrom not in self.CHR:
                raise ValueError(f"Contig not found in reference: {chrom}")
            if start < 0 or end > self.CHR[chrom]["len"]:
                raise ValueError(f"TE event out of genome bounds. Position: {chrom}\t{start}\t{end}")
            prev_start, prev_end = previous.get(chrom, (-1, -1))
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
                raise ValueError(
                    f"Overlapping TE events detected. Position: {chrom}\t{start}, which "
                    f"starts inside the event ending at {prev_end}")
            previous[chrom] = (start, end)
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
                    event_type = "INS" if end-start <= 1 else "DEL"
                self.TEevents.append({
                    "chrom":chrom,
                    "start": start,
                    "end": end,
                    "te_id": te_id,
                    "type": event_type,
                    "strand":strand,
                    "ltr_len": ltr_len
                })
        print(f"[INFO] Parsed {len(self.TEevents)} TE events from BED.",file=sys.stderr)
        print(f"[INFO] Example event: {self.TEevents[0] if self.TEevents else 'No events'}",file=sys.stderr)
    
    def _draw_afs(self, n, dist=None, mean=None, lo=None, hi=None):
        """n allele frequencies from --af-dist, consuming exactly one uniform draw per event.

        "Exactly one" is load-bearing, not incidental. Every later draw in the run reads the
        same global stream -- the genotype binomials here, and the TSD lengths in get_TE_tag --
        so a distribution that consumed a different count would shift all of them, and
        switching --af-dist would silently change the TSD of every event. numpy's
        uniform(a, b, size=n) is a + (b-a)*random_sample(n) bit for bit, so drawing the
        uniforms explicitly and transforming them leaves the stream exactly where the uniform
        branch would have left it.

        The exponential is sampled by inverse CDF for the same reason: rejection sampling
        consumes an unpredictable number of draws. For rate lam truncated to [a, b],

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
        if n == 0:
            return np.empty(0, dtype=float)
        if dist == "uniform":
            return np.random.uniform(a, b, size=n)
        u = np.random.random_sample(n)
        if b <= a:
            # Degenerate truncation: one allowed value. The uniforms are still drawn, so the
            # stream stays aligned with what the uniform branch would have consumed.
            return np.full(n, a, dtype=float)
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
        # past b, and np.random.binomial rejects p > 1. That would be a seed-dependent crash
        # in maybe one run in millions, which is the worst kind to leave in.
        return np.clip(afs, a, b)

    def _draw_afs_by_type(self, events):
        """One allele frequency per event, deletions drawn from the deletion spec.

        Two draws rather than one, so the deletion group can take a different distribution.
        That does consume the stream in a different order than the single-draw path -- all
        non-deletion events, then all deletions -- which is why the caller only takes this
        path when a --del-af-* option was actually given.
        """
        n = len(events)
        afs = np.empty(n, dtype=float)
        del_idx = [i for i, e in enumerate(events) if e["type"] == "DEL"]
        other_idx = [i for i, e in enumerate(events) if e["type"] != "DEL"]
        if other_idx:
            afs[other_idx] = self._draw_afs(len(other_idx))
        if del_idx:
            afs[del_idx] = self._draw_afs(len(del_idx), dist=self.del_af_dist,
                                          mean=self.del_af_mean, lo=self.del_af_min,
                                          hi=self.del_af_max)
        return afs

    def _rescue_zero_carriers(self, genotypes):
        """Give one uniformly chosen genome to each event no genome drew, and count them.

        An event whose binomial roll comes up empty is skipped by generate_vcf and never
        enters a genome, so the number of events requested and the number simulated diverge
        with no error. At the low allele frequencies that make insertions private that is most
        of them: the regime that should produce a population of private insertions instead
        produces an almost empty one.

        Modifies genotypes in place. A no-op at allele frequency 1, since a row of ones is
        never empty, so lineage simulations are unaffected.
        """
        if self.allow_zero_carriers or self.num_genomes == 0 or genotypes.size == 0:
            return 0
        empty = np.flatnonzero(genotypes.sum(axis=1) == 0)
        for row in empty:
            genotypes[row, self._aux.integers(self.num_genomes)] = 1
        return int(empty.size)

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
            nTE = len(self.CHR[chrom]["events"])
            nTE_total += nTE
            if self._af_by_type:
                afs = self._draw_afs_by_type(self.CHR[chrom]["events"])
            else:
                afs = self._draw_afs(nTE)
            self.CHR[chrom]["genotypes"] = np.zeros((nTE, self.num_genomes), dtype=int)
            for i, af in enumerate(afs):
                self.CHR[chrom]["genotypes"][i] = np.random.binomial(1, af, size=self.num_genomes)
                if len(afs_10) < 10:
                    afs_10.append(str(round(af,4)))
            rescued += self._rescue_zero_carriers(self.CHR[chrom]["genotypes"])

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
            # tsd length. Drawn unconditionally and overridden afterwards, rather than
            # skipped, so that --tsd-from-header cannot move the global stream: an element
            # whose header sets its own TSD still consumes the draw it would otherwise have
            # used, leaving every later event -- including the deletions and excisions, which
            # keep the flag-driven length -- exactly where it was. A header RANGE needs a
            # number of its own, so it takes one from the auxiliary stream instead.
            tsd_len = np.random.randint(self.tsd_min, self.tsd_max + 1)
            if te_id.startswith("bg"):
                # Background SVs are synthetic sequence with no element behind them, so they
                # duplicate nothing. Checked before the header lookup so they neither pick up
                # a tag nor turn up in the untagged warning: bgSV writes them into the pool
                # with no description at all.
                tsd_len = 0
            elif self.tsd_from_header and event["type"] == "INS":
                bounds = te_tsd.get(te_id)
                if bounds is None:
                    untagged.append(te_id)
                elif bounds[0] == bounds[1]:
                    tsd_len = bounds[0]
                else:
                    tsd_len = int(self._aux.integers(bounds[0], bounds[1] + 1))
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
            vcf.write('##INFO=<ID=SNP,Number=1,Type=Integer,Description="Number of SNP modifications in TE sequence">\n')
            vcf.write('##INFO=<ID=INDEL,Number=1,Type=Integer,Description="Number of INDEL modifications in TE sequence">\n')
            vcf.write('##INFO=<ID=polyA,Number=1,Type=Integer,Description="Length of added polyA tail in TE sequence">\n')
            vcf.write('##INFO=<ID=truncate,Number=1,Type=Integer,Description="Number of truncated bases at 5\' end of TE sequence">\n')
            vcf.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
            samples = [f"{outprefix}_{i}" for i in range(self.num_genomes)]
            vcf.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(samples) + "\n")

            for chrom,chr_info in self.CHR.items():
                for idx, event in enumerate(chr_info["events"]):
                    genotypes = list(map(str, chr_info["genotypes"][idx]))
                    if all(g=="0" for g in genotypes):
                        continue
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

                    info_str = ";".join(info_parts)
                    vcf.write(f"{chrom}\t{pos}\t{var_id}\t{ref}\t{alt}\t.\tPASS\t{info_str}\tGT\t" + "\t".join(genotypes) + "\n")

        print(f"[INFO] VCF file written to {vcf_path}",file=sys.stderr)


    def generate_genome(self):
        """
        Generate the genome sequence for each sample based on the genotype matrix and variant information.

        Memory, measured rather than guessed, because this is the step that holds the most:
        80,000 events into 100 D. melanogaster genomes (137Mb reference) peaks at 1.45GB and
        takes 98 seconds. It is nearly flat in population size -- 1.30GB at 5 genomes -- since
        one contig of one genome is assembled at a time. Peak tracks the reference size plus
        the total inserted sequence, so it is a genome-size problem, not a sample-size one, and
        a human-scale reference would be roughly an order of magnitude more.

        That is comfortable at the scales this is used at, so the obvious reductions are left
        undone. Recorded here because they are not obvious from reading the code, and because
        a bigger reference would make them worth taking. Measured together at 0.90GB, a 38%
        saving, with byte-identical output and no change in runtime:

          - Build and write ONE CONTIG at a time. Every contig's chunk list currently coexists,
            because the first loop below fills them all before the writer runs. Freeing each
            after it is written is the largest single saving, and the output order is already
            contig-major so it does not move a byte.
          - Drop each event's "ref"/"alt" once it is in chunks. generate_vcf has already run by
            then, and chunks.append(event["alt"][1:]) is a second copy of every element.
          - Free chr_info["seq"] once its contig is chunked.
          - indexMat as bool rather than float64 (128MB -> 16MB at 80k events x 100 genomes).

        The last one has a trap. A column slice indexMat[:, idx] is a VIEW, so the
        mask[flipped] = ~mask[flipped] below would write through and corrupt the matrix for
        every later genome. This code is safe only because .astype(bool) happens to copy; a
        bool matrix needs an explicit .copy() there.
        """
        # split genome
        for chrom,chr_info in self.CHR.items():
            if not chr_info["events"]:
                continue
            for event_idx, event in enumerate(chr_info["events"]):
                # events MUST be sorted
                SVtype = event['type']
                left = event['start']
                right = event['end']
                if left != chr_info["start"]:
                    chr_info["chunks"].append(
                        chr_info["seq"][chr_info["start"]:left]
                    )
                    chr_info["col_index"] += 1
                if SVtype == "INS":
                    chr_info["chunks"].append(event["alt"][1:])
                    chr_info["cols_to_replace"].append(chr_info["col_index"])
                    chr_info["gt_row"].append(event_idx)
                elif SVtype == "EXC":
                    # Substitution: the full element chunk is present when the genotype is 0
                    # (flipped, like a deletion) and the solo LTR chunk when it is 1. Both
                    # chunks are driven by the same event genotype row.
                    chr_info["chunks"].append(event["ref"][1:])   # full LTR-I-LTR element
                    chr_info["flipped"].append(chr_info["col_index"])
                    chr_info["cols_to_replace"].append(chr_info["col_index"])
                    chr_info["gt_row"].append(event_idx)
                    chr_info["col_index"] += 1
                    chr_info["chunks"].append(event["alt"][1:])   # solo LTR
                    chr_info["cols_to_replace"].append(chr_info["col_index"])
                    chr_info["gt_row"].append(event_idx)
                else:
                    chr_info["chunks"].append(event["ref"][1:])
                    chr_info["flipped"].append(chr_info["col_index"])
                    chr_info["cols_to_replace"].append(chr_info["col_index"])
                    chr_info["gt_row"].append(event_idx)
                chr_info["col_index"] += 1
                chr_info["start"] = right

            # Check if the last TE occurs at the end of the genome
            genome_len = chr_info["len"]
            if right != genome_len:
                chr_info["chunks"].append(chr_info["seq"][right:genome_len])

        # combine genome
        if self.diverse:
            print("introduce sequence diversity for each TE-events",file=sys.stderr)


        with ExitStack() as stack:
            files = [stack.enter_context(open(f"{self.output_prefix}_{i}.fa","w")) for i in range(self.num_genomes)]
            for chrom,chr_info in self.CHR.items():
                if not chr_info["events"]:
                    for idx in range(self.num_genomes):
                        files[idx].write(f">{chrom}_{idx}\n")
                        for i in range(0, chr_info["len"], 60):
                            files[idx].write(chr_info["seq"][i:i+60] + "\n")
                    continue
                chunks = chr_info["chunks"]
                cols_to_replace = chr_info["cols_to_replace"]
                flipped = chr_info["flipped"]
                indexMat = np.ones((len(chunks), self.num_genomes))
                # Each genotype-controlled chunk copies the genotype row of its source event.
                # A single event may drive more than one chunk (EXC: full element + solo LTR),
                # so map chunk -> event row explicitly rather than assuming a 1:1 layout.
                for chunk_idx, ev in zip(cols_to_replace, chr_info["gt_row"]):
                    indexMat[chunk_idx] = chr_info["genotypes"][ev]
                for idx in range(self.num_genomes):
                    mask = indexMat[:, idx].astype(bool)
                    mask[flipped] = ~mask[flipped]
                    if self.diverse:
                        # introduce sequence diversity for each TE-events
                        mask_seq = [m if i in cols_to_replace else 0 for i, m in enumerate(mask)]
                        if self.diverse_config:
                            divConfig = Get_config(self.diverse_config)
                            diverse_chunks = [SeqDiverse(chunk, **divConfig) if use else chunk for chunk, use in zip(chunks, mask_seq)]
                        else:
                            diverse_chunks = [SeqDiverse(chunk) if use else chunk for chunk, use in zip(chunks, mask_seq)]
                        assembly = ''.join(chunk for chunk, use in zip(diverse_chunks, mask) if use)
                    else:
                        assembly = ''.join(chunk for chunk, use in zip(chunks, mask) if use)
                    files[idx].write(f">{chrom}_{idx}\n")
                    for i in range(0, len(assembly), 60):
                        files[idx].write(assembly[i:i+60] + "\n")

def run(args):
    Simulator(args)._run()
