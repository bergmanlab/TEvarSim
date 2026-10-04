"""
Tests for background variation: Simulate --bg-pi (SNPs and short indels shared by descent) and
TErandom --nSV (bgSV's synthetic insertions and deletions).

No external test runner is required::

    PYTHONPATH=. python tests/test_background.py

The test functions are also discoverable by pytest.
"""
import os
import random
import sys
import tempfile

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from TEvarSim import simulate, utils  # noqa: E402


class _Args:
    def __init__(self, **kw):
        self.__dict__.update(kw)


def _random_seq(n, seed):
    rnd = random.Random(seed)
    return "".join(rnd.choice("ACGT") for _ in range(n))


CONTIGS = {"chrA": _random_seq(30000, 1), "chrB": _random_seq(20000, 2)}
TE_SEQ = _random_seq(400, 3)
# Two insertions and one deletion per contig, well apart.
TE_EVENTS = [("chrA", 5000, 5000, "TE1_0SNP"), ("chrA", 15000, 15600, "DEL"), ("chrA", 25000, 25000, "TE1_0SNP"),
             ("chrB", 4000, 4000, "TE1_0SNP"), ("chrB", 12000, 12500, "DEL"), ("chrB", 17000, 17000, "TE1_0SNP")]


def _write_inputs(d, events=TE_EVENTS):
    ref = os.path.join(d, "ref.fa")
    with open(ref, "w") as f:
        for name, seq in CONTIGS.items():
            f.write(f">{name}\n")
            for i in range(0, len(seq), 60):
                f.write(seq[i:i + 60] + "\n")
    pool = os.path.join(d, "pool.fa")
    with open(pool, "w") as f:
        f.write(f">TE1_0SNP\n{TE_SEQ}\n")
    bed = os.path.join(d, "te.bed")
    with open(bed, "w") as f:
        for chrom, s, e, name in events:
            f.write(f"{chrom}\t{s}\t{e}\t{name}\t.\t+\n")
    return ref, pool, bed


def _args(ref, pool, bed, out, **kw):
    a = dict(ref=ref, pool=pool, bed=bed, outprefix=out, num=12, af_min=0.2, af_max=0.6, tsd_min=5,
             tsd_max=10, sense_strand_ratio=0.5, diverse=False, diverse_config=None, seed=7,
             af_dist="uniform", af_mean=None, allow_zero_carriers=False, tsd_from_header=False,
             del_af_dist=None, del_af_mean=None, del_af_min=None, del_af_max=None)
    a.update(kw)
    return _Args(**a)


def _read_vcf(path):
    samples, recs = None, []
    if not os.path.exists(path):
        return samples, recs
    with open(path) as f:
        for line in f:
            if line.startswith("#CHROM"):
                samples = line.rstrip("\n").split("\t")[9:]
            elif not line.startswith("#"):
                x = line.rstrip("\n").split("\t")
                recs.append((x[0], int(x[1]), x[3], x[4], x[2], x[9:]))
    return samples, recs


def _read_fasta(path):
    seqs, name = {}, None
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                name = line[1:].split()[0]
                seqs[name] = []
            elif name:
                seqs[name].append(line)
    return {k: "".join(v) for k, v in seqs.items()}


def _apply(contig_seq, variants):
    """Apply (pos, ref, alt) records (VCF 1-based POS, REF anchored) to a contig, right to left."""
    s = contig_seq
    for pos, ref, alt in sorted(variants, reverse=True):
        i = pos - 1
        assert s[i:i + len(ref)] == ref, f"REF mismatch at {pos}"
        s = s[:i] + alt + s[i + len(ref):]
    return s


def test_background_genomes_are_reference_plus_both_truths():
    """Every simulated genome is exactly the reference with its TE truth and its background variants applied."""
    with tempfile.TemporaryDirectory() as d:
        ref, pool, bed = _write_inputs(d)
        out = os.path.join(d, "Sim")
        simulate.Simulator(_args(ref, pool, bed, out, bg_pi=0.01))._run()
        samples, te = _read_vcf(out + ".vcf")
        _, bg = _read_vcf(out + ".background.vcf")
        assert bg, "no background variants written"
        assert all(not r[4].startswith("bg") for r in te), "background leaked into the TE truth"
        assert {r[4] for r in bg} <= {"bgSNP", "bgINS", "bgDEL"}
        for gi, sample in enumerate(samples):
            genome = _read_fasta(f"{out}_{gi}.fa")
            for chrom, seq in CONTIGS.items():
                carried = [(p, r, a) for c, p, r, a, _, g in te + bg if c == chrom and g[gi] == "1"]
                assert genome[f"{chrom}_{gi}"] == _apply(seq, carried), f"{sample} {chrom}"


def test_background_is_shared_by_descent_and_clear_of_te_events():
    with tempfile.TemporaryDirectory() as d:
        ref, pool, bed = _write_inputs(d)
        out = os.path.join(d, "Sim")
        simulate.Simulator(_args(ref, pool, bed, out, bg_pi=0.01, bg_margin=30))._run()
        _, bg = _read_vcf(out + ".background.vcf")
        carriers = [sum(g == "1" for g in r[5]) for r in bg]
        assert max(carriers) > 1, "no background variant is shared"
        assert min(carriers) >= 1 and max(carriers) < 12, "a background variant is fixed or absent"
        for chrom, pos, r, a, _, _ in bg:
            lo, hi = pos - 1, pos - 1 + len(r)
            for c, s, e, _ in TE_EVENTS:
                if c == chrom:
                    assert hi + 30 <= s or e + 30 <= lo + 1, f"background at {chrom}:{pos} within the margin of {s}-{e}"


def test_background_off_changes_nothing():
    """--bg-pi 0 (or no option at all) writes byte-identical output: the background draws its own stream."""
    outs = []
    for kw in ({}, {"bg_pi": 0.0}):
        with tempfile.TemporaryDirectory() as d:
            ref, pool, bed = _write_inputs(d)
            out = os.path.join(d, "Sim")
            simulate.Simulator(_args(ref, pool, bed, out, **kw))._run()
            assert not os.path.exists(out + ".background.vcf")
            outs.append((open(out + ".vcf").read(), open(f"{out}_0.fa").read()))
    assert outs[0] == outs[1]


def test_background_leaves_te_events_unchanged():
    """The TE truth (and so its genotypes and TSDs) is the same with and without background."""
    vcfs = []
    for kw in ({}, {"bg_pi": 0.01}):
        with tempfile.TemporaryDirectory() as d:
            ref, pool, bed = _write_inputs(d)
            out = os.path.join(d, "Sim")
            simulate.Simulator(_args(ref, pool, bed, out, **kw))._run()
            vcfs.append(open(out + ".vcf").read())
    assert vcfs[0] == vcfs[1]


def test_bgsv_per_contig_points_and_margins():
    """bgSV draws per contig, writes insertions as points, and keeps clear of TE events and itself."""
    with tempfile.TemporaryDirectory() as d:
        _, pool, bed = _write_inputs(d)
        out_bed, out_fa = os.path.join(d, "bg.bed"), os.path.join(d, "bg.fa")
        random.seed(3)
        utils.bgSV(bed, out_bed, 40, 0.5, pool, out_fa)
        rows = [l.rstrip("\n").split("\t") for l in open(out_bed)]
        bg = [(r[0], int(r[1]), int(r[2]), r[3]) for r in rows if r[3].startswith("bg")]
        assert len(bg) == 40
        assert {c for c, *_ in bg} == {"chrA", "chrB"}, "background SVs not drawn per contig"
        for c, s, e, name in bg:
            if name.startswith("bgINS"):
                assert s == e, f"insertion {name} is not a point"
            assert 0 < s <= e <= len(CONTIGS[c])
            for tc, ts, te, _ in TE_EVENTS:
                if tc == c:
                    assert e + 30 <= ts or te + 30 <= s, f"{name} within the margin of a TE event"
        spans = sorted((c, s, max(e, s + 1)) for c, s, e, _ in bg)
        for (c1, s1, e1), (c2, s2, e2) in zip(spans, spans[1:]):
            assert c1 != c2 or e1 + 30 <= s2, "background SVs overlap"
        ins_names = {n for *_, n in bg if n.startswith("bgINS")}
        assert ins_names <= set(_read_fasta(out_fa)), "insertion sequence missing from the pool"


def test_bgsv_events_go_to_background_vcf():
    """Simulate keeps --nSV's events out of the TE truth and writes them to the background VCF."""
    with tempfile.TemporaryDirectory() as d:
        ref, pool, bed = _write_inputs(d)
        out_bed, out_fa = os.path.join(d, "bg.bed"), os.path.join(d, "bg.fa")
        random.seed(5)
        utils.bgSV(bed, out_bed, 10, 0.5, pool, out_fa)
        out = os.path.join(d, "Sim")
        simulate.Simulator(_args(ref, out_fa, out_bed, out))._run()
        samples, te = _read_vcf(out + ".vcf")
        _, bg = _read_vcf(out + ".background.vcf")
        assert bg and all(r[4].startswith("bg") for r in bg)
        assert not any(r[4].startswith("bg") for r in te)
        for gi in range(len(samples)):
            genome = _read_fasta(f"{out}_{gi}.fa")
            for chrom, seq in CONTIGS.items():
                carried = [(p, r, a) for c, p, r, a, _, g in te + bg if c == chrom and g[gi] == "1"]
                assert genome[f"{chrom}_{gi}"] == _apply(seq, carried)


if __name__ == "__main__":
    for name, fn in list(globals().items()):
        if name.startswith("test_") and callable(fn):
            fn()
            print(f"PASS {name}")
