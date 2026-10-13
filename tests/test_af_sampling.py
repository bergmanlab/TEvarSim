"""
Tests for how Simulate turns an event pool into a population: the allele frequency draw,
and the guarantee that every event reaches at least one genome.

Two things drive these.

The allele frequency was drawn flat between --af-min and --af-max. Over a hundred genomes a
flat spectrum makes most insertions shared, which is the opposite of a real TE landscape --
insertions are overwhelmingly private, with a thin tail of common ones. --af-dist
exponential produces that shape.

And an event whose binomial roll came up empty was skipped by generate_vcf and never entered
a genome, silently, so the number of events requested and the number simulated diverged. At
the low frequencies that make insertions private that is most of them: asking for a private
population produced an almost empty one. Every event now gets a carrier unless
--allow-zero-carriers says otherwise.

Every event draws from streams of its own, seeded by the run's seed and the event, so its
allele frequency, carriers and TSD are the same in every run of the seed that has it, whatever
else the BED holds. The exponential is sampled by inverse CDF, one uniform per event, so that
--af-dist reshapes the frequencies without reordering which events are rare.

No external test runner is required::

    PYTHONPATH=. python tests/test_af_sampling.py

The test functions are also discoverable by pytest.
"""
import io
import os
import sys
import tempfile
from contextlib import redirect_stderr

import numpy as np

# Import the working-tree package regardless of any installed copy.
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from TEvarSim import simulate  # noqa: E402


class _Args:
    def __init__(self, **kw):
        self.__dict__.update(kw)


CHR_A = "ACGT" * 500   # 2000 bp
ELEMENT = "TTTTTGGGGGAAAAACCCCC"


def _args(d, bed_rows, **overrides):
    ref = os.path.join(d, "ref.fa")
    with open(ref, "w") as f:
        f.write(">chrA\n")
        for i in range(0, len(CHR_A), 60):
            f.write(CHR_A[i:i + 60] + "\n")
    pool = os.path.join(d, "pool.fa")
    with open(pool, "w") as f:
        f.write(f">elem#LTR/Copia\n{ELEMENT}\n")
    bed = os.path.join(d, "events.bed")
    with open(bed, "w") as f:
        for row in bed_rows:
            f.write("\t".join(str(field) for field in row) + "\n")
    kw = dict(ref=ref, pool=pool, bed=bed, outprefix=os.path.join(d, "Sim"),
              num=8, af_min=0.1, af_max=0.9, tsd_min=5, tsd_max=20,
              sense_strand_ratio=0.5, diverse=False, diverse_config=None, seed=1,
              af_dist="uniform", af_mean=None, allow_zero_carriers=False,
              del_af_dist=None, del_af_mean=None, del_af_min=None, del_af_max=None,
              tsd_from_header=False)
    kw.update(overrides)
    return _Args(**kw)


def _sampler(**overrides):
    """A Simulator built only far enough to exercise _draw_afs -- no files touched."""
    with tempfile.TemporaryDirectory() as d:
        return simulate.Simulator(_args(d, [], **overrides))


def _events(n, spacing=25):
    """n point insertions, spaced far enough apart to clear the overlap check."""
    return [("chrA", 30 + i * spacing, 30 + i * spacing, "elem#LTR/Copia", ".", "+")
            for i in range(n)]


def _genotype_rows(sim):
    return np.vstack([info["genotypes"] for info in sim.CHR.values() if info["events"]])


def _run_genotypes(d, n, bed_rows=None, **overrides):
    sim = simulate.Simulator(_args(d, _events(n) if bed_rows is None else bed_rows, **overrides))
    sim._parse_bed()
    sim._check_bed()
    sim._random_sample_genotypes()
    return sim


# ---------------------------------------------------------------- the allele frequency draw

def _u(n, seed):
    return np.random.default_rng(seed).random(n)


def test_uniform_is_the_uniform_draw_rescaled():
    """The default spreads each event's own uniform draw linearly over [af_min, af_max]."""
    sim = _sampler(af_dist="uniform", af_min=0.1, af_max=0.9)
    u = _u(50, 7)
    assert np.allclose(sim._draw_afs(u), 0.1 + 0.8 * u)


def test_a_frequency_depends_on_its_own_draw_alone():
    """--af-dist reshapes every frequency without reordering which events are the rare ones:
    each is a monotone function of the same uniform, whichever distribution."""
    uni = _sampler(af_dist="uniform", af_min=0.0, af_max=1.0)
    exp = _sampler(af_dist="exponential", af_mean=0.05, af_min=0.0, af_max=1.0)
    u = _u(50, 11)
    assert np.array_equal(np.argsort(uni._draw_afs(u)), np.argsort(exp._draw_afs(u)))


def test_exponential_matches_the_closed_form():
    """Pins the parameterisation: lam is 1/af_mean, and the inversion is the truncated one."""
    a, b, mean = 0.0, 1.0, 0.05
    sim = _sampler(af_dist="exponential", af_mean=mean, af_min=a, af_max=b)
    u = _u(200, 3)
    got = sim._draw_afs(u)
    lam = 1.0 / mean
    want = np.clip(a - np.log1p(u * np.expm1(-lam * (b - a))) / lam, a, b)
    assert np.array_equal(got, want)


def test_exponential_mean_matches_the_analytic_truncated_mean():
    """The shape check. Fixed seed, so the tolerance is a margin and not a coin flip."""
    a, b, mean = 0.0, 1.0, 0.05
    sim = _sampler(af_dist="exponential", af_mean=mean, af_min=a, af_max=b)
    afs = sim._draw_afs(_u(20000, 19))
    lam = 1.0 / mean
    analytic = a + 1.0 / lam - (b - a) / np.expm1(lam * (b - a))
    assert abs(afs.mean() - analytic) < 0.005, (afs.mean(), analytic)
    # A check on the distribution, not just its mean: an exponential puts 1-1/e of its mass
    # below its own mean, where a uniform over [0,1] would put 0.05 there.
    assert abs((afs < mean).mean() - (1 - 1 / np.e)) < 0.02, (afs < mean).mean()


def test_exponential_respects_the_truncation_bounds():
    """--af-min and --af-max bound the exponential, and nothing escapes into p > 1."""
    sim = _sampler(af_dist="exponential", af_mean=0.05, af_min=0.3, af_max=0.6)
    afs = sim._draw_afs(_u(5000, 5))
    assert afs.min() >= 0.3 and afs.max() <= 0.6, (afs.min(), afs.max())


def test_degenerate_bounds_return_the_value():
    """af_min == af_max leaves one legal frequency."""
    sim = _sampler(af_dist="exponential", af_mean=0.05, af_min=0.4, af_max=0.4)
    afs = sim._draw_afs(_u(20, 13))
    assert np.all(afs == 0.4), afs[:3]


def test_no_events_draws_nothing():
    """A contig with no events gets no frequencies."""
    sim = _sampler(af_dist="exponential", af_mean=0.05)
    assert sim._draw_afs(np.empty(0)).size == 0


def test_an_event_comes_out_the_same_whatever_else_is_in_the_bed():
    """The point of drawing per event: a run that adds events to a BED -- deletions to a set
    of insertions, say -- leaves the carriers of every event the two runs share as they were."""
    with tempfile.TemporaryDirectory() as d1, tempfile.TemporaryDirectory() as d2:
        kw = dict(num=8, af_dist="exponential", af_mean=0.2, af_min=0.0, af_max=1.0,
                  allow_zero_carriers=True)
        every = _events(40)
        alone = _genotype_rows(_run_genotypes(d1, 0, bed_rows=every[::2], **kw))
        mixed = _genotype_rows(_run_genotypes(d2, 0, bed_rows=every, **kw))
        assert np.array_equal(alone, mixed[::2])
        assert alone.any() and not alone.all(), "fixture exercises neither branch"


def test_the_seed_still_decides():
    """Per-event streams are the run's seed and the event, not the event alone."""
    with tempfile.TemporaryDirectory() as d1, tempfile.TemporaryDirectory() as d2:
        kw = dict(num=8, af_min=0.5, af_max=0.5, allow_zero_carriers=True)
        one = _genotype_rows(_run_genotypes(d1, 40, seed=1, **kw))
        two = _genotype_rows(_run_genotypes(d2, 40, seed=2, **kw))
        assert not np.array_equal(one, two)


# ---------------------------------------------------------------- the zero-carrier rescue

def test_every_event_reaches_a_genome_by_default():
    """The point of the rescue: a low mean must give private insertions, not an empty pool."""
    with tempfile.TemporaryDirectory() as d:
        sim = _run_genotypes(d, 60, num=8, af_dist="exponential",
                             af_mean=0.01, af_min=0.0, af_max=1.0)
        rows = _genotype_rows(sim)
        assert rows.shape == (60, 8)
        assert rows.sum(axis=1).min() >= 1, "an event reached no genome"


def test_rescued_events_get_exactly_one_carrier():
    """The rescue hands out one carrier, never a second -- it must not inflate frequencies."""
    with tempfile.TemporaryDirectory() as d:
        sim = _run_genotypes(d, 60, num=8, af_dist="exponential",
                             af_mean=0.001, af_min=0.0, af_max=1.0)
        counts = _genotype_rows(sim).sum(axis=1)
        assert counts.max() == 1, counts.max()


def test_rescue_count_is_reported():
    """A silent rescue would be as misleading as the silent drop it replaces."""
    with tempfile.TemporaryDirectory() as d:
        err = io.StringIO()
        with redirect_stderr(err):
            _run_genotypes(d, 60, num=8, af_dist="exponential",
                           af_mean=0.001, af_min=0.0, af_max=1.0)
        lines = [ln for ln in err.getvalue().splitlines() if "drew no carrier" in ln]
        assert len(lines) == 1, err.getvalue()
        assert "--allow-zero-carriers" in lines[0], lines[0]


def test_rescue_is_a_no_op_at_frequency_one():
    """Lineage simulations pin af to 1, where a row of ones is never empty."""
    with tempfile.TemporaryDirectory() as d:
        err = io.StringIO()
        with redirect_stderr(err):
            sim = _run_genotypes(d, 20, num=8, af_min=1.0, af_max=1.0)
        assert (_genotype_rows(sim) == 1).all()
        assert "drew no carrier" not in err.getvalue()


def test_allow_zero_carriers_restores_the_old_behaviour():
    """The escape hatch really does leave empty rows alone."""
    with tempfile.TemporaryDirectory() as d:
        sim = _run_genotypes(d, 60, num=8, af_dist="exponential", af_mean=0.001,
                             af_min=0.0, af_max=1.0, allow_zero_carriers=True)
        counts = _genotype_rows(sim).sum(axis=1)
        assert counts.min() == 0, "expected some events to reach no genome"


def test_rescue_draws_from_its_own_stream():
    """Two runs differing only in the rescue must differ only in the rows it touched.

    The rescue reads a stream of the event's own, apart from its genotype draw, so the genotype
    rows of every event it did not touch are unchanged.
    """
    with tempfile.TemporaryDirectory() as d1, tempfile.TemporaryDirectory() as d2:
        kw = dict(num=8, af_dist="exponential", af_mean=0.01, af_min=0.0, af_max=1.0)
        rescued = _genotype_rows(_run_genotypes(d1, 60, **kw))
        dropped = _genotype_rows(_run_genotypes(d2, 60, allow_zero_carriers=True, **kw))
        untouched = dropped.sum(axis=1) > 0
        assert untouched.any() and not untouched.all(), "fixture exercises neither branch"
        assert np.array_equal(rescued[untouched], dropped[untouched])


# ---------------------------------------------------------------- deletions get their own AF

def test_reflected_exponential_mirrors_the_exponential():
    """The deletion shape is the insertion shape end for end, so its mean sits at b - mean.

    An ordinary exponential cannot express "most genomes lack this element": asking for a
    high mean just flattens it towards uniform instead of concentrating mass at the top.
    """
    sim = _sampler(af_dist="exponential", af_mean=0.02, af_min=0.0, af_max=1.0)
    u = _u(20000, 31)
    fwd = sim._draw_afs(u)
    rev = sim._draw_afs(u, dist="reflected_exponential")
    assert np.allclose(rev, 0.0 + 1.0 - fwd)
    assert abs(fwd.mean() - 0.02) < 0.005, fwd.mean()
    assert abs(rev.mean() - 0.98) < 0.005, rev.mean()
    assert rev.min() >= 0.0 and rev.max() <= 1.0


def test_deletion_spec_inherits_part_by_part():
    """An unset --del-af-* falls back to the insertion's, so a partial spec is meaningful."""
    with tempfile.TemporaryDirectory() as d:
        sim = simulate.Simulator(_args(d, [], af_dist="exponential", af_mean=0.02,
                                       af_min=0.05, af_max=0.8,
                                       del_af_dist="reflected_exponential"))
        assert sim.del_af_dist == "reflected_exponential"
        assert sim.del_af_mean == 0.02 and sim.del_af_min == 0.05 and sim.del_af_max == 0.8
        assert sim._af_by_type


def test_no_deletion_spec_leaves_the_draw_untouched():
    """Without a --del-af-* option nothing reports a separate deletion spectrum."""
    with tempfile.TemporaryDirectory() as d:
        sim = simulate.Simulator(_args(d, []))
        assert not sim._af_by_type


def test_deletions_draw_high_and_insertions_low():
    """The point of the whole feature, end to end on a mixed BED."""
    # Laid out to fit inside the 2000bp test contig, alternating so the two types are
    # interleaved in BED order, the arrangement most likely to expose a mix-up between them.
    rows = []
    for i in range(12):
        p = 60 + i * 150
        rows.append(("chrA", p, p, "elem#LTR/Copia", ".", "+"))          # INS
        rows.append(("chrA", p + 40, p + 100, "-", ".", "+"))            # DEL
    with tempfile.TemporaryDirectory() as d:
        sim = _run_genotypes(d, 0, num=50, af_dist="exponential", af_mean=0.02,
                             af_min=0.01, af_max=0.99,
                             del_af_dist="reflected_exponential",
                             del_af_min=0.01, del_af_max=0.99, bed_rows=rows)
        rows_gt = _genotype_rows(sim)
        events = [e for info in sim.CHR.values() if info["events"] for e in info["events"]]
        ins = np.array([rows_gt[i].sum() for i, e in enumerate(events) if e["type"] == "INS"])
        dels = np.array([rows_gt[i].sum() for i, e in enumerate(events) if e["type"] == "DEL"])
        assert len(ins) and len(dels), (len(ins), len(dels))
        assert ins.mean() < 10, f"insertions should be rare, got {ins.mean()}"
        assert dels.mean() > 40, f"deletions should be near-fixed, got {dels.mean()}"


if __name__ == "__main__":
    failures = 0
    for name, fn in sorted(globals().items()):
        if name.startswith("test_") and callable(fn):
            try:
                fn()
                print(f"PASS {name}")
            except AssertionError as exc:
                failures += 1
                print(f"FAIL {name}: {exc}")
    print("\nall passed" if not failures else f"\n{failures} failure(s)")
    sys.exit(1 if failures else 0)
