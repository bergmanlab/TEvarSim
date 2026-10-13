"""
Tests for the ordering and overlap checks Simulate applies to its BED.

Both checks existed but never fired: the variables they compared against were initialised
once and never updated, so every event was measured against -1 and passed. An unsorted or
overlapping BED was therefore built into genomes in silence, and the genomes were wrong --
generate_genome walks a contig's events with a single forward cursor, so an event that
starts before the previous one ends sends that cursor backwards and duplicates sequence.

The checks are per contig, because that is the scope generate_genome works in. A BED sorted
by contig and then by position -- which is what TErandom writes -- must pass.

No external test runner is required::

    PYTHONPATH=. python tests/test_bed_checks.py

The test functions are also discoverable by pytest.
"""
import os
import sys
import tempfile

# Import the working-tree package regardless of any installed copy.
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from TEvarSim import simulate  # noqa: E402


class _Args:
    def __init__(self, **kw):
        self.__dict__.update(kw)


CHR_A = "ACGT" * 250   # 1000 bp
CHR_B = "TGCA" * 250   # 1000 bp


def _fixture(d, bed_rows):
    ref = os.path.join(d, "ref.fa")
    with open(ref, "w") as f:
        for name, seq in (("chrA", CHR_A), ("chrB", CHR_B)):
            f.write(f">{name}\n")
            for i in range(0, len(seq), 60):
                f.write(seq[i:i + 60] + "\n")
    pool = os.path.join(d, "pool.fa")
    with open(pool, "w") as f:
        f.write(">elem#LTR/Copia\nACGTACGTAC\n")
    bed = os.path.join(d, "events.bed")
    with open(bed, "w") as f:
        for row in bed_rows:
            f.write("\t".join(str(field) for field in row) + "\n")
    return ref, pool, bed


def _simulator(d, bed_rows):
    ref, pool, bed = _fixture(d, bed_rows)
    args = _Args(ref=ref, pool=pool, bed=bed, outprefix=os.path.join(d, "Sim"),
                 num=4, af_min=0.5, af_max=0.5, tsd_min=0, tsd_max=0,
                 sense_strand_ratio=0.5, diverse=False, diverse_config=None, seed=1,
                 af_dist="uniform", af_mean=None, allow_zero_carriers=False,
              del_af_dist=None, del_af_mean=None, del_af_min=None, del_af_max=None, tsd_from_header=False)
    sim = simulate.Simulator(args)
    sim._parse_bed()
    return sim


def _ins(chrom, pos):
    """A point insertion, as TErandom writes one: start == end."""
    return (chrom, pos, pos, "elem#LTR/Copia", ".", "+")


def _del(chrom, start, end):
    return (chrom, start, end, "-", ".", "+")


def _expect_error(sim, wanted):
    try:
        sim._check_bed()
    except ValueError as exc:
        assert wanted in str(exc), f"expected {wanted!r} in error, got: {exc}"
        return str(exc)
    raise AssertionError(f"expected a ValueError mentioning {wanted!r}, none raised")


def test_sorted_non_overlapping_bed_passes():
    """The shape TErandom writes: sorted by contig, then by position."""
    with tempfile.TemporaryDirectory() as d:
        sim = _simulator(d, [_ins("chrA", 100), _del("chrA", 200, 300), _ins("chrA", 400),
                             _ins("chrB", 50), _del("chrB", 600, 700)])
        sim._check_bed()
        assert len(sim.CHR["chrA"]["events"]) == 3
        assert len(sim.CHR["chrB"]["events"]) == 2


def test_contig_boundary_is_not_an_ordering_error():
    """
    The regression the per-contig scope exists to avoid. chrB's first event starts before
    chrA's last one ends, which is fine -- they are different sequences -- and a check held
    globally would reject it.
    """
    with tempfile.TemporaryDirectory() as d:
        sim = _simulator(d, [_del("chrA", 800, 900), _ins("chrB", 10)])
        sim._check_bed()
        assert len(sim.CHR["chrB"]["events"]) == 1


def test_unsorted_within_a_contig_is_rejected():
    with tempfile.TemporaryDirectory() as d:
        sim = _simulator(d, [_ins("chrA", 400), _ins("chrA", 100)])
        _expect_error(sim, "not sorted")


def test_an_insertion_inside_a_deletion_is_nested():
    """An insertion strictly inside a deletion's span landed in the deleted element: allowed,
    and linked to its host, which decides which genomes have its site."""
    with tempfile.TemporaryDirectory() as d:
        sim = _simulator(d, [_del("chrA", 200, 300), _ins("chrA", 250)])
        sim._check_bed()
        dele, ins = sim.CHR["chrA"]["events"]
        assert ins["host"] is dele


def test_a_deletion_overlapping_an_insertion_is_still_rejected():
    """Only an insertion may nest, and only strictly inside: a deletion over an earlier
    insertion's point, or an insertion on a deletion's own breakpoint, is an overlap."""
    with tempfile.TemporaryDirectory() as d:
        sim = _simulator(d, [_del("chrA", 200, 300), _ins("chrA", 250), _del("chrA", 260, 280)])
        _expect_error(sim, "Overlapping")


def test_overlapping_spans_are_rejected():
    with tempfile.TemporaryDirectory() as d:
        sim = _simulator(d, [_del("chrA", 200, 300), _del("chrA", 250, 350)])
        _expect_error(sim, "Overlapping")


def test_two_insertions_at_one_position_are_adjacent_not_overlapping():
    """
    A point event ends where it starts, so a second insertion at the same position abuts the
    first rather than overlapping it, and the genome built from the pair is coherent. The
    comparison stays strict so this keeps passing.
    """
    with tempfile.TemporaryDirectory() as d:
        sim = _simulator(d, [_ins("chrA", 100), _ins("chrA", 100)])
        sim._check_bed()
        assert len(sim.CHR["chrA"]["events"]) == 2


def test_event_touching_the_previous_end_is_allowed():
    """[200,300) then a deletion starting at 300: adjacent, sharing no base."""
    with tempfile.TemporaryDirectory() as d:
        sim = _simulator(d, [_del("chrA", 200, 300), _del("chrA", 300, 400)])
        sim._check_bed()
        assert len(sim.CHR["chrA"]["events"]) == 2


def test_interleaved_contigs_are_allowed():
    """Per-contig scope means a bed need not group its contigs, only order each one."""
    with tempfile.TemporaryDirectory() as d:
        sim = _simulator(d, [_ins("chrA", 100), _ins("chrB", 100),
                             _ins("chrA", 200), _ins("chrB", 200)])
        sim._check_bed()
        assert len(sim.CHR["chrA"]["events"]) == 2
        assert len(sim.CHR["chrB"]["events"]) == 2


def test_error_names_the_contig():
    """Per-contig tracking makes a bare position ambiguous, so the contig is in the message."""
    with tempfile.TemporaryDirectory() as d:
        sim = _simulator(d, [_ins("chrB", 400), _ins("chrB", 100)])
        assert "chrB" in _expect_error(sim, "not sorted")


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
