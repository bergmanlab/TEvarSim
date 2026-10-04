"""
Tests for the target site duplication Simulate writes into an insertion.

A TSD is the stretch of host sequence the integration machinery cut twice, so the copy
carried into the ALT allele must be the bases immediately 5' of the insertion point, and
INFO/TSD must be how many of them there are. Neither held.

The slice read seq[start-tsd_len+1 : start+1] -- one base too far right. It duplicated
tsd_len-1 of the left flank's bases and then re-copied seq[start], the first base of the
right flank, so every insertion carried a spurious 1 bp tandem repeat at its 3' junction.
The net inserted length was still tsd_len, which is why nothing caught it: no length check
could see it, and the VCF and the genome were built from the same wrong string, so they
agreed with each other. It was only wrong against the reference.

The same slice went silently empty for an insertion closer to the contig start than
tsd_len, because a negative slice start reads from the far end of the contig.

The second half covers --tsd-from-header. TSD length is a property of the clade, not a
global constant: a cut-and-paste transposon's TSD is set by the stagger between its
transposase's two cuts, so hAT duplicates 8 bp and Tc1/mariner 2 bp every time, while a
non-LTR element inserting by target-primed reverse transcription produces a spread, and a
Helitron duplicates nothing. One --tsd-min/--tsd-max range for a whole library cannot say
that, so the length is read from a TSD= tag on the element's own FASTA header.

No external test runner is required::

    PYTHONPATH=. python tests/test_tsd.py

The test functions are also discoverable by pytest.
"""
import io
import os
import sys
import tempfile
from contextlib import redirect_stderr

# Import the working-tree package regardless of any installed copy.
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from TEvarSim import simulate  # noqa: E402
from TEvarSim.utils import parse_tsd_tag  # noqa: E402


class _Args:
    def __init__(self, **kw):
        self.__dict__.update(kw)


# A non-repetitive contig, so a duplicated run is locatable and unambiguous. "ACGT"*n would
# make every candidate TSD match every other one and the assertions below would pass on the
# broken slice too.
CHR_T = ("GATTACACCGTAAGCTTGGCAATCGTTAAGGCCTATGACTGCAAGTTCCGGATCAGTTAACGGCTTAAGCCT"
         "GTCAAGGTTCCAAGGTTACGCATGGCCTTAAGGCATTCCGGAATTCCGGTTAACCGGTTAACCGGAATTCC"
         "AGGCTTAACCGGATCGATCGATTCCGGAAGGCCTTAACCGGTTAACCGGAATTCCGGAATTCCGGAATTCC")
ELEMENT = "TTTTTGGGGGAAAAACCCCC"

# The insertion point used throughout. Chosen so that CHR_T[INS_POS-1] != CHR_T[INS_POS]:
# the old slice re-copied CHR_T[INS_POS] onto the end of the duplication, and that doubling
# is only detectable where the two colliding bases differ.
INS_POS = 99


def _fixture(d, bed_rows, pool_records=((">elem#LTR/Copia", ELEMENT),)):
    ref = os.path.join(d, "ref.fa")
    with open(ref, "w") as f:
        f.write(">chrT\n")
        for i in range(0, len(CHR_T), 60):
            f.write(CHR_T[i:i + 60] + "\n")
    pool = os.path.join(d, "pool.fa")
    with open(pool, "w") as f:
        for header, seq in pool_records:
            f.write(f"{header}\n{seq}\n")
    bed = os.path.join(d, "events.bed")
    with open(bed, "w") as f:
        for row in bed_rows:
            f.write("\t".join(str(field) for field in row) + "\n")
    return ref, pool, bed


def _args(d, bed_rows, pool_records=((">elem#LTR/Copia", ELEMENT),), **overrides):
    ref, pool, bed = _fixture(d, bed_rows, pool_records)
    kw = dict(ref=ref, pool=pool, bed=bed, outprefix=os.path.join(d, "Sim"),
              num=4, af_min=1.0, af_max=1.0, tsd_min=6, tsd_max=6,
              sense_strand_ratio=0.5, diverse=False, diverse_config=None, seed=1,
                 af_dist="uniform", af_mean=None, allow_zero_carriers=False,
              del_af_dist=None, del_af_mean=None, del_af_min=None, del_af_max=None,
              tsd_from_header=False)
    kw.update(overrides)
    return _Args(**kw)


def _run(d, bed_rows, pool_records=((">elem#LTR/Copia", ELEMENT),), **overrides):
    """Run a whole simulation and hand back the VCF records and genome 0's sequence."""
    args = _args(d, bed_rows, pool_records, **overrides)
    simulate.Simulator(args)._run()
    records = []
    with open(args.outprefix + ".vcf") as f:
        for line in f:
            if not line.startswith("#"):
                records.append(line.rstrip("\n").split("\t"))
    with open(f"{args.outprefix}_0.fa") as f:
        genome = "".join(line.strip() for line in f if not line.startswith(">"))
    return records, genome


def _info(record):
    return dict(field.partition("=")[::2] for field in record[7].split(";"))


def _ins(pos, te_id="elem#LTR/Copia"):
    """A point insertion, as TErandom writes one: start == end."""
    return ("chrT", pos, pos, te_id, ".", "+")


def test_realised_duplication_equals_the_requested_length():
    """The ALT allele ends with the tsd_len bases immediately 5' of the insertion point."""
    pos = INS_POS
    with tempfile.TemporaryDirectory() as d:
        records, genome = _run(d, [_ins(pos)], tsd_min=6, tsd_max=6)
        assert len(records) == 1, records
        alt, info = records[0][4], _info(records[0])
        assert info["TSD"] == "6", info
        assert alt == CHR_T[pos - 1] + ELEMENT + CHR_T[pos - 6:pos], alt
        assert genome == CHR_T[:pos] + ELEMENT + CHR_T[pos - 6:pos] + CHR_T[pos:], genome


def test_no_spurious_tandem_duplication_at_the_junction():
    """The regression guard: the base at the insertion point must not appear twice.

    The old slice re-copied seq[start] on the end of the duplication, so the genome read
    ...<TE><tsd-1 bases><seq[start]><seq[start]>... The fixture position is chosen so that
    the two bases that would collide differ, making the doubling detectable.
    """
    pos = INS_POS
    assert CHR_T[pos - 1] != CHR_T[pos], "fixture must not already repeat at the junction"
    with tempfile.TemporaryDirectory() as d:
        _, genome = _run(d, [_ins(pos)], tsd_min=6, tsd_max=6)
        junction = pos + len(ELEMENT) + 6
        assert genome[junction] != genome[junction - 1], genome[junction - 8:junction + 4]
        assert genome == CHR_T[:pos] + ELEMENT + CHR_T[pos - 6:pos] + CHR_T[pos:]


def test_zero_tsd_inserts_only_the_element():
    """tsd_len 0 means no duplication -- the case the bg-prefix override already relied on."""
    pos = INS_POS
    with tempfile.TemporaryDirectory() as d:
        records, genome = _run(d, [_ins(pos)], tsd_min=0, tsd_max=0)
        assert records[0][4] == CHR_T[pos - 1] + ELEMENT, records[0][4]
        assert genome == CHR_T[:pos] + ELEMENT + CHR_T[pos:]


def test_tsd_near_the_contig_start_is_clamped():
    """An insertion closer to the start than tsd_len duplicates what there is, not nothing.

    A negative slice start reads from the far end of the contig, so the unclamped form
    produced an empty duplication while INFO/TSD still claimed tsd_len bases.
    """
    pos = 3
    with tempfile.TemporaryDirectory() as d:
        records, genome = _run(d, [_ins(pos)], tsd_min=10, tsd_max=10)
        alt = records[0][4]
        assert alt == CHR_T[pos - 1] + ELEMENT + CHR_T[0:pos], alt
        assert genome == CHR_T[:pos] + ELEMENT + CHR_T[0:pos] + CHR_T[pos:]
        # Specifically not the tail of the contig, which is what the negative slice reached for.
        assert CHR_T[-10:] not in alt


# ---------------------------------------------------------------- TSD from the header

COPIA = (">copia#LTR/Copia TSD=4", ELEMENT)
JOCKEY = (">jockey#LINE/Jockey TSD=5-15", ELEMENT)
HELITRON = (">helitron#RC/Helitron TSD=0", ELEMENT)
BARE = (">bare#LTR/Gypsy", ELEMENT)


def test_parse_tsd_tag_grammar():
    """A tag is a fixed length or an inclusive range; a broken one is a typo, not a default."""
    for tail, want in [("TSD=4", (4, 4)), ("TSD=5-15", (5, 15)), ("TSD=0", (0, 0)),
                       ("TSD=0-15", (0, 15)), ("note TSD=8 more", (8, 8)),
                       ("", None), ("some note", None), ("xTSD=3", None)]:
        assert parse_tsd_tag(tail) == want, (tail, parse_tsd_tag(tail))
    for tail in ["TSD=", "TSD=x", "TSD=5-", "TSD=-3", "TSD=9-2"]:
        try:
            parse_tsd_tag(tail, "rec")
        except ValueError:
            continue
        raise AssertionError(f"{tail!r} should not have parsed")


def test_header_value_overrides_the_flags():
    """A tagged element takes its own length, not the one --tsd-min/--tsd-max would give."""
    pos = INS_POS
    with tempfile.TemporaryDirectory() as d:
        records, genome = _run(d, [_ins(pos, "copia#LTR/Copia")], [COPIA],
                               tsd_min=9, tsd_max=9, tsd_from_header=True)
        assert _info(records[0])["TSD"] == "4", _info(records[0])
        assert records[0][4] == CHR_T[pos - 1] + ELEMENT + CHR_T[pos - 4:pos]
        assert genome == CHR_T[:pos] + ELEMENT + CHR_T[pos - 4:pos] + CHR_T[pos:]


def test_header_range_is_drawn_within_its_bounds():
    """A range means a spread, and every drawn length must land inside it."""
    rows = [_ins(30 + i * 25, "jockey#LINE/Jockey") for i in range(6)]
    with tempfile.TemporaryDirectory() as d:
        records, _ = _run(d, rows, [JOCKEY], tsd_min=0, tsd_max=0, tsd_from_header=True)
        lengths = {int(_info(r)["TSD"]) for r in records}
        assert lengths and all(5 <= n <= 15 for n in lengths), lengths
        assert lengths != {0}, "the global flags leaked past the tag"


def test_header_zero_means_no_duplication():
    """TSD=0 is a real value -- a Helitron duplicates nothing -- not a missing one."""
    pos = INS_POS
    with tempfile.TemporaryDirectory() as d:
        records, genome = _run(d, [_ins(pos, "helitron#RC/Helitron")], [HELITRON],
                               tsd_min=8, tsd_max=8, tsd_from_header=True)
        assert _info(records[0])["TSD"] == "0"
        assert records[0][4] == CHR_T[pos - 1] + ELEMENT
        assert genome == CHR_T[:pos] + ELEMENT + CHR_T[pos:]


def test_untagged_element_falls_back_and_warns_once():
    """A mixed library must keep working, and must say which elements it could not read."""
    rows = [_ins(60, "copia#LTR/Copia"), _ins(120, "bare#LTR/Gypsy")]
    with tempfile.TemporaryDirectory() as d:
        err = io.StringIO()
        with redirect_stderr(err):
            records, _ = _run(d, rows, [COPIA, BARE],
                              tsd_min=7, tsd_max=7, tsd_from_header=True)
        warnings = [ln for ln in err.getvalue().splitlines() if "--tsd-from-header" in ln]
        assert len(warnings) == 1, err.getvalue()
        assert "bare#LTR/Gypsy" in warnings[0], warnings[0]
        assert "copia#LTR/Copia" not in warnings[0], warnings[0]
        by_id = {r[2]: int(_info(r)["TSD"]) for r in records}
        assert by_id["copia#LTR/Copia"] == 4, by_id
        assert by_id["bare#LTR/Gypsy"] == 7, by_id


def test_header_tsd_does_not_move_the_global_stream():
    """Reading a tag must not shift the draws of events that do not use one.

    The length is drawn and then overridden rather than skipped, so a deletion -- which has
    no pool record and keeps the flag-driven length -- gets the same number either way. This
    is the assertion that fails if someone later "optimises" that draw away.
    """
    rows = [_ins(60, "copia#LTR/Copia"), ("chrT", 100, 130, "-", ".", "+")]
    with tempfile.TemporaryDirectory() as d1, tempfile.TemporaryDirectory() as d2:
        with redirect_stderr(io.StringIO()):
            on, _ = _run(d1, rows, [COPIA], tsd_min=5, tsd_max=20, tsd_from_header=True)
            off, _ = _run(d2, rows, [COPIA], tsd_min=5, tsd_max=20, tsd_from_header=False)
        deletions = [[r for r in recs if _info(r)["TYPE"] == "DEL"] for recs in (on, off)]
        assert deletions[0] and deletions[1], (on, off)
        assert _info(deletions[0][0])["TSD"] == _info(deletions[1][0])["TSD"]


def test_background_insertions_keep_zero_tsd_and_do_not_warn():
    """bgSV writes its synthetic sequence into the pool untagged; it must not be flagged. It
    duplicates nothing, and it is background, not TE truth: it goes to <prefix>.background.vcf,
    with an ALT of the anchor base and the inserted sequence alone."""
    with tempfile.TemporaryDirectory() as d:
        err = io.StringIO()
        with redirect_stderr(err):
            records, _ = _run(d, [_ins(INS_POS, "bgINS_0_20")],
                              [(">bgINS_0_20", ELEMENT)],
                              tsd_min=8, tsd_max=8, tsd_from_header=True)
        assert records == [], records
        with open(os.path.join(d, "Sim.background.vcf")) as f:
            background = [line.rstrip("\n").split("\t") for line in f if not line.startswith("#")]
        assert len(background) == 1 and background[0][4] == background[0][3] + ELEMENT, background
        assert "--tsd-from-header" not in err.getvalue(), err.getvalue()


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
