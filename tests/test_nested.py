"""
Insertions that land inside a reference element a deletion, excision or non-mobilizing SV takes
in some genomes (NESTED).

Each draws its carriers alone, and a genome without the host sequence loses it and has no site
for it at all: GT "." in the truth, not 0. These tests build small populations and check the
genomes base for base against what each genome's alleles say it should hold.

    PYTHONPATH=. python tests/test_nested.py
"""
import io
import os
import random
import sys
import tempfile
from contextlib import redirect_stderr

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from TEvarSim import simulate  # noqa: E402


class _Args:
    def __init__(self, **kw):
        self.__dict__.update(kw)


def _seq(n, seed):
    rnd = random.Random(seed)
    return "".join(rnd.choice("ACGT") for _ in range(n))


REF = _seq(4000, 1)
ELEM = _seq(300, 2)        # the pool element every insertion uses
CHILD = _seq(120, 3)       # a second family
N = 12


def _run(d, rows, seed=5, **overrides):
    ref = os.path.join(d, "ref.fa")
    with open(ref, "w") as f:
        f.write(">chrA\n" + REF + "\n")
    pool = os.path.join(d, "pool.fa")
    with open(pool, "w") as f:
        f.write(f">elem#LTR/Copia TSD=5\n{ELEM}\n>kid#LINE/L1 TSD=7\n{CHILD}\n")
    bed = os.path.join(d, "events.bed")
    with open(bed, "w") as f:
        for r in rows:
            f.write("\t".join(map(str, r)) + "\n")
    kw = dict(ref=ref, pool=pool, bed=bed, outprefix=os.path.join(d, "Sim"), num=N,
              af_min=0.5, af_max=0.5, tsd_min=None, tsd_max=None, sense_strand_ratio=0.5,
              diverse=False, diverse_config=None, seed=seed, af_dist="uniform", af_mean=None,
              allow_zero_carriers=True, del_af_dist=None, del_af_mean=None, del_af_min=None,
              del_af_max=None, tsd_from_header=True)
    kw.update(overrides)
    with redirect_stderr(io.StringIO()):
        simulate.Simulator(_Args(**kw))._run()
    recs = [l.rstrip("\n").split("\t") for l in open(kw["outprefix"] + ".vcf") if not l.startswith("#")]
    genomes = []
    for i in range(N):
        with open(f"{kw['outprefix']}_{i}.fa") as f:
            genomes.append("".join(l.strip() for l in f if not l.startswith(">")))
    return recs, genomes


def _info(rec):
    return dict(x.partition("=")[::2] for x in rec[7].split(";"))


def _by_id(recs):
    return {r[2]: r for r in recs}


def test_insertion_nested_in_a_deletion():
    rows = [("chrA", 1000, 2000, "DEL-chrA-1000-2000-LTR/Copia-TY1", ".", "+"),
            ("chrA", 1400, 1400, "elem#LTR/Copia", ".", "+")]
    with tempfile.TemporaryDirectory() as d:
        recs, genomes = _run(d, rows)
        dele, ins = _by_id(recs)["DEL-chrA-1000-2000-LTR/Copia-TY1"], _by_id(recs)["elem#LTR/Copia"]
        assert _info(ins)["HOST"] == "chrA:1000:DEL-chrA-1000-2000-LTR/Copia-TY1"
        inserted = ELEM + REF[1395:1400]
        seen = set()
        for g in range(N):
            d_gt, i_gt = dele[9 + g], ins[9 + g]
            seen.add((d_gt, i_gt))
            if d_gt == "1":
                assert i_gt == ".", (g, d_gt, i_gt)
                assert genomes[g] == REF[:1000] + REF[2000:]
            elif i_gt == "1":
                assert genomes[g] == REF[:1400] + inserted + REF[1400:]
            else:
                assert i_gt == "0" and genomes[g] == REF
        assert {("1", "."), ("0", "1"), ("0", "0")} <= seen, seen


def test_insertion_in_an_excised_element_survives_only_in_the_solo_ltr():
    # A 600 bp element whose 5' LTR is its first 100 bp; the excision leaves REF[1000:1100].
    rows = [("chrA", 1000, 1600, "EXC-chrA-1000-1600-100-LTR/Copia-TY1-FULL", ".", "+", 100),
            ("chrA", 1050, 1050, "elem#LTR/Copia", ".", "+"),     # in the 5' LTR: survives
            ("chrA", 1300, 1300, "kid#LINE/L1", ".", "+")]        # internal: excised with it
    with tempfile.TemporaryDirectory() as d:
        recs, genomes = _run(d, rows, seed=11)
        exc, ltr_ins, int_ins = (_by_id(recs)[k] for k in
                                 ("EXC-chrA-1000-1600-100-LTR/Copia-TY1-FULL", "elem#LTR/Copia", "kid#LINE/L1"))
        a = ELEM + REF[1045:1050]
        b = CHILD + REF[1293:1300]
        cases = set()
        for g in range(N):
            e, x, y = exc[9 + g], ltr_ins[9 + g], int_ins[9 + g]
            cases.add((e, x, y))
            assert x != "."                               # the 5' LTR is there in every genome
            if e == "1":
                assert y == "."
                body = REF[1000:1050] + (a if x == "1" else "") + REF[1050:1100]
                assert genomes[g] == REF[:1000] + body + REF[1600:], g
            else:
                body = (REF[1000:1050] + (a if x == "1" else "") + REF[1050:1300] +
                        (b if y == "1" else "") + REF[1300:1600])
                assert genomes[g] == REF[:1000] + body + REF[1600:], g
        assert any(c[0] == "1" and c[1] == "1" for c in cases), cases


def test_a_host_changes_only_the_genomes_it_has_to():
    """Draw, then mask: the nested insertion's carriers with the host's deletion are its carriers
    without it, less the genomes the deletion takes the site from."""
    ins = ("chrA", 1400, 1400, "elem#LTR/Copia", ".", "+")
    with tempfile.TemporaryDirectory() as d1, tempfile.TemporaryDirectory() as d2:
        alone, _ = _run(d1, [ins])
        hosted, _ = _run(d2, [("chrA", 1000, 2000, "DEL-chrA-1000-2000-LTR/Copia-TY1", ".", "+"), ins])
        a, h = _by_id(alone)["elem#LTR/Copia"], _by_id(hosted)["elem#LTR/Copia"]
        assert all(hg == ag for ag, hg in zip(a[9:], h[9:]) if hg != ".")
        assert "." in h[9:]


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
