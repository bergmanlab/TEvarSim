"""
TErandom makes its events in layers -- insertions, then deletions, then excisions, then
non-mobilizing SVs -- each from a stream of its own, each keeping clear of the layers before it.
So runs that share a seed and differ in one layer share every layer before it, and Simulate,
drawing each event from streams of its own, gives a shared event the same carriers and TSD.
That is what lets a series of simulations vary one thing at a time: the same insertions with
and without deletions, with and without excisions.

    PYTHONPATH=. python tests/test_layers.py
"""
import os
import random
import subprocess
import sys
import tempfile

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
REF_LEN = 400000
LTR_LEN, INTERNAL_LEN = 300, 1500
ELEN = 2 * LTR_LEN + INTERNAL_LEN


def _seq(n, seed):
    rnd = random.Random(seed)
    return "".join(rnd.choice("ACGT") for _ in range(n))


def _tevarsim(*args):
    run = subprocess.run([sys.executable, "-m", "TEvarSim", *map(str, args)], cwd=ROOT,
                         capture_output=True, text=True, env={**os.environ, "PYTHONPATH": ROOT})
    assert run.returncode == 0, run.stderr[-2000:]


def _fixture(d):
    ref = os.path.join(d, "ref.fa")
    with open(ref, "w") as f:
        f.write(">chrT\n" + _seq(REF_LEN, 3) + "\n")
    with open(ref + ".fai", "w") as f:
        f.write(f"chrT\t{REF_LEN}\t6\t{REF_LEN}\t{REF_LEN + 1}\n")
    lib = os.path.join(d, "lib.fa")
    with open(lib, "w") as f:
        f.write(f">TY1-LTR#LTR/Copia TSD=5\n{_seq(LTR_LEN, 1)}\n>TY1-I#LTR/Copia\n{_seq(INTERNAL_LEN, 2)}\n"
                f">jockey#LINE/Jockey TSD=5-15\n{_seq(2000, 4)}\n")
    # Reference elements every 10 kb: deletion candidates span element + 5 bp TSD copy,
    # excision candidates the element alone.
    dels, excs = os.path.join(d, "del.bed"), os.path.join(d, "exc.bed")
    with open(dels, "w") as fd, open(excs, "w") as fe:
        for start in range(10000, REF_LEN - 10000, 10000):
            fd.write(f"chrT\t{start}\t{start + ELEN + 5}\tTY1-FULL#LTR/Copia\t.\t+\n")
            fe.write(f"chrT\t{start}\t{start + ELEN}\tTY1-FULL#LTR/Copia\t.\t+\t{LTR_LEN}\n")
    return ref, lib, dels, excs


def _pool(d, name, n_del=0, n_exc=0, n_nm=0, seed=7, extra=()):
    ref, lib, dels, excs = _fixture(d)
    out = os.path.join(d, name)
    _tevarsim("TErandom", "--ref", ref, "--consensus", lib, "--existingTEs", dels,
              "--excCandidates", excs, "--nINS", 20, "--nDEL", n_del, "--nEXC", n_exc,
              "--nNM", n_nm, "-A", 0, "-T", 0, "-S", 0, "-I", 0, "--seed", seed, "--outprefix", out,
              *extra)
    rows = [l.rstrip("\n").split("\t") for l in open(out + ".bed")]
    return out, rows


def _of(rows, kind):
    if kind == "INS":
        return [r[:4] + r[5:6] for r in rows if int(r[2]) - int(r[1]) <= 1 and not r[3].startswith("NM-")]
    return [r[:4] for r in rows if r[3].startswith(kind)]


def test_insertions_are_the_same_whatever_else_is_asked_for():
    with tempfile.TemporaryDirectory() as d:
        base = _of(_pool(d, "base")[1], "INS")
        assert len(base) == 20, base
        for kw in (dict(n_del=5), dict(n_exc=4), dict(n_nm=6), dict(n_del=5, n_exc=4, n_nm=6)):
            got = _of(_pool(d, "v" + "_".join(f"{k}{v}" for k, v in kw.items()), **kw)[1], "INS")
            assert got == base, (kw, got, base)


def test_strand_targets_do_not_tie_the_insertions_to_the_deletions():
    with tempfile.TemporaryDirectory() as d:
        ratio = ("--sense-strand-ratio", 0.7)
        base = _of(_pool(d, "a", extra=ratio)[1], "INS")
        assert sum(r[4] == "+" for r in base) == 14, base
        assert _of(_pool(d, "b", n_del=5, extra=ratio)[1], "INS") == base


def test_more_deletions_keep_the_ones_fewer_made():
    with tempfile.TemporaryDirectory() as d:
        three = _of(_pool(d, "a", n_del=3)[1], "DEL")
        six = _of(_pool(d, "b", n_del=6)[1], "DEL")
        assert len(three) == 3 and len(six) == 6
        assert set(map(tuple, three)) <= set(map(tuple, six)), (three, six)


def test_excisions_leave_the_deletions_alone_and_nest():
    with tempfile.TemporaryDirectory() as d:
        dels = _of(_pool(d, "a", n_del=5)[1], "DEL")
        _, rows2 = _pool(d, "b", n_del=5, n_exc=2)
        _, rows4 = _pool(d, "c", n_del=5, n_exc=4)
        assert _of(rows2, "DEL") == dels and _of(rows4, "DEL") == dels
        two, four = _of(rows2, "EXC"), _of(rows4, "EXC")
        assert len(two) == 2 and len(four) == 4
        assert set(map(tuple, two)) <= set(map(tuple, four)), (two, four)
        for e in four:     # never the element a deletion took
            assert not any(e[0] == x[0] and int(e[1]) < int(x[2]) and int(x[1]) < int(e[2]) for x in dels)


def test_later_layers_keep_clear_of_insertions():
    with tempfile.TemporaryDirectory() as d:
        _, rows = _pool(d, "a", n_del=10, n_exc=10, n_nm=6)
        points = [int(r[1]) for r in _of(rows, "INS")]
        spans = [(int(r[1]), int(r[2])) for r in rows if r[3][:3] in ("DEL", "EXC", "NM-")]
        for s, e in spans:
            assert not any(s - 1 <= p < e + 1 for p in points), (s, e, points)


def test_simulate_gives_a_shared_insertion_the_same_carriers_and_tsd():
    """The whole pipeline: adding deletions and excisions changes no insertion's record."""
    with tempfile.TemporaryDirectory() as d:
        ref = _fixture(d)[0]
        recs = []
        for name, kw in (("a", {}), ("b", dict(n_del=5, n_exc=3))):
            out, _ = _pool(d, name, **kw)
            _tevarsim("Simulate", "--ref", ref, "--bed", out + ".bed", "--pool", out + ".fa",
                      "--num", 10, "--af-dist", "exponential", "--af-mean", 0.2, "--af-min", 0.1,
                      "--allow-zero-carriers", "--tsd-from-header", "--seed", 7, "-O", out + "_sim")
            with open(out + "_sim.vcf") as f:
                recs.append({tuple(r[:2]): (r[3], r[4], r[7], r[9:]) for r in
                             (l.rstrip("\n").split("\t") for l in f if not l.startswith("#"))
                             if "SVTYPE=INS" in r[7] or "TYPE=INS" in r[7]})
        alone, mixed = recs
        assert alone, "no insertion reached a genome"
        assert alone == {k: v for k, v in mixed.items() if k in alone}, "an insertion changed"
        assert set(alone) == set(mixed), "the same insertions should be the only ones"


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
