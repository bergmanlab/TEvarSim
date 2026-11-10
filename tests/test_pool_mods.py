"""
Poly-A tails and 5' truncation are pool modifications: TErandom gives a pool element a tail
(-A/-M/-X) or cuts its 5' end (-T/-K). Their lengths come from a stream per element, and the
global stream that picks families and modification masks is drawn the same with them on or off,
so a run that adds them keeps every family and position of a run without. --mod-classes limits
them to the named classes.

    PYTHONPATH=. python tests/test_pool_mods.py
"""
import os
import re
import sys
import tempfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from test_layers import _pool  # noqa: E402

SUFFIX = re.compile(r"_(?:\d+(?:SNP|INDEL|polyA|truncate))*$")


def _ins(rows):
    return [(r[0], r[1], r[5], SUFFIX.sub("", r[3])) for r in rows]


def _fasta(path):
    seqs, name = {}, None
    for line in open(path):
        line = line.rstrip("\n")
        if line.startswith(">"):
            name = line[1:].split()[0]
            seqs[name] = []
        else:
            seqs[name].append(line)
    return {k: "".join(v) for k, v in seqs.items()}


def test_modifications_keep_families_and_positions():
    with tempfile.TemporaryDirectory() as d:
        _, off = _pool(d, "off", n_del=5, n_exc=4)
        mods = ("-A", 1, "-T", 1, "-K", 0.5, "-M", 5, "-X", 20)
        out_on, on = _pool(d, "on", n_del=5, n_exc=4, extra=mods)
        out_cls, cls = _pool(d, "cls", n_del=5, n_exc=4, extra=mods + ("--mod-classes", "LINE,SINE"))
        assert _ins(off) == _ins(on) == _ins(cls), "modifications moved an event"
        on_ids = [r[3] for r in on if r[3].endswith(("polyA", "truncate"))]
        assert on_ids, "-A 1 -T 1 modified nothing"
        # With --mod-classes LINE,SINE the LTR elements are untouched, the jockeys modified.
        for r in cls:
            if "#LTR/" in r[3] and int(r[2]) - int(r[1]) <= 1:
                assert r[3].endswith("_"), r[3]
        lines = [r[3] for r in cls if "#LINE/" in r[3]]
        assert lines and all(re.search(r"\d+truncate\d+polyA$", x) for x in lines), lines
        # The same element gets the same tail and cut whatever other classes are allowed.
        on_by_site = {(r[0], r[1]): r[3] for r in on}
        for r in cls:
            if "#LINE/" in r[3]:
                assert on_by_site[(r[0], r[1])] == r[3], (r[3], on_by_site[(r[0], r[1])])
        # A poly-A tail is a 3' run of A; a truncation removes the 5' end.
        seqs = _fasta(out_cls + ".fa")
        base = _fasta(os.path.join(d, "lib.fa"))["jockey#LINE/Jockey"]
        for name, seq in seqs.items():
            m = re.search(r"_(\d+)truncate(\d+)polyA$", name)
            if m:
                cut, tail = int(m.group(1)), int(m.group(2))
                assert seq == base[cut:] + "A" * tail, name


def test_no_modifications_is_byte_identical_to_ratio_zero():
    with tempfile.TemporaryDirectory() as d:
        out_a, a = _pool(d, "a", n_del=5)
        out_b, b = _pool(d, "b", n_del=5, extra=("--mod-classes", "LINE"))
        assert a == b
        assert open(out_a + ".fa").read() == open(out_b + ".fa").read()


def test_simulate_gives_a_modified_insertion_the_same_carriers():
    """Poly-A tails and truncation change an insertion's ALT, never its carriers or TSD."""
    from test_layers import _tevarsim, _fixture
    with tempfile.TemporaryDirectory() as d:
        ref = _fixture(d)[0]
        recs = []
        for name, extra in (("off", ()), ("on", ("-A", 1, "-T", 1, "--mod-classes", "LINE"))):
            out, _ = _pool(d, name, n_del=5, extra=extra + ("--nINS", 60))
            _tevarsim("Simulate", "--ref", ref, "--bed", out + ".bed", "--pool", out + ".fa",
                      "--num", 10, "--af-dist", "exponential", "--af-mean", 0.2, "--af-min", 0.1,
                      "--allow-zero-carriers", "--tsd-from-header", "--seed", 7, "-O", out + "_sim")
            with open(out + "_sim.vcf") as f:
                rows = [l.rstrip("\n").split("\t") for l in f if not l.startswith("#")]
            recs.append({(r[0], r[1], SUFFIX.sub("", r[2])): (r[4], r[9:]) for r in rows})
        off, on = recs
        assert set(off) == set(on), "a modification added or removed a record"
        changed = [k for k in off if off[k][0] != on[k][0]]
        assert changed and all("#LINE/" in k[2] for k in changed), changed
        assert all(off[k][1] == on[k][1] for k in off), "a modification changed carriers"


if __name__ == "__main__":
    for name, fn in list(globals().items()):
        if name.startswith("test_") and callable(fn):
            fn()
            print("ok", name)
