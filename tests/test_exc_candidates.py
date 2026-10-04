"""
TErandom --excCandidates: excisions from a BED of exactly known full-length LTR elements, and a
deletion candidate overlapping a chosen element is not deleted as well.

    PYTHONPATH=. python tests/test_exc_candidates.py
"""
import os
import random
import subprocess
import sys
import tempfile

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def _seq(n, seed):
    rnd = random.Random(seed)
    return "".join(rnd.choice("ACGT") for _ in range(n))


def test_excisions_from_candidates_and_no_double_use():
    ltr, internal = _seq(300, 1), _seq(1500, 2)
    with tempfile.TemporaryDirectory() as d:
        ref = os.path.join(d, "ref.fa")
        with open(ref, "w") as f:
            f.write(">chrT\n" + _seq(60000, 3) + "\n")
        with open(ref + ".fai", "w") as f:
            f.write("chrT\t60000\t6\t60000\t60001\n")
        lib = os.path.join(d, "lib.fa")
        with open(lib, "w") as f:
            f.write(f">TY1-LTR#LTR/Copia\n{ltr}\n>TY1-I#LTR/Copia\n{internal}\n")
        # Two elements: deletion candidates span element + 5 bp TSD copy, excision candidates the element alone.
        dels, excs = os.path.join(d, "del.bed"), os.path.join(d, "exc.bed")
        with open(dels, "w") as fd, open(excs, "w") as fe:
            for start in (10000, 30000):
                elen = 2 * len(ltr) + len(internal)
                fd.write(f"chrT\t{start}\t{start + elen + 5}\tTY1-FULL#LTR/Copia\t.\t+\n")
                fe.write(f"chrT\t{start}\t{start + elen}\tTY1-FULL#LTR/Copia\t.\t+\t{len(ltr)}\n")
        out = os.path.join(d, "pool")
        run = subprocess.run([sys.executable, "-m", "TEvarSim", "TErandom", "--ref", ref, "--consensus", lib,
                              "--existingTEs", dels, "--excCandidates", excs, "--nINS", "2", "--nDEL", "1",
                              "--nEXC", "1", "--outprefix", out, "--seed", "4"],
                             cwd=ROOT, capture_output=True, text=True, env={**os.environ, "PYTHONPATH": ROOT})
        assert run.returncode == 0, run.stderr[-2000:]
        rows = [l.rstrip("\n").split("\t") for l in open(out + ".bed")]
        exc = [r for r in rows if r[3].startswith("EXC")]
        dele = [r for r in rows if r[3].startswith("DEL")]
        assert len(exc) == 1 and len(exc[0]) == 7 and int(exc[0][6]) == len(ltr), exc
        assert len(dele) == 1, dele
        es, ee, ds, de = int(exc[0][1]), int(exc[0][2]), int(dele[0][1]), int(dele[0][2])
        assert de <= es or ee <= ds, "the excised element was also deleted"


if __name__ == "__main__":
    test_excisions_from_candidates_and_no_double_use()
    print("PASS test_excisions_from_candidates_and_no_double_use")
