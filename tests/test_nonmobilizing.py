"""
TErandom --nNM, Simulate's non-mobilizing VCF and Evaluate --nonmobilizing: deletions that cut into or
swallow a reference element by no TE mechanism, in the genomes, out of the TE truth, and charged to a
caller that calls them.

    PYTHONPATH=. python tests/test_nonmobilizing.py
"""
import json
import os
import random
import subprocess
import sys
import tempfile

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ELEN = 2000


def _seq(n, seed):
    rnd = random.Random(seed)
    return "".join(rnd.choice("ACGT") for _ in range(n))


def _tevarsim(*args):
    run = subprocess.run([sys.executable, "-m", "TEvarSim", *map(str, args)], cwd=ROOT,
                         capture_output=True, text=True, env={**os.environ, "PYTHONPATH": ROOT})
    assert run.returncode == 0, run.stderr[-3000:]
    return run


def _fasta(path):
    seqs, name = {}, None
    for line in open(path):
        if line.startswith(">"):
            name = line[1:].split()[0]
            seqs[name] = []
        else:
            seqs[name].append(line.strip())
    return {k: "".join(v) for k, v in seqs.items()}


def _vcf(path):
    return [l.rstrip("\n").split("\t") for l in open(path) if not l.startswith("#")]


def test_nonmobilizing_svs():
    with tempfile.TemporaryDirectory() as d:
        ref_seq = _seq(400000, 3)
        ref = os.path.join(d, "ref.fa")
        with open(ref, "w") as f:
            f.write(">chrT\n" + ref_seq + "\n")
        with open(ref + ".fai", "w") as f:
            f.write(f"chrT\t{len(ref_seq)}\t6\t{len(ref_seq)}\t{len(ref_seq) + 1}\n")
        lib = os.path.join(d, "lib.fa")
        with open(lib, "w") as f:
            f.write(f">TY1#LTR/Copia TSD=5\n{_seq(1500, 1)}\n")
        hosts = os.path.join(d, "hosts.bed")
        with open(hosts, "w") as f:
            for start in range(10000, 390000, 10000):
                f.write(f"chrT\t{start}\t{start + ELEN}\tTY1#LTR/Copia\t.\t+\n")
        out = os.path.join(d, "pool")
        _tevarsim("TErandom", "--ref", ref, "--consensus", lib, "--existingTEs", hosts, "--nINS", 4,
                  "--nDEL", 3, "--nNM", 15, "--outprefix", out, "--seed", 4)

        rows = [l.rstrip("\n").split("\t") for l in open(out + ".bed")]
        nm = [r for r in rows if r[3].startswith("NM-")]
        assert len(nm) == 15, nm
        kinds = {}
        for r in nm:
            _, kind, chrom, hs, he, *_ = r[3].split("-")
            s, e, hs, he = int(r[1]), int(r[2]), int(hs), int(he)
            kinds[kind] = kinds.get(kind, 0) + 1
            if kind == "internal":
                assert hs + 50 <= s and e <= he - 50 and 100 <= e - s <= ELEN // 2, r
            elif kind == "left":
                assert hs - 2000 <= s <= hs - 25 and hs + 50 <= e <= he - 50, r
            elif kind == "right":
                assert hs + 50 <= s <= he - 50 and he + 25 <= e <= he + 2000, r
            elif kind == "tight":
                assert hs - 250 <= s <= hs - 25 and he + 25 <= e <= he + 250, r
            else:
                assert kind == "wide" and hs - 5000 <= s <= hs - 250 and he + 250 <= e <= he + 5000, r
        assert kinds == {k: 3 for k in ("internal", "left", "right", "tight", "wide")}, kinds
        spans = sorted((int(r[1]), int(r[2])) for r in rows)
        assert all(a[1] < b[0] or (a[0] == a[1] and a[1] <= b[0]) for a, b in zip(spans, spans[1:])), \
            "events overlap"

        sim = os.path.join(d, "sim")
        _tevarsim("Simulate", "--ref", ref, "--bed", out + ".bed", "--pool", out + ".fa", "--num", 8,
                  "--af-dist", "uniform", "--af-min", 0.5, "--af-max", 0.5, "--tsd-from-header",
                  "--bg-pi", 0.002, "-O", sim, "-D", 5)
        truth = _vcf(sim + ".vcf")
        nmv = _vcf(sim + ".nonmobilizing.vcf")
        assert truth and not any(r[2].startswith("NM-") for r in truth)
        assert nmv and all(r[2].startswith("NM-") and "NMKIND=" in r[7] for r in nmv)
        # Each genome is the reference plus every event, TE, non-mobilizing and background, it carries.
        events = truth + nmv + _vcf(sim + ".background.vcf")
        for g in range(8):
            genome = _fasta(f"{sim}_{g}.fa")[f"chrT_{g}"]
            delta = sum(len(r[4]) - len(r[3]) for r in events if r[9 + g] == "1")
            assert len(genome) == len(ref_seq) + delta, g
            # And a carrier's non-mobilizing deletion is the join of its two flanks.
            for r in nmv:
                if r[9 + g] == "1":
                    pos, end = int(r[1]), int(r[1]) + len(r[3]) - 1
                    assert ref_seq[pos - 30:pos] + ref_seq[end:end + 30] in genome, (g, r[2])

        # A caller that calls every non-mobilizing SV, and one truth event, is charged for each.
        pred = os.path.join(d, "pred.vcf")
        header = [l for l in open(sim + ".vcf") if l.startswith("#")]
        with open(pred, "w") as f:
            f.writelines(header)
            for r in [truth[0]] + nmv:
                f.write("\t".join(r[:7] + ["TYPE=DEL"] + r[8:]) + "\n")
        evp = os.path.join(d, "ev")
        _tevarsim("Evaluate", "--truth", sim + ".vcf", "--pred", pred, "--nonmobilizing",
                  sim + ".nonmobilizing.vcf", "-O", evp)
        summary = json.load(open(evp + ".json"))["summary"]
        got = summary["nonmobilizing"]
        assert got["sites"] == len(nmv) and got["sites_called"] == len(nmv), got
        assert got["calls"] == len(nmv) and summary["predictions"]["matched"] == 1, got


if __name__ == "__main__":
    test_nonmobilizing_svs()
    print("PASS test_nonmobilizing_svs")
