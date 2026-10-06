"""
Simulate --bg-sv-rate and Evaluate --background_sv: background deletions and tandem duplications
shared by descent, in the genomes and their own VCF, out of the TE truth, and charged to a caller
that calls them.

    PYTHONPATH=. python tests/test_background_sv.py
"""
import json
import os
import subprocess
import sys
import tempfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from test_nonmobilizing import _fasta, _seq, _tevarsim, _vcf  # noqa: E402

N_GENOMES = 30
RATE = 40          # per genome per Mb
REF_LEN = 1_000_000


def test_background_svs():
    with tempfile.TemporaryDirectory() as d:
        ref_seq = _seq(REF_LEN, 3)
        ref = os.path.join(d, "ref.fa")
        with open(ref, "w") as f:
            f.write(">chrT\n" + ref_seq + "\n")
        with open(ref + ".fai", "w") as f:
            f.write(f"chrT\t{REF_LEN}\t6\t{REF_LEN}\t{REF_LEN + 1}\n")
        lib = os.path.join(d, "lib.fa")
        with open(lib, "w") as f:
            f.write(f">TY1#LTR/Copia TSD=5\n{_seq(1500, 1)}\n")
        pool = os.path.join(d, "pool")
        _tevarsim("TErandom", "--ref", ref, "--consensus", lib, "--nINS", 20, "--outprefix", pool, "--seed", 4)
        sim = os.path.join(d, "sim")
        _tevarsim("Simulate", "--ref", ref, "--bed", pool + ".bed", "--pool", pool + ".fa", "--num", N_GENOMES,
                  "--af-dist", "uniform", "--af-min", 0.3, "--af-max", 0.3, "--tsd-from-header",
                  "--bg-pi", 0.002, "--bg-sv-rate", RATE, "-O", sim, "-D", 5)
        truth, bg, sv = _vcf(sim + ".vcf"), _vcf(sim + ".background.vcf"), _vcf(sim + ".background_sv.vcf")
        assert sv and not any(r[2].startswith("bgSV") for r in truth + bg)
        types = {r[7].split(";")[0] for r in sv}
        assert types == {"TYPE=DEL", "TYPE=DUP"}, types
        sizes = [abs(int(r[7].split("SVLEN=")[1].split(";")[0])) for r in sv]
        assert 50 <= min(sizes) and max(sizes) <= 10000, (min(sizes), max(sizes))
        # The rate is per genome against the reference; allow for Poisson noise and the genealogies.
        per_genome = sum(r[9:].count("1") for r in sv) / N_GENOMES
        assert 0.4 * RATE <= per_genome <= 1.8 * RATE, per_genome
        # Every genome is the reference plus everything it carries, and each SV is what it says.
        for g in range(N_GENOMES):
            genome = _fasta(f"{sim}_{g}.fa")[f"chrT_{g}"]
            delta = sum(len(r[4]) - len(r[3]) for r in truth + bg + sv if r[9 + g] == "1")
            assert len(genome) == REF_LEN + delta, g
            for r in sv:
                if r[9 + g] != "1":
                    continue
                pos = int(r[1])
                if "TYPE=DUP" in r[7]:
                    lo = int(r[7].split("DUPSTART=")[1]) - 1
                    copy = ref_seq[lo:pos]
                    assert ref_seq[lo - 25:lo] + copy + copy + ref_seq[pos:pos + 25] in genome, (g, r[2])
                else:
                    end = pos + len(r[3]) - 1
                    assert ref_seq[pos - 25:pos] + ref_seq[end:end + 25] in genome, (g, r[2])

        # A caller that calls every background SV and one truth event is charged for each.
        pred = os.path.join(d, "pred.vcf")
        header = [l for l in open(sim + ".vcf") if l.startswith("#")]
        with open(pred, "w") as f:
            f.writelines(header)
            for r in [truth[0]] + sv:
                f.write("\t".join(r[:7] + ["TYPE=SV"] + r[8:]) + "\n")
        evp = os.path.join(d, "ev")
        run = _tevarsim("Evaluate", "--truth", sim + ".vcf", "--pred", pred, "--background_sv",
                        sim + ".background_sv.vcf", "-O", evp)
        assert "Background SVs" in run.stdout
        summary = json.load(open(evp + ".json"))["summary"]
        got = summary["background_sv"]
        # A background SV record can still match a simulated event close enough by position and length;
        # every other one is an unmatched prediction, and every one of those is charged to its SV.
        unmatched = summary["predictions"]["unmatched"]
        assert summary["predictions"]["matched"] >= 1 and unmatched >= len(sv) - 2, summary["predictions"]
        assert got["sites"] == len(sv) and got["calls"] == unmatched and got["sites_called"] == unmatched, \
            (got["sites"], got["sites_called"], got["calls"], unmatched)


if __name__ == "__main__":
    test_background_svs()
    print("PASS test_background_svs")
