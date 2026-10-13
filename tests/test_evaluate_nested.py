"""
Evaluate on a truth with nested insertions (Simulate's INFO/HOST): a genome without the host
sequence (GT ".") has no site and is scored only at the host's locus, and an insertion into an
insertion may be called as a second ALT of its parent's record.

    PYTHONPATH=. python tests/test_evaluate_nested.py
"""
import os
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
sys.path.insert(0, HERE)

from TEvarSim import evaluate  # noqa: E402
import test_evaluate as te  # noqa: E402

HEADER = te.TRUTH_HEADER_INFO + '##INFO=<ID=HOST,Number=1,Type=String,Description="Host event">\n'
KID = "GATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACA"
SAMPLES = ["S0", "S1", "S2"]


def _evaluate(d, truth_records, pred_records):
    truth, pred = os.path.join(d, "truth.vcf"), os.path.join(d, "pred.vcf")
    te._write_vcf(truth, SAMPLES, truth_records, HEADER)
    te._write_vcf(pred, SAMPLES, pred_records)
    args = te._Args(truth=truth, pred=pred, predType="VCF", TEtype=None, INSonly=False,
                    nHap=1, max_dist=100, gt_len_tol=50, sample_map=None,
                    size_bins=None, af_bins=None, outprefix=os.path.join(d, "out"))
    return te._quiet(lambda: evaluate.Evaluator(args)._run())


def _scored(ev):
    return {s["id"]: s for s in ev.loci}


def test_an_insertion_into_an_insertion_called_as_a_second_alt():
    parent = te.ANCHOR + te.ELEMENT
    k = len(te.ELEMENT) // 2
    child = te.ANCHOR + te.ELEMENT[:k] + KID + te.ELEMENT[k:]
    truth = [(1000, "TY1-FULL#LTR/Copia_", te.ANCHOR, [parent], "TYPE=INS;EVENTTYPE=INS", ["1", "1", "0"]),
             (1000, "kid#LINE/L1_", te.ANCHOR, [child],
              "TYPE=INS;EVENTTYPE=INS;HOST=chrT:1000:TY1-FULL#LTR/Copia_", ["0", "1", "."])]
    pred = [(1000, "call.1", te.ANCHOR, [parent, child], ".", ["1", "2", "0"])]
    with tempfile.TemporaryDirectory() as d:
        ev = _evaluate(d, truth, pred)
        s = _scored(ev)
        p, c = s["TY1-FULL#LTR/Copia_"], s["kid#LINE/L1_"]
        assert p["detected"] and c["detected"], (p["match"], c["match"])
        assert (p["carriers"]["tp"], p["carriers"]["fp"], p["carriers"]["fn"]) == (2, 0, 0), p["carriers"]
        assert (c["carriers"]["tp"], c["carriers"]["fp"], c["carriers"]["fn"]) == (1, 0, 0), c["carriers"]
        assert c["host"] == "TY1-FULL#LTR/Copia_" and c["no_site"] == ["S2"], c
        assert c["alleles"]["noncarrier"]["total"] == 1      # S0, which has the parent but not the child


def test_a_genome_without_the_host_has_no_site_at_the_nested_locus():
    host_ref = te.ANCHOR + te.ELEMENT * 3
    truth = [(1000, "DEL-chrT-1000-1600-LTR/Copia-TY1", host_ref, [te.ANCHOR], "TYPE=DEL;EVENTTYPE=DEL",
              ["1", "0", "0"]),
             (1300, "kid#LINE/L1_", "A", ["A" + KID],
              "TYPE=INS;EVENTTYPE=INS;HOST=chrT:1000:DEL-chrT-1000-1600-LTR/Copia-TY1", [".", "1", "0"])]
    # The call writes S0 as 0 at the nested site, which is no error: S0 has no such site.
    pred = [(1000, "call.1", host_ref, [te.ANCHOR], ".", ["1", "0", "0"]),
            (1300, "call.2", "A", ["A" + KID], ".", ["0", "1", "0"])]
    with tempfile.TemporaryDirectory() as d:
        s = _scored(_evaluate(d, truth, pred))
        n = s["kid#LINE/L1_"]
        assert n["no_site"] == ["S0"], n
        assert (n["carriers"]["tp"], n["carriers"]["fp"], n["carriers"]["fn"]) == (1, 0, 0), n["carriers"]
        assert n["alleles"]["carrier"]["total"] + n["alleles"]["noncarrier"]["total"] == 2
        h = s["DEL-chrT-1000-1600-LTR/Copia-TY1"]
        assert (h["carriers"]["tp"], h["carriers"]["fp"], h["carriers"]["fn"]) == (1, 0, 0), h["carriers"]


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
