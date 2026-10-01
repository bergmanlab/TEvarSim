"""
Tests for ``tevarsim Evaluate``.

Evaluate scores a prediction against every simulated event at once: the unit of analysis is
the event, not the genome, and the simulated genomes are paired with the prediction's samples
by name rather than named on the command line. The tests below pin the behaviour that makes
that trustworthy -- the one-to-one event/prediction matching, the carrier and genotype
arithmetic over all paired genomes, the stratification, and the refusal to guess a pairing.

No external test runner is required::

    PYTHONPATH=. python tests/test_evaluate.py

The test functions are also discoverable by pytest.
"""
import json
import os
import sys
import tempfile

# Import the working-tree package regardless of any installed copy.
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from TEvarSim import evaluate  # noqa: E402


class _Args:
    def __init__(self, **kw):
        self.__dict__.update(kw)


# ---- shared fixture sequences ------------------------------------------------
ANCHOR = "A"
ELEMENT = "GATTACA" * 847          # 5929 bp, stands in for a full-length Ty1
ELEMENT_JITTERED = ELEMENT[:5923]  # the same element seen 6 bp shorter across a TSD
SOLO_LTR = ELEMENT[:336]           # what an excision leaves behind
STACKED = ELEMENT + ELEMENT        # two elements at one site

TRUTH_HEADER_INFO = (
    '##INFO=<ID=TYPE,Number=1,Type=String,Description="Variant type">\n'
    '##INFO=<ID=EVENTTYPE,Number=A,Type=String,Description="Event per ALT allele">\n'
    '##INFO=<ID=MEPRESENT,Number=A,Type=Integer,Description="Element present per ALT">\n'
    '##INFO=<ID=LTRLEN,Number=1,Type=Integer,Description="Solo LTR length">\n'
)


def _write_vcf(path, samples, records, info_headers=""):
    """records: list of (pos, varID, ref, [alts], info, [gt per sample])"""
    with open(path, "w") as f:
        f.write("##fileformat=VCFv4.2\n")
        f.write("##contig=<ID=chrT,length=1000000>\n")
        f.write(info_headers)
        f.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
                + "\t".join(samples) + "\n")
        for pos, varID, ref, alts, info, gts in records:
            f.write("\t".join(["chrT", str(pos), varID, ref, ",".join(alts), ".", "PASS",
                               info or ".", "GT"] + list(gts)) + "\n")


def _evaluate(d, truth_samples, truth_records, pred_samples, pred_records, **kw):
    """Run one Evaluate over freshly written fixtures and return the Evaluator."""
    truth = os.path.join(d, "truth.vcf")
    pred = os.path.join(d, "pred.vcf")
    _write_vcf(truth, truth_samples, truth_records, TRUTH_HEADER_INFO)
    if kw.get("predType") == "BED":
        pred = os.path.join(d, "pred.bed")
        with open(pred, "w") as f:
            for chrom, start, end, name in pred_records:
                f.write(f"{chrom}\t{start}\t{end}\t{name}\n")
    else:
        _write_vcf(pred, pred_samples, pred_records)
    args = _Args(truth=truth, pred=pred, predType="VCF", TEtype=None, INSonly=False,
                 nHap=1, max_dist=100, gt_len_tol=50, sample_map=None,
                 size_bins=None, af_bins=None, outprefix=os.path.join(d, "out"))
    args.__dict__.update(kw)
    args.pred = pred
    return evaluate.Evaluator(args)._run()


def _ins(pos, gts, varID="TY1-FULL#LTR/Copia_1INDEL", alt=None):
    return (pos, varID, ANCHOR, [alt or (ANCHOR + ELEMENT)], "TYPE=INS;EVENTTYPE=INS", gts)


def _quiet(fn, *a, **kw):
    """Run something with its printed report swallowed."""
    import contextlib
    import io
    with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
        return fn(*a, **kw)


# ---- detection ---------------------------------------------------------------


def test_every_event_is_scored_without_naming_a_genome():
    """Three events over three genomes, evaluated in one pass with no --truthID."""
    samples = ["S0", "S1", "S2"]
    with tempfile.TemporaryDirectory() as d:
        ev = _quiet(_evaluate, d, samples,
                    [_ins(1000, ["1", "1", "0"]),
                     _ins(5000, ["0", "1", "1"]),
                     _ins(9000, ["1", "0", "0"])],
                    samples,
                    [(1002, "1.1", ANCHOR, [ANCHOR + ELEMENT_JITTERED], ".", ["1", "1", "0"]),
                     (5000, "1.2", ANCHOR, [ANCHOR + ELEMENT], ".", ["0", "1", "1"])])
        overall = ev.summary["overall"]
        assert overall["n_loci"] == 3, overall
        assert overall["n_detected"] == 2, overall
        assert overall["recovery_rate"] == round(2 / 3, 4), overall
        # the event nobody called contributes its carrier to FN, not to a genotype error
        assert overall["carriers"] == {"tp": 4, "fp": 0, "fn": 1,
                                       "recall": 0.8, "precision": 1.0, "f1": 0.8889}, \
            overall["carriers"]
        assert overall["genotypes"]["compared"] == 6, overall["genotypes"]
        assert overall["genotypes"]["concordance"] == 1.0, overall["genotypes"]
    print("PASS test_every_event_is_scored_without_naming_a_genome")


def test_breakpoint_and_length_error_are_measured_per_event():
    with tempfile.TemporaryDirectory() as d:
        ev = _quiet(_evaluate, d, ["S0"], [_ins(1000, ["1"])], ["S0"],
                    [(1007, "1.1", ANCHOR, [ANCHOR + ELEMENT_JITTERED], ".", ["1"])])
        match = ev.loci[0]["match"]
        assert match["pos_offset"] == 7, match
        assert match["length_error"] == len(ELEMENT_JITTERED) - len(ELEMENT), match
        assert match["allele_match"] is True, match
    print("PASS test_breakpoint_and_length_error_are_measured_per_event")


def test_offsets_are_signed_so_a_bias_does_not_cancel_against_jitter():
    """A call short of the simulated value reads negative, past it positive.

    Absolute values would give a caller consistently 4bp downstream and one scattering
    +-4bp at random the same mean, and those are different problems.
    """
    with tempfile.TemporaryDirectory() as d:
        # both calls land 4bp downstream: a bias, so the mean is +4 and the SD 0
        ev = _quiet(_evaluate, d, ["S0"], [_ins(1000, ["1"]), _ins(9000, ["1"])], ["S0"],
                    [(1004, "1.1", ANCHOR, [ANCHOR + ELEMENT], ".", ["1"]),
                     (9004, "1.2", ANCHOR, [ANCHOR + ELEMENT], ".", ["1"])])
        offsets = ev.summary["overall"]["breakpoint_offset_bp"]
        assert (offsets["n"], offsets["mean"], offsets["sd"]) == (2, 4.0, 0.0), offsets
    with tempfile.TemporaryDirectory() as d:
        # one short, one long by the same amount: no bias, but the scatter is visible
        ev = _quiet(_evaluate, d, ["S0"], [_ins(1000, ["1"]), _ins(9000, ["1"])], ["S0"],
                    [(996, "1.1", ANCHOR, [ANCHOR + ELEMENT], ".", ["1"]),
                     (9004, "1.2", ANCHOR, [ANCHOR + ELEMENT], ".", ["1"])])
        offsets = ev.summary["overall"]["breakpoint_offset_bp"]
        assert (offsets["mean"], offsets["sd"]) == (0.0, 4.0), offsets
    print("PASS test_offsets_are_signed_so_a_bias_does_not_cancel_against_jitter")


def test_a_short_allele_reads_negative_and_a_long_one_positive():
    with tempfile.TemporaryDirectory() as d:
        ev = _quiet(_evaluate, d, ["S0"], [_ins(1000, ["1"])], ["S0"],
                    [(1000, "1.1", ANCHOR, [ANCHOR + ELEMENT_JITTERED], ".", ["1"])])
        errors = ev.summary["overall"]["allele_length_error_bp"]
        # the prediction is 6bp shorter than what was simulated
        assert errors["mean"] == float(len(ELEMENT_JITTERED) - len(ELEMENT)), errors
        assert errors["mean"] < 0, errors
    print("PASS test_a_short_allele_reads_negative_and_a_long_one_positive")


def test_a_prediction_beyond_max_dist_is_not_a_detection():
    with tempfile.TemporaryDirectory() as d:
        ev = _quiet(_evaluate, d, ["S0"], [_ins(1000, ["1"])], ["S0"],
                    [(1101, "1.1", ANCHOR, [ANCHOR + ELEMENT], ".", ["1"])])
        assert ev.summary["overall"]["n_detected"] == 0, ev.summary["overall"]
        assert ev.summary["predictions"]["unmatched"] == 1, ev.summary["predictions"]
        assert ev.summary["predictions"]["precision"] == 0.0, ev.summary["predictions"]
    print("PASS test_a_prediction_beyond_max_dist_is_not_a_detection")


def test_one_prediction_cannot_cover_two_stacked_events():
    """Two elements at one locus need two calls; one call recovers one of them."""
    with tempfile.TemporaryDirectory() as d:
        ev = _quiet(_evaluate, d, ["S0"],
                    [_ins(1000, ["1"]), _ins(1010, ["1"])],
                    ["S0"],
                    [(1000, "1.1", ANCHOR, [ANCHOR + ELEMENT], ".", ["1"])])
        assert ev.summary["overall"]["n_detected"] == 1, ev.summary["overall"]
        # the prediction is claimed by the event it sits on, not by both
        detected = [e for e in ev.loci if e["detected"]]
        assert detected[0]["pos"] == 1000, detected[0]
    print("PASS test_one_prediction_cannot_cover_two_stacked_events")


def test_allele_agreement_wins_over_proximity_when_pairing():
    """A nearer call of the wrong allele must not steal the event from the right one."""
    with tempfile.TemporaryDirectory() as d:
        ev = _quiet(_evaluate, d, ["S0"], [_ins(1000, ["1"])], ["S0"],
                    [(1001, "solo", ANCHOR, [ANCHOR + SOLO_LTR], ".", ["1"]),
                     (1040, "elem", ANCHOR, [ANCHOR + ELEMENT_JITTERED], ".", ["1"])])
        match = ev.loci[0]["match"]
        assert match["id"] == "elem", match
        assert match["allele_match"] is True, match
    print("PASS test_allele_agreement_wins_over_proximity_when_pairing")


def test_a_matched_call_of_the_wrong_allele_is_detected_but_not_concordant():
    """A solo LTR called where a full element was simulated: right locus, wrong allele."""
    with tempfile.TemporaryDirectory() as d:
        ev = _quiet(_evaluate, d, ["S0"], [_ins(1000, ["1"])], ["S0"],
                    [(1000, "solo", ANCHOR, [ANCHOR + SOLO_LTR], ".", ["1"])])
        overall = ev.summary["overall"]
        assert overall["n_detected"] == 1, overall
        assert overall["n_allele_concordant"] == 0, overall
        assert overall["genotypes"]["concordance"] == 0.0, overall["genotypes"]
        # the locus was found, so the carrier is still a carrier TP
        assert overall["carriers"]["tp"] == 1, overall["carriers"]
    print("PASS test_a_matched_call_of_the_wrong_allele_is_detected_but_not_concordant")


def test_a_stacked_call_is_not_the_same_allele_as_a_single_element():
    with tempfile.TemporaryDirectory() as d:
        ev = _quiet(_evaluate, d, ["S0"], [_ins(1000, ["1"])], ["S0"],
                    [(1000, "1.1", ANCHOR, [ANCHOR + STACKED], ".", ["1"])])
        assert ev.summary["overall"]["n_allele_concordant"] == 0, ev.summary["overall"]
    print("PASS test_a_stacked_call_is_not_the_same_allele_as_a_single_element")


# ---- carriers and genotypes --------------------------------------------------


def test_a_genome_called_as_a_carrier_that_is_not_one_is_a_carrier_fp():
    samples = ["S0", "S1"]
    with tempfile.TemporaryDirectory() as d:
        ev = _quiet(_evaluate, d, samples, [_ins(1000, ["1", "0"])], samples,
                    [(1000, "1.1", ANCHOR, [ANCHOR + ELEMENT], ".", ["1", "1"])])
        carriers = ev.summary["overall"]["carriers"]
        assert (carriers["tp"], carriers["fp"], carriers["fn"]) == (1, 1, 0), carriers
        assert ev.loci[0]["carriers"]["predicted"] == ["S0", "S1"], ev.loci[0]["carriers"]
        assert ev.loci[0]["carriers"]["truth"] == ["S0"], ev.loci[0]["carriers"]
    print("PASS test_a_genome_called_as_a_carrier_that_is_not_one_is_a_carrier_fp")


def test_the_no_variant_allele_does_not_make_a_genome_a_carrier():
    """Some callers write <*> for "no variant in this genome"; it is not a call."""
    samples = ["S0", "S1"]
    with tempfile.TemporaryDirectory() as d:
        ev = _quiet(_evaluate, d, samples, [_ins(1000, ["1", "0"])], samples,
                    [(1000, "1.1", ANCHOR, ["<*>", ANCHOR + ELEMENT], ".", ["2", "1"])])
        carriers = ev.summary["overall"]["carriers"]
        assert (carriers["tp"], carriers["fp"], carriers["fn"]) == (1, 0, 0), carriers
    print("PASS test_the_no_variant_allele_does_not_make_a_genome_a_carrier")


def test_a_prediction_record_no_genome_calls_is_not_counted_against_precision():
    with tempfile.TemporaryDirectory() as d:
        ev = _quiet(_evaluate, d, ["S0"], [_ins(1000, ["1"])], ["S0"],
                    [(1000, "1.1", ANCHOR, [ANCHOR + ELEMENT], ".", ["1"]),
                     (7000, "1.2", ANCHOR, ["<*>"], ".", ["1"])])
        assert ev.summary["predictions"]["total"] == 1, ev.summary["predictions"]
        assert ev.summary["predictions"]["precision"] == 1.0, ev.summary["predictions"]
    print("PASS test_a_prediction_record_no_genome_calls_is_not_counted_against_precision")


def test_a_haploid_truth_is_comparable_to_a_diploid_call():
    with tempfile.TemporaryDirectory() as d:
        ev = _quiet(_evaluate, d, ["S0"], [_ins(1000, ["1"])], ["S0"],
                    [(1000, "1.1", ANCHOR, [ANCHOR + ELEMENT_JITTERED], ".", ["1/1"])])
        assert ev.summary["overall"]["genotypes"]["concordance"] == 1.0, ev.summary["overall"]
        assert ev.meta["genotype_label"] == "genotype", ev.meta["genotype_label"]
    print("PASS test_a_haploid_truth_is_comparable_to_a_diploid_call")


def test_haploid_data_is_reported_as_haplotype_concordance():
    with tempfile.TemporaryDirectory() as d:
        ev = _quiet(_evaluate, d, ["S0"], [_ins(1000, ["1"])], ["S0"],
                    [(1000, "1.1", ANCHOR, [ANCHOR + ELEMENT], ".", ["1"])])
        assert ev.meta["genotype_label"] == "haplotype", ev.meta["genotype_label"]
    print("PASS test_haploid_data_is_reported_as_haplotype_concordance")


# ---- sample pairing ----------------------------------------------------------


def test_samples_are_paired_by_name():
    truth_samples = ["A", "B", "C"]
    with tempfile.TemporaryDirectory() as d:
        ev = _quiet(_evaluate, d, truth_samples, [_ins(1000, ["1", "0", "1"])],
                    ["B", "C"], [(1000, "1.1", ANCHOR, [ANCHOR + ELEMENT], ".", ["0", "1"])])
        assert ev.meta["sample_pairs"] == [("B", "B"), ("C", "C")], ev.meta["sample_pairs"]
        # A is unpaired, so its carrier status is out of scope: only B and C are scored
        carriers = ev.summary["overall"]["carriers"]
        assert (carriers["tp"], carriers["fp"], carriers["fn"]) == (1, 0, 0), carriers
        # allele frequency still comes from every simulated genome, paired or not
        assert ev.loci[0]["allele_frequency"] == round(2 / 3, 4), ev.loci[0]
    print("PASS test_samples_are_paired_by_name")


def test_a_lone_sample_on_each_side_is_paired_whatever_it_is_called():
    with tempfile.TemporaryDirectory() as d:
        ev = _quiet(_evaluate, d, ["Hap1"], [_ins(1000, ["1"])], ["callers_name"],
                    [(1000, "1.1", ANCHOR, [ANCHOR + ELEMENT], ".", ["1"])])
        assert ev.meta["sample_pairs"] == [("Hap1", "callers_name")], ev.meta["sample_pairs"]
    print("PASS test_a_lone_sample_on_each_side_is_paired_whatever_it_is_called")


def test_unpairable_samples_still_report_detection_but_no_carriers():
    """Guessing a pairing by column order would report confident nonsense; refuse instead."""
    with tempfile.TemporaryDirectory() as d:
        ev = _quiet(_evaluate, d, ["A", "B"], [_ins(1000, ["1", "1"])], ["X", "Y"],
                    [(1000, "1.1", ANCHOR, [ANCHOR + ELEMENT], ".", ["1", "1"])])
        assert ev.meta["sample_pairs"] == [], ev.meta["sample_pairs"]
        assert ev.summary["overall"]["n_detected"] == 1, ev.summary["overall"]
        carriers = ev.summary["overall"]["carriers"]
        assert (carriers["tp"], carriers["fp"], carriers["fn"]) == (0, 0, 0), carriers
    print("PASS test_unpairable_samples_still_report_detection_but_no_carriers")


def test_sample_map_pairs_differently_named_genomes():
    with tempfile.TemporaryDirectory() as d:
        mapping = os.path.join(d, "map.tsv")
        with open(mapping, "w") as f:
            f.write("# truth\tpred\nA\tX\nB\tY\n")
        ev = _quiet(_evaluate, d, ["A", "B"], [_ins(1000, ["1", "0"])], ["X", "Y"],
                    [(1000, "1.1", ANCHOR, [ANCHOR + ELEMENT], ".", ["1", "0"])],
                    sample_map=mapping)
        assert ev.meta["sample_pairs"] == [("A", "X"), ("B", "Y")], ev.meta["sample_pairs"]
        carriers = ev.summary["overall"]["carriers"]
        assert (carriers["tp"], carriers["fp"], carriers["fn"]) == (1, 0, 0), carriers
    print("PASS test_sample_map_pairs_differently_named_genomes")


def test_sample_map_naming_an_absent_sample_is_an_error():
    with tempfile.TemporaryDirectory() as d:
        mapping = os.path.join(d, "map.tsv")
        with open(mapping, "w") as f:
            f.write("A\tNOPE\n")
        try:
            _quiet(_evaluate, d, ["A"], [_ins(1000, ["1"])], ["X"],
                   [(1000, "1.1", ANCHOR, [ANCHOR + ELEMENT], ".", ["1"])],
                   sample_map=mapping)
        except ValueError as e:
            assert "NOPE" in str(e), str(e)
        else:
            raise AssertionError("expected ValueError for an unknown prediction sample")
    print("PASS test_sample_map_naming_an_absent_sample_is_an_error")


def test_nhap_merges_haplotype_columns_before_pairing():
    """A diploid simulation writes one column per haplotype; --nHap 2 makes individuals."""
    with tempfile.TemporaryDirectory() as d:
        ev = _quiet(_evaluate, d, ["Hap1", "Hap2"], [_ins(1000, ["1", "0"])],
                    ["Hap1_Hap2"],
                    [(1000, "1.1", ANCHOR, [ANCHOR + ELEMENT], ".", ["1|0"])],
                    nHap=2)
        assert ev.meta["sample_pairs"] == [("Hap1_Hap2", "Hap1_Hap2")], ev.meta["sample_pairs"]
        assert ev.meta["n_truth_samples"] == 1, ev.meta
        assert ev.summary["overall"]["genotypes"]["concordance"] == 1.0, ev.summary["overall"]
    print("PASS test_nhap_merges_haplotype_columns_before_pairing")


# ---- filtering and stratification --------------------------------------------


def test_stratification_splits_events_by_type_family_size_and_frequency():
    samples = ["S0", "S1", "S2", "S3"]
    exc = (2000, "EXC-chrT-2000-7929-336-LTR/Copia-TY1-FULL", ANCHOR + ELEMENT,
           [ANCHOR + SOLO_LTR], "TYPE=EXC;EVENTTYPE=EXC", ["1", "0", "0", "0"])
    with tempfile.TemporaryDirectory() as d:
        ev = _quiet(_evaluate, d, samples,
                    [_ins(1000, ["1", "1", "1", "1"]), exc,
                     _ins(9000, ["1", "0", "0", "0"], varID="TY2-FULL#LTR/Copia_2SNP")],
                    samples,
                    [(1000, "1.1", ANCHOR, [ANCHOR + ELEMENT], ".", ["1", "1", "1", "1"]),
                     (2000, "1.2", ANCHOR + ELEMENT, [ANCHOR + SOLO_LTR], ".",
                      ["1", "0", "0", "0"])])
        summary = ev.summary
        assert set(summary["by_event_type"]) == {"INS", "EXC"}, list(summary["by_event_type"])
        assert summary["by_event_type"]["INS"]["n_loci"] == 2, summary["by_event_type"]
        assert set(summary["by_family"]) == {"TY1-FULL", "TY2-FULL"}, list(summary["by_family"])
        # the excision shortens the locus by ~5.6 kb, the insertions lengthen it by ~5.9 kb
        assert set(summary["by_size"]) == {"5kb-10kb"}, list(summary["by_size"])
        assert summary["by_size"]["5kb-10kb"]["n_loci"] == 3, summary["by_size"]
        assert summary["by_allele_frequency"]["0.75-1"]["n_loci"] == 1, \
            summary["by_allele_frequency"]
        assert summary["by_carrier_count"]["1"]["n_loci"] == 2, summary["by_carrier_count"]
        assert summary["by_carrier_count"]["4"]["n_loci"] == 1, summary["by_carrier_count"]
    print("PASS test_stratification_splits_events_by_type_family_size_and_frequency")


# A merged record: the element was inserted, then excised to a solo LTR on a later allele.
# INFO/TYPE is Number=1 and names only the event that created the record, so it reads INS.
_MERGED = (3000, "TY1-FULL#LTR/Copia_1INDEL", ANCHOR,
           [ANCHOR + ELEMENT, ANCHOR + SOLO_LTR], "TYPE=INS;EVENTTYPE=INS,EXC")
# A standalone excision of an element the reference already carried.
_STANDALONE_EXC = (5000, "EXC-chrT-5000-10929-336-LTR/Copia-TY1-FULL", ANCHOR + ELEMENT,
                   [ANCHOR + SOLO_LTR], "TYPE=EXC;EVENTTYPE=EXC")


def test_an_excised_insertion_counts_under_both_ins_and_exc():
    """INFO/TYPE names only the event that created the record, so counting by it put every
    excision of a simulated element in the INS row and left the EXC row counting only the
    excisions of elements the reference already carried."""
    with tempfile.TemporaryDirectory() as d:
        ev = _quiet(_evaluate, d, ["S0"],
                    [_ins(1000, ["1"]), _MERGED + (["2"],), _STANDALONE_EXC + (["1"],)],
                    ["S0"], [])
        merged = [e for e in ev.loci if e["pos"] == 3000][0]
        assert merged["type"] == "INS", merged["type"]
        assert merged["event_classes"] == ["EXC", "INS"], merged["event_classes"]
        by_type = ev.summary["by_event_type"]
        # 3 events, but 4 rows-worth: the merged one is in both
        assert by_type["INS"]["n_loci"] == 2, by_type["INS"]
        assert by_type["EXC"]["n_loci"] == 2, by_type["EXC"]
        assert sum(s["n_loci"] for s in by_type.values()) == 4, by_type
        assert ev.summary["overall"]["n_loci"] == 3, ev.summary["overall"]
        # the history table still counts each record once
        assert sum(s["n_loci"] for s in ev.summary["by_event_history"].values()) == 3
    print("PASS test_an_excised_insertion_counts_under_both_ins_and_exc")


def test_an_event_is_only_counted_where_its_allele_has_the_right_shape():
    """An INS allele has to add sequence and an EXC allele has to take it away."""
    # declares INS,EXC but the second allele is a stacked pair -- it grows, not shrinks
    bad = (3000, "TY1-FULL#LTR/Copia_1INDEL", ANCHOR,
           [ANCHOR + ELEMENT, ANCHOR + STACKED], "TYPE=INS;EVENTTYPE=INS,EXC", ["2"])
    with tempfile.TemporaryDirectory() as d:
        ev = _quiet(_evaluate, d, ["S0"], [bad], ["S0"], [])
        event = ev.loci[0]
        assert event["event_classes"] == ["INS"], event["event_classes"]
        assert event["unsupported_events"], event
        assert "EXC allele 2" in event["unsupported_events"][0], event["unsupported_events"]
        assert "EXC" not in ev.summary["by_event_type"], list(ev.summary["by_event_type"])
        assert len(ev.summary["unsupported_events"]) == 1, ev.summary["unsupported_events"]
    print("PASS test_an_event_is_only_counted_where_its_allele_has_the_right_shape")


def test_an_unsupported_event_is_named_in_the_report():
    import contextlib
    import io
    bad = (3000, "TY1-FULL#LTR/Copia_1INDEL", ANCHOR,
           [ANCHOR + ELEMENT, ANCHOR + STACKED], "TYPE=INS;EVENTTYPE=INS,EXC", ["2"])
    with tempfile.TemporaryDirectory() as d:
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf), contextlib.redirect_stderr(io.StringIO()):
            _evaluate(d, ["S0"], [bad], ["S0"], [])
        out = buf.getvalue()
        assert "no allele bears out" in out, out
        assert "chrT:3000" in out, out
    print("PASS test_an_unsupported_event_is_named_in_the_report")


def test_a_length_change_within_tolerance_counts_as_neither():
    """Below --gt_len_tol two alleles are the same length everywhere else, so a change that
    small is no evidence of an insertion or an excision."""
    tiny = (3000, "TY1-FULL#LTR/Copia_1INDEL", ANCHOR, [ANCHOR + "A" * 10],
            "TYPE=INS;EVENTTYPE=INS", ["1"])
    with tempfile.TemporaryDirectory() as d:
        ev = _quiet(_evaluate, d, ["S0"], [tiny], ["S0"], [], gt_len_tol=50)
        assert ev.loci[0]["event_classes"] == [], ev.loci[0]["event_classes"]
        assert ev.loci[0]["unsupported_events"], ev.loci[0]
        # it still appears, under a label that says it demonstrated nothing
        assert list(ev.summary["by_event_type"]) == ["?"], list(ev.summary["by_event_type"])
    print("PASS test_a_length_change_within_tolerance_counts_as_neither")


# A standalone excision: REF holds the full element, the ALT the anchor plus the solo LTR.
# INFO/LTRLEN says how long that solo LTR is, which is how far downstream a caller anchoring
# past it would put the call.
_LTRLEN = len(SOLO_LTR)
_EXC = (5000, "EXC-chrT-5000-10929-336-LTR/Copia-TY1-FULL", ANCHOR + ELEMENT,
        [ANCHOR + SOLO_LTR], f"TYPE=EXC;EVENTTYPE=EXC;LTRLEN={_LTRLEN}")


def test_an_excision_called_past_its_solo_ltr_counts_as_found_but_not_as_correct():
    """A caller anchoring where the alleles start to differ puts an excision at the far end
    of the surviving LTR, one whole LTR from the simulated POS. The locus was found, so it
    is not a miss -- but the call is in the wrong place, so it is not a correct one."""
    with tempfile.TemporaryDirectory() as d:
        ev = _quiet(_evaluate, d, ["S0"], [_EXC + (["1"],)], ["S0"],
                    [(5000 + _LTRLEN, "6.1", ANCHOR + ELEMENT, [ANCHOR + SOLO_LTR], ".", ["1"])])
        overall = ev.summary["overall"]
        # counts toward recovery ...
        assert (overall["n_recovered"], overall["recovery_rate"]) == (1, 1.0), overall
        assert (overall["n_detected"], overall["n_displaced"]) == (0, 1), overall
        # ... and toward nothing else
        assert ev.summary["predictions"]["unmatched"] == 1, ev.summary["predictions"]
        assert ev.summary["predictions"]["displaced"] == 1, ev.summary["predictions"]
        assert ev.summary["predictions"]["precision"] == 0.0, ev.summary["predictions"]
        assert overall["carriers"]["tp"] == 0, overall["carriers"]
        assert overall["carriers"]["fn"] == 1, overall["carriers"]
        assert overall["breakpoint_offset_bp"]["n"] == 0, overall["breakpoint_offset_bp"]
        assert overall["n_allele_concordant"] == 0, overall
        locus = ev.loci[0]
        assert locus["recovered"] is True and locus["detected"] is False, locus
        assert locus["displaced"]["pos_offset"] == _LTRLEN, locus["displaced"]
    print("PASS test_an_excision_called_past_its_solo_ltr_counts_as_found_but_not_as_correct")


def test_a_correctly_placed_excision_call_still_matches_normally():
    with tempfile.TemporaryDirectory() as d:
        ev = _quiet(_evaluate, d, ["S0"], [_EXC + (["1"],)], ["S0"],
                    [(5000, "6.1", ANCHOR + ELEMENT, [ANCHOR + SOLO_LTR], ".", ["1"])])
        overall = ev.summary["overall"]
        assert (overall["n_detected"], overall["n_displaced"]) == (1, 0), overall
        assert ev.summary["predictions"]["precision"] == 1.0, ev.summary["predictions"]
        assert overall["carriers"]["tp"] == 1, overall["carriers"]
    print("PASS test_a_correctly_placed_excision_call_still_matches_normally")


def test_a_displaced_call_is_recorded_as_such_in_its_locus_file():
    with tempfile.TemporaryDirectory() as d:
        ev = _quiet(_evaluate, d, ["S0"], [_EXC + (["1"],)], ["S0"],
                    [(5000 + _LTRLEN, "6.1", ANCHOR + ELEMENT, [ANCHOR + SOLO_LTR], ".", ["1"])])
        payloads = [json.load(open(os.path.join(ev.locus_dir, n)))
                    for n in sorted(os.listdir(ev.locus_dir))]
        call = [p for p in payloads if p["kind"] == "unmatched_prediction"][0]
        assert call["displaced_match_for"]["pos"] == 5000, call["displaced_match_for"]
        assert call["displaced_match_for"]["pos_offset"] == _LTRLEN, call["displaced_match_for"]
    print("PASS test_a_displaced_call_is_recorded_as_such_in_its_locus_file")


def test_the_ltr_shift_falls_back_to_the_excision_allele_length():
    """Without INFO/LTRLEN the excision allele is the anchor plus the solo LTR, so its own
    length says how far a displaced call would sit."""
    from TEvarSim.evaluate import ltr_shift

    class _L:
        alt_events = ["EXC"]
        alts = (ANCHOR + SOLO_LTR,)
        info = {}
    assert ltr_shift(_L()) == len(SOLO_LTR), ltr_shift(_L())
    _L.info = {"LTRLEN": 300}
    assert ltr_shift(_L()) == 300, ltr_shift(_L())
    # an insertion-only locus has no solo LTR and so no shift to try
    class _I:
        alt_events = ["INS"]
        alts = (ANCHOR + ELEMENT,)
        info = {}
    assert ltr_shift(_I()) is None, ltr_shift(_I())
    print("PASS test_the_ltr_shift_falls_back_to_the_excision_allele_length")


def test_a_far_off_call_is_not_rescued_as_displaced():
    """The shift is one solo LTR, not a licence to match anything downstream."""
    with tempfile.TemporaryDirectory() as d:
        ev = _quiet(_evaluate, d, ["S0"], [_EXC + (["1"],)], ["S0"],
                    [(5000 + _LTRLEN + 5000, "6.1", ANCHOR + ELEMENT, [ANCHOR + SOLO_LTR],
                      ".", ["1"])])
        overall = ev.summary["overall"]
        assert (overall["n_recovered"], overall["n_displaced"]) == (0, 0), overall
        assert ev.summary["predictions"]["displaced"] == 0, ev.summary["predictions"]
    print("PASS test_a_far_off_call_is_not_rescued_as_displaced")


def test_the_superfamily_table_appears_only_when_it_groups_something():
    """TY1 and TY2 are both LTR/Copia, so the superfamily is a coarser cut; one family alone
    would make the table a copy of the family table, so it is left out."""
    ty2 = _ins(9000, ["1"], varID="TY2-FULL#LTR/Copia_2SNP")
    with tempfile.TemporaryDirectory() as d:
        ev = _quiet(_evaluate, d, ["S0"], [_ins(1000, ["1"]), ty2], ["S0"], [])
        assert list(ev.summary["by_superfamily"]) == ["Copia"], ev.summary["by_superfamily"]
    with tempfile.TemporaryDirectory() as d:
        ev = _quiet(_evaluate, d, ["S0"], [_ins(1000, ["1"])], ["S0"], [])
        assert "by_superfamily" not in ev.summary, list(ev.summary)
    print("PASS test_the_superfamily_table_appears_only_when_it_groups_something")


def test_size_bins_are_reported_in_bin_order_not_by_size_of_group():
    with tempfile.TemporaryDirectory() as d:
        ev = _quiet(_evaluate, d, ["S0"],
                    [_ins(1000, ["1"], alt=ANCHOR + SOLO_LTR),
                     _ins(5000, ["1"]),
                     _ins(9000, ["1"]),
                     _ins(20000, ["1"], alt=ANCHOR + STACKED)],
                    ["S0"], [])
        assert list(ev.summary["by_size"]) == ["100bp-500bp", "5kb-10kb", ">=10kb"], \
            list(ev.summary["by_size"])
    print("PASS test_size_bins_are_reported_in_bin_order_not_by_size_of_group")


def test_custom_bin_edges_are_honoured():
    with tempfile.TemporaryDirectory() as d:
        ev = _quiet(_evaluate, d, ["S0", "S1"],
                    [_ins(1000, ["1", "0"]), _ins(5000, ["1", "1"])], ["S0", "S1"], [],
                    size_bins=[1000], af_bins=[0.5])
        assert list(ev.summary["by_size"]) == [">=1kb"], list(ev.summary["by_size"])
        assert set(ev.summary["by_allele_frequency"]) == {"0-0.5", "0.5-1"}, \
            list(ev.summary["by_allele_frequency"])
    print("PASS test_custom_bin_edges_are_honoured")


def test_event_history_separates_an_excised_insertion_from_a_plain_one():
    """A lineage truth VCF merges an element's later history onto its record."""
    history = (1000, "TY1-FULL#LTR/Copia_1INDEL", ANCHOR,
               [ANCHOR + ELEMENT, ANCHOR + SOLO_LTR],
               "TYPE=INS;EVENTTYPE=INS,EXC;MEPRESENT=1,1", ["1", "2"])
    with tempfile.TemporaryDirectory() as d:
        ev = _quiet(_evaluate, d, ["S0", "S1"], [history, _ins(9000, ["1", "1"])],
                    ["S0", "S1"], [])
        assert ev.loci[0]["history"] == "INS,EXC", ev.loci[0]
        assert ev.loci[1]["history"] == "INS", ev.loci[1]
        assert set(ev.summary["by_event_history"]) == {"INS", "INS,EXC"}, \
            list(ev.summary["by_event_history"])
        # INFO/TYPE reads INS for both records, but the excision allele on the first is
        # still counted as an excision rather than being absorbed into the INS row
        assert set(ev.summary["by_event_type"]) == {"INS", "EXC"}, \
            list(ev.summary["by_event_type"])
        assert ev.summary["by_event_type"]["INS"]["n_loci"] == 2, ev.summary["by_event_type"]
        assert ev.summary["by_event_type"]["EXC"]["n_loci"] == 1, ev.summary["by_event_type"]
    print("PASS test_event_history_separates_an_excised_insertion_from_a_plain_one")


def test_a_nested_element_is_labelled_as_such():
    nested = (1000, "TY1-FULL#LTR/Copia_1INDEL", ANCHOR,
              [ANCHOR + ELEMENT, ANCHOR + STACKED],
              "TYPE=INS;EVENTTYPE=INS,INS;MEPRESENT=0,1", ["1", "2"])
    with tempfile.TemporaryDirectory() as d:
        ev = _quiet(_evaluate, d, ["S0", "S1"], [nested], ["S0", "S1"], [])
        assert ev.loci[0]["history"] == "nested INS", ev.loci[0]
    print("PASS test_a_nested_element_is_labelled_as_such")


def test_insonly_and_tetype_filter_the_events_that_are_scored():
    exc = (2000, "EXC-chrT-2000-7929-336-LTR/Copia-TY1-FULL", ANCHOR + ELEMENT,
           [ANCHOR + SOLO_LTR], "TYPE=EXC;EVENTTYPE=EXC", ["1"])
    ty2 = _ins(9000, ["1"], varID="TY2-FULL#LTR/Copia_2SNP")
    with tempfile.TemporaryDirectory() as d:
        ev = _quiet(_evaluate, d, ["S0"], [_ins(1000, ["1"]), exc, ty2], ["S0"], [],
                    INSonly=True)
        assert ev.summary["overall"]["n_loci"] == 2, ev.summary["overall"]
    with tempfile.TemporaryDirectory() as d:
        ev = _quiet(_evaluate, d, ["S0"], [_ins(1000, ["1"]), exc, ty2], ["S0"], [],
                    TEtype="ty1-full")
        assert ev.summary["overall"]["n_loci"] == 2, ev.summary["overall"]
        assert list(ev.summary["by_family"]) == ["TY1-FULL"], list(ev.summary["by_family"])
    print("PASS test_insonly_and_tetype_filter_the_events_that_are_scored")


def test_unknown_tetype_raises_with_the_available_families():
    with tempfile.TemporaryDirectory() as d:
        try:
            _quiet(_evaluate, d, ["S0"], [_ins(1000, ["1"])], ["S0"], [], TEtype="Copia")
        except ValueError as e:
            assert "TY1-FULL" in str(e), str(e)
        else:
            raise AssertionError("expected ValueError for an unknown --TEtype")
    print("PASS test_unknown_tetype_raises_with_the_available_families")


# ---- output ------------------------------------------------------------------


def test_json_holds_one_object_per_event_and_is_reloadable():
    samples = ["S0", "S1"]
    with tempfile.TemporaryDirectory() as d:
        _quiet(_evaluate, d, samples, [_ins(1000, ["1", "0"]), _ins(9000, ["0", "1"])],
               samples, [(1000, "1.1", ANCHOR, [ANCHOR + ELEMENT], ".", ["1", "0"])])
        with open(os.path.join(d, "out.json")) as f:
            report = json.load(f)
        assert set(report) == {"meta", "summary", "loci"}, list(report)
        assert len(report["loci"]) == 2, report["loci"]
        first = report["loci"][0]
        assert first["pos"] == 1000 and first["detected"] is True, first
        assert first["carriers"]["truth"] == ["S0"], first["carriers"]
        assert report["loci"][1]["match"] is None, report["loci"][1]
        assert report["meta"]["pairing"] == "matching names", report["meta"]
    print("PASS test_json_holds_one_object_per_event_and_is_reloadable")


def test_the_printed_report_names_the_headline_numbers():
    import contextlib
    import io
    samples = ["S0", "S1"]
    with tempfile.TemporaryDirectory() as d:
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf), contextlib.redirect_stderr(io.StringIO()):
            _evaluate(d, samples, [_ins(1000, ["1", "1"])], samples,
                      [(1000, "1.1", ANCHOR, [ANCHOR + ELEMENT], ".", ["1", "1"])])
        out = buf.getvalue()
        for expected in ("Locus detection", "recall", "By event type", "By TE family",
                         "By event size", "By allele frequency", "haplotype concordance"):
            assert expected in out, (expected, out)
    print("PASS test_the_printed_report_names_the_headline_numbers")


def test_one_json_is_written_per_locus_either_side_knows_about():
    samples = ["S0", "S1"]
    with tempfile.TemporaryDirectory() as d:
        ev = _quiet(_evaluate, d, samples,
                    [_ins(1000, ["1", "0"]), _ins(9000, ["0", "1"])], samples,
                    [(1000, "1.1", ANCHOR, [ANCHOR + ELEMENT], ".", ["1", "0"]),
                     (40000, "9.9", ANCHOR, [ANCHOR + ELEMENT], ".", ["1", "1"])])
        written = sorted(os.listdir(ev.locus_dir))
        # two simulated events (one recovered, one missed) and one unmatched prediction
        assert written == ["chrT_1000.json", "chrT_40000.json", "chrT_9000.json"], written
        kinds = {}
        for name in written:
            with open(os.path.join(ev.locus_dir, name)) as f:
                kinds[name] = json.load(f)["kind"]
        assert kinds["chrT_1000.json"] == "simulated_locus", kinds
        assert kinds["chrT_9000.json"] == "simulated_locus", kinds
        assert kinds["chrT_40000.json"] == "unmatched_prediction", kinds
    print("PASS test_one_json_is_written_per_locus_either_side_knows_about")


def test_a_locus_file_is_its_summary_entry_plus_run_context():
    with tempfile.TemporaryDirectory() as d:
        ev = _quiet(_evaluate, d, ["S0"], [_ins(1000, ["1"])], ["S0"],
                    [(1003, "1.1", ANCHOR, [ANCHOR + ELEMENT], ".", ["1"])])
        with open(os.path.join(d, "out.json")) as f:
            summary = json.load(f)
        entry = summary["loci"][0]
        assert entry["locus_file"] == "out_loci/chrT_1000.json", entry["locus_file"]
        with open(os.path.join(d, entry["locus_file"])) as f:
            locus = json.load(f)
        assert locus.pop("kind") == "simulated_locus"
        run = locus.pop("run")
        assert locus == entry, (locus, entry)
        # the file stands alone: it says which run produced it and where the rest lives
        assert run["max_dist"] == 100 and run["n_paired_genomes"] == 1, run
        assert run["summary_file"].endswith("out.json"), run
    print("PASS test_a_locus_file_is_its_summary_entry_plus_run_context")


def test_stacked_events_at_one_position_get_distinct_files():
    with tempfile.TemporaryDirectory() as d:
        ev = _quiet(_evaluate, d, ["S0"],
                    [_ins(1000, ["1"]), _ins(1000, ["1"])], ["S0"],
                    [(1000, "1.1", ANCHOR, [ANCHOR + ELEMENT], ".", ["1"])])
        written = sorted(os.listdir(ev.locus_dir))
        assert written == ["chrT_1000-2.json", "chrT_1000.json"], written
        # one of the two was recovered and the other was not; neither overwrote the other
        detected = sorted(json.load(open(os.path.join(ev.locus_dir, n)))["detected"]
                          for n in written)
        assert detected == [False, True], detected
    print("PASS test_stacked_events_at_one_position_get_distinct_files")


def test_an_unmatched_prediction_records_its_nearest_simulated_locus():
    with tempfile.TemporaryDirectory() as d:
        ev = _quiet(_evaluate, d, ["S0"], [_ins(1000, ["1"])], ["S0"],
                    [(1400, "1.1", ANCHOR, [ANCHOR + ELEMENT], ".", ["1"])])
        with open(os.path.join(ev.locus_dir, "chrT_1400.json")) as f:
            locus = json.load(f)
        assert locus["kind"] == "unmatched_prediction", locus
        assert locus["id"] == "1.1" and locus["carriers"]["predicted"] == ["S0"], locus
        near = locus["nearest_simulated_locus"]
        # 400 bp away: outside --max_dist, but plainly a near miss rather than a novel call
        assert near["distance_bp"] == 400 and near["pos"] == 1000, near
        assert near["detected"] is False, near
    print("PASS test_an_unmatched_prediction_records_its_nearest_simulated_locus")


def test_rerunning_clears_locus_files_left_by_the_previous_run():
    """A file from an earlier truth VCF would read as this run's result."""
    with tempfile.TemporaryDirectory() as d:
        ev = _quiet(_evaluate, d, ["S0"], [_ins(1000, ["1"]), _ins(9000, ["1"])],
                    ["S0"], [])
        assert len(os.listdir(ev.locus_dir)) == 2, os.listdir(ev.locus_dir)
        ev = _quiet(_evaluate, d, ["S0"], [_ins(1000, ["1"])], ["S0"], [])
        assert os.listdir(ev.locus_dir) == ["chrT_1000.json"], os.listdir(ev.locus_dir)
    print("PASS test_rerunning_clears_locus_files_left_by_the_previous_run")


def test_a_contig_name_unsafe_in_a_filename_is_sanitised():
    from TEvarSim.evaluate import locus_filenames

    class _Rec:
        chrom, pos = "scaffold|7/a", 42
    events = [{"chrom": "chr:1*", "pos": 100}]
    assert locus_filenames(events, [_Rec()]) == ["chr_1__100.json", "scaffold_7_a_42.json"]
    print("PASS test_a_contig_name_unsafe_in_a_filename_is_sanitised")


def test_a_bed_prediction_scores_detection_only():
    with tempfile.TemporaryDirectory() as d:
        ev = _quiet(_evaluate, d, ["S0"], [_ins(1000, ["1"])], [],
                    [("chrT", 1003, 1004, "call1")], predType="BED")
        assert ev.summary["overall"]["n_detected"] == 1, ev.summary["overall"]
        assert ev.loci[0]["match"]["pos_offset"] == 3, ev.loci[0]["match"]
        # no allele sequences and no samples in a BED: nothing else is claimed
        assert ev.loci[0]["match"]["allele_match"] is None, ev.loci[0]["match"]
        assert ev.loci[0]["match"]["length_error"] is None, ev.loci[0]["match"]
        assert ev.meta["sample_pairs"] == [], ev.meta["sample_pairs"]
    print("PASS test_a_bed_prediction_scores_detection_only")


def test_an_empty_prediction_scores_zero_without_dividing_by_zero():
    with tempfile.TemporaryDirectory() as d:
        ev = _quiet(_evaluate, d, ["S0"], [_ins(1000, ["1"])], ["S0"], [])
        overall = ev.summary["overall"]
        assert overall["recovery_rate"] == 0.0, overall
        assert overall["allele_concordance_rate"] is None, overall
        assert overall["breakpoint_offset_bp"]["mean"] is None, overall
        assert ev.summary["predictions"]["precision"] is None, ev.summary["predictions"]
        assert overall["carriers"]["recall"] == 0.0, overall["carriers"]
    print("PASS test_an_empty_prediction_scores_zero_without_dividing_by_zero")


# ---- counts by locus and by allele -------------------------------------------


def _counts(ev, key):
    c = ev.summary["overall"]["counts"][key]
    return (c["total"], c["genotyped"], c["tp"], c["fp"], c["fn"])


def test_each_allele_class_is_scored_as_the_positive_class_in_turn():
    """A noncarrier called a carrier is a carrier FP and a noncarrier FN; a carrier called a
    noncarrier the other way round. "all" sums the classes."""
    samples = ["S0", "S1", "S2", "S3"]
    with tempfile.TemporaryDirectory() as d:
        ev = _quiet(_evaluate, d, samples, [_ins(1000, ["1", "1", "0", "0"])], samples,
                    [(1000, "1.1", ANCHOR, [ANCHOR + ELEMENT], ".", ["1", "0", "1", "0"])])
        assert _counts(ev, "carrier") == (2, 2, 1, 1, 1), _counts(ev, "carrier")
        assert _counts(ev, "noncarrier") == (2, 2, 1, 1, 1), _counts(ev, "noncarrier")
        assert _counts(ev, "all") == (4, 4, 2, 2, 2), _counts(ev, "all")
        c = ev.summary["overall"]["counts"]["carrier"]
        assert (c["recall"], c["precision"], c["f1"]) == (0.5, 0.5, 0.5), c
    print("PASS test_each_allele_class_is_scored_as_the_positive_class_in_turn")


def test_a_missing_genotype_is_an_fn_of_its_class_and_an_fp_of_neither():
    """A "." is no call: before this, a carrier genotyped "." read as a noncarrier call."""
    samples = ["S0", "S1", "S2"]
    with tempfile.TemporaryDirectory() as d:
        ev = _quiet(_evaluate, d, samples, [_ins(1000, ["1", "1", "0"])], samples,
                    [(1000, "1.1", ANCHOR, [ANCHOR + ELEMENT], ".", ["1", ".", "."])])
        assert _counts(ev, "carrier") == (2, 1, 1, 0, 1), _counts(ev, "carrier")
        assert _counts(ev, "noncarrier") == (1, 0, 0, 0, 1), _counts(ev, "noncarrier")
        assert ev.loci[0]["alleles"]["untyped"] == ["S1", "S2"], ev.loci[0]["alleles"]
    print("PASS test_a_missing_genotype_is_an_fn_of_its_class_and_an_fp_of_neither")


def test_a_missed_locus_leaves_every_allele_ungenotyped():
    samples = ["S0", "S1"]
    with tempfile.TemporaryDirectory() as d:
        ev = _quiet(_evaluate, d, samples, [_ins(1000, ["1", "0"])], samples,
                    [(9000, "2.1", ANCHOR, [ANCHOR + ELEMENT], ".", ["0", "1"])])
        assert _counts(ev, "loci") == (1, 0, 0, 1, 1), _counts(ev, "loci")
        assert _counts(ev, "carrier") == (1, 0, 0, 1, 1), _counts(ev, "carrier")
        assert _counts(ev, "noncarrier") == (1, 0, 0, 0, 1), _counts(ev, "noncarrier")
    print("PASS test_a_missed_locus_leaves_every_allele_ungenotyped")


def test_the_carriers_of_a_call_on_no_simulated_locus_are_carrier_fps():
    """Evaluate's per-locus carrier FP cannot see them; the run's counts do."""
    samples = ["S0", "S1"]
    with tempfile.TemporaryDirectory() as d:
        ev = _quiet(_evaluate, d, samples, [_ins(1000, ["1", "0"])], samples,
                    [(1000, "1.1", ANCHOR, [ANCHOR + ELEMENT], ".", ["1", "0"]),
                     (9000, "2.1", ANCHOR, [ANCHOR + ELEMENT], ".", ["1", "1"])])
        assert _counts(ev, "loci") == (1, 1, 1, 1, 0), _counts(ev, "loci")
        assert _counts(ev, "carrier") == (1, 1, 1, 2, 0), _counts(ev, "carrier")
        assert ev.summary["overall"]["carriers"]["fp"] == 0, "the per-locus carrier block is unchanged"
        # a stratum has no locus FPs to count
        row = ev.summary["by_event_type"]["INS"]["counts"]["loci"]
        assert row["fp"] is None and row["precision"] is None, row
    print("PASS test_the_carriers_of_a_call_on_no_simulated_locus_are_carrier_fps")


def test_a_displaced_call_is_recalled_but_genotypes_nothing_and_is_no_fp():
    with tempfile.TemporaryDirectory() as d:
        ev = _quiet(_evaluate, d, ["S0"], [_EXC + (["1"],)], ["S0"],
                    [(5000 + _LTRLEN, "6.1", ANCHOR + ELEMENT, [ANCHOR + SOLO_LTR], ".", ["1"])])
        assert _counts(ev, "loci") == (1, 1, 1, 0, 0), _counts(ev, "loci")
        assert _counts(ev, "carrier") == (1, 0, 0, 0, 1), _counts(ev, "carrier")
    print("PASS test_a_displaced_call_is_recalled_but_genotypes_nothing_and_is_no_fp")


if __name__ == "__main__":
    test_every_event_is_scored_without_naming_a_genome()
    test_breakpoint_and_length_error_are_measured_per_event()
    test_offsets_are_signed_so_a_bias_does_not_cancel_against_jitter()
    test_a_short_allele_reads_negative_and_a_long_one_positive()
    test_a_prediction_beyond_max_dist_is_not_a_detection()
    test_one_prediction_cannot_cover_two_stacked_events()
    test_allele_agreement_wins_over_proximity_when_pairing()
    test_a_matched_call_of_the_wrong_allele_is_detected_but_not_concordant()
    test_a_stacked_call_is_not_the_same_allele_as_a_single_element()
    test_a_genome_called_as_a_carrier_that_is_not_one_is_a_carrier_fp()
    test_the_no_variant_allele_does_not_make_a_genome_a_carrier()
    test_a_prediction_record_no_genome_calls_is_not_counted_against_precision()
    test_a_haploid_truth_is_comparable_to_a_diploid_call()
    test_haploid_data_is_reported_as_haplotype_concordance()
    test_samples_are_paired_by_name()
    test_a_lone_sample_on_each_side_is_paired_whatever_it_is_called()
    test_unpairable_samples_still_report_detection_but_no_carriers()
    test_sample_map_pairs_differently_named_genomes()
    test_sample_map_naming_an_absent_sample_is_an_error()
    test_nhap_merges_haplotype_columns_before_pairing()
    test_stratification_splits_events_by_type_family_size_and_frequency()
    test_an_excised_insertion_counts_under_both_ins_and_exc()
    test_an_event_is_only_counted_where_its_allele_has_the_right_shape()
    test_an_unsupported_event_is_named_in_the_report()
    test_a_length_change_within_tolerance_counts_as_neither()
    test_an_excision_called_past_its_solo_ltr_counts_as_found_but_not_as_correct()
    test_a_correctly_placed_excision_call_still_matches_normally()
    test_a_displaced_call_is_recorded_as_such_in_its_locus_file()
    test_the_ltr_shift_falls_back_to_the_excision_allele_length()
    test_a_far_off_call_is_not_rescued_as_displaced()
    test_the_superfamily_table_appears_only_when_it_groups_something()
    test_size_bins_are_reported_in_bin_order_not_by_size_of_group()
    test_custom_bin_edges_are_honoured()
    test_event_history_separates_an_excised_insertion_from_a_plain_one()
    test_a_nested_element_is_labelled_as_such()
    test_insonly_and_tetype_filter_the_events_that_are_scored()
    test_unknown_tetype_raises_with_the_available_families()
    test_json_holds_one_object_per_event_and_is_reloadable()
    test_the_printed_report_names_the_headline_numbers()
    test_one_json_is_written_per_locus_either_side_knows_about()
    test_a_locus_file_is_its_summary_entry_plus_run_context()
    test_stacked_events_at_one_position_get_distinct_files()
    test_an_unmatched_prediction_records_its_nearest_simulated_locus()
    test_rerunning_clears_locus_files_left_by_the_previous_run()
    test_a_contig_name_unsafe_in_a_filename_is_sanitised()
    test_a_bed_prediction_scores_detection_only()
    test_an_empty_prediction_scores_zero_without_dividing_by_zero()
    test_each_allele_class_is_scored_as_the_positive_class_in_turn()
    test_a_missing_genotype_is_an_fn_of_its_class_and_an_fp_of_neither()
    test_a_missed_locus_leaves_every_allele_ungenotyped()
    test_the_carriers_of_a_call_on_no_simulated_locus_are_carrier_fps()
    test_a_displaced_call_is_recalled_but_genotypes_nothing_and_is_no_fp()
