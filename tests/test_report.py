"""
Tests for the HTML report ``tevarsim Evaluate`` writes.

The report exists because the summary cannot show a *shape*: a mean and an SD of the
breakpoint error do not separate a caller that is exact everywhere but badly wrong twice
from one that jitters a base either way on every locus. The tests below pin the things
that make the figures readable rather than merely present -- an exact call getting a bar
of its own, an axis that stays symmetric so a bias cannot hide, the outcome segments
partitioning the simulated distribution instead of double-counting it -- and the two
properties that make the page usable where these simulations actually run: it opens
offline, and it degrades rather than lies when the prediction carries no genotypes.

No external test runner is required::

    PYTHONPATH=. python tests/test_report.py

The test functions are also discoverable by pytest.
"""
import os
import re
import sys
import tempfile

# Import the working-tree package regardless of any installed copy.
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from TEvarSim import report  # noqa: E402

from test_evaluate import (  # noqa: E402
    ANCHOR,
    ELEMENT,
    ELEMENT_JITTERED,
    SOLO_LTR,
    _evaluate,
    _ins,
    _quiet,
)


def _render(ev):
    """The page Evaluate would write for this run."""
    return report.render(ev.meta, ev.summary, ev.loci,
                         truth_rows=ev.truth_rows, pred_rows=ev.pred_rows,
                         locus_rows=ev.locus_rows)


def _viewer_row(page, dom_id):
    """The one viewer row with this DOM id, as raw markup."""
    found = re.search(r'<div class="vcf-row [^"]*" id="%s"[^>]*>' % re.escape(dom_id), page)
    return found.group(0) if found else None


# ---- binning -----------------------------------------------------------------


def test_an_exact_call_gets_a_bar_of_its_own():
    """
    At single-nucleotide resolution the zero bar is the whole point of the figure, so it
    must not be shared with the calls that were a base out or split across two bars.
    """
    centres, width, counts = report.signed_bins([0, 0, 0, 1, -1])
    assert width == 1, width
    assert 0 in centres, centres
    assert counts[centres.index(0)] == 3, list(zip(centres, counts))
    assert counts[centres.index(1)] == 1, list(zip(centres, counts))
    print("PASS test_an_exact_call_gets_a_bar_of_its_own")


def test_the_error_axis_stays_symmetric_so_a_bias_cannot_hide():
    """
    A caller that is only ever late is a different problem from one that jitters. Cropping
    the axis to the data would put both masses in the middle of the plot and hide it.
    """
    centres, _width, counts = report.signed_bins([3, 4, 5, 5, 6])
    assert centres[0] == -centres[-1], centres
    assert sum(counts[:centres.index(0)]) == 0, list(zip(centres, counts))
    print("PASS test_the_error_axis_stays_symmetric_so_a_bias_cannot_hide")


def test_every_value_lands_in_exactly_one_bin():
    values = [-97, -3, 0, 0, 12, 55, 100]
    _centres, _width, counts = report.signed_bins(values)
    assert sum(counts) == len(values), counts
    print("PASS test_every_value_lands_in_exactly_one_bin")


def test_a_handful_of_element_sizes_keep_their_own_bars():
    """
    A simulation drawing from three elements produces three sizes. Binning those to round
    widths scatters real classes over a mostly empty axis, so below the threshold each
    distinct size is its own bar.
    """
    centres, width, key = report.value_bins([-5587, 340, 5930, 5930])
    assert width == 0, width
    assert centres == [-5587, 340, 5930], centres
    assert key(5930) == 5930, key(5930)
    print("PASS test_a_handful_of_element_sizes_keep_their_own_bars")


def test_many_distinct_sizes_are_binned_to_round_widths():
    values = list(range(0, 6000, 37))
    centres, width, key = report.value_bins(values)
    assert width > 0, width
    assert len(centres) <= 30, len(centres)
    # The bins run continuously, so an empty stretch of the distribution stays visible as
    # a gap rather than being compressed out of the axis.
    assert centres == [centres[0] + i * width for i in range(len(centres))], centres
    assert all(abs(key(v) - v) <= width / 2.0 for v in values)
    print("PASS test_many_distinct_sizes_are_binned_to_round_widths")


# ---- the cumulative resolution readout ---------------------------------------


def test_the_resolution_readout_is_cumulative_and_stops_at_max_dist():
    rows = report.resolution_rows([0, 0, 1, 4, 30], max_dist=15)
    counts = [row[1] for row in rows]
    assert counts == sorted(counts), rows
    assert rows[0][0] == "exact (0 bp)" and rows[0][1] == 2, rows[0]
    # Nothing beyond --max_dist could have been matched, so no step past it is offered.
    assert rows[-1][0] == "within 15 bp", rows[-1]
    assert all("within 25 bp" not in row[0] for row in rows), rows
    print("PASS test_the_resolution_readout_is_cumulative_and_stops_at_max_dist")


# ---- the figures, over a real run --------------------------------------------


def test_the_breakpoint_figure_holds_one_bar_per_placed_call():
    samples = ["S0", "S1"]
    with tempfile.TemporaryDirectory() as d:
        ev = _quiet(_evaluate, d, samples,
                    [_ins(1000, ["1", "1"]), _ins(5000, ["1", "0"]),
                     _ins(9000, ["0", "1"])],
                    samples,
                    [(1000, "p1", ANCHOR, [ANCHOR + ELEMENT], ".", ["1", "1"]),
                     (5003, "p2", ANCHOR, [ANCHOR + ELEMENT], ".", ["1", "0"]),
                     (9000, "p3", ANCHOR, [ANCHOR + ELEMENT], ".", ["0", "1"])])
        offsets = [locus["match"]["pos_offset"] for locus in ev.loci if locus["detected"]]
        assert sorted(offsets) == [0, 0, 3], offsets
        svg, rows = report.error_histogram(
            [(offset, f"chrT:{1000 + i}") for i, offset in enumerate(offsets)],
            "breakpoint error", "loci")
        assert sum(row[1] for row in rows) == 3, rows
        assert [row[:2] for row in rows if row[0] == "0"] == [["0", 2]], rows
        assert svg.count("<title>") == 2, svg.count("<title>")
        # every bar names the loci that made it, and names each of them exactly once
        named = [str(ref) for row in rows for ref in row[2]]
        assert sorted(named) == ["chrT:1000", "chrT:1001", "chrT:1002"], rows
        assert [row[1] for row in rows] == [len(row[2]) for row in rows], rows
    print("PASS test_the_breakpoint_figure_holds_one_bar_per_placed_call")


def test_a_displaced_call_is_kept_out_of_the_breakpoint_figure():
    """
    A call anchored past the surviving solo LTR found the locus but not its breakpoint.
    Evaluate gives it detection credit and nothing else; plotting its offset would put a
    whole LTR of error into a figure about single-nucleotide resolution.
    """
    excision = (2000, "TY1-FULL#LTR/Copia_2", ANCHOR,
                [ANCHOR + ELEMENT, ANCHOR + SOLO_LTR],
                "TYPE=INS;EVENTTYPE=INS,EXC;MEPRESENT=1,1;LTRLEN=336", ["2", "2"])
    samples = ["S0", "S1"]
    with tempfile.TemporaryDirectory() as d:
        ev = _quiet(_evaluate, d, samples, [excision], samples,
                    [(2336, "p1", ANCHOR, [ANCHOR + SOLO_LTR], ".", ["1", "1"])])
        assert ev.summary["overall"]["n_displaced"] == 1, ev.summary["overall"]
        page = _render(ev)
        assert "not plotted here" in page, page[:400]
        offsets = [locus["match"]["pos_offset"] for locus in ev.loci if locus["detected"]]
        assert offsets == [], offsets
    print("PASS test_a_displaced_call_is_kept_out_of_the_breakpoint_figure")


def test_the_outcome_segments_partition_the_simulated_distribution():
    """
    The column height is every locus simulated at that size and the segments are a
    partition of it, so "how many were missed here" is read off directly. If a locus could
    land in two segments the column would overstate the simulation.
    """
    samples = ["S0", "S1"]
    with tempfile.TemporaryDirectory() as d:
        ev = _quiet(_evaluate, d, samples,
                    [_ins(1000, ["1", "1"]), _ins(5000, ["1", "0"]),
                     _ins(9000, ["0", "1"])],
                    samples,
                    [(1000, "p1", ANCHOR, [ANCHOR + ELEMENT_JITTERED], ".", ["1", "1"])])
        centres, width, key = report.value_bins(
            [locus["size_bp"] for locus in ev.loci])
        placed, missed = {}, {}
        for locus, label in zip(ev.loci, report.locus_labels(ev.loci)):
            bucket = placed if locus["detected"] else missed
            bucket.setdefault(key(locus["size_bp"]), []).append(label)
        svg, rows = report.outcome_columns(
            centres, width,
            [("correctly placed", "var(--step-near)", placed),
             ("missed", "var(--step-far)", missed)],
            "size", "loci")
        assert svg is not None
        # column 1 is the total; the outcome columns between it and the member list must
        # sum back to it
        for row in rows:
            assert row[1] == sum(row[2:-1]), row
        assert sum(row[1] for row in rows) == len(ev.loci), rows
    print("PASS test_the_outcome_segments_partition_the_simulated_distribution")


def test_the_size_table_names_the_missed_loci_apart_from_the_recovered_ones():
    """
    A bin holding both recovered and missed loci is the case this figure exists to show,
    so its member list is grouped by outcome. A flat list of everything at that size could
    not say which of them were missed, which is the only reason to open the table.
    """
    samples = ["S0", "S1"]
    with tempfile.TemporaryDirectory() as d:
        ev = _quiet(_evaluate, d, samples,
                    [_ins(1000, ["1", "1"]), _ins(5000, ["1", "0"])],
                    samples,
                    [(1000, "p1", ANCHOR, [ANCHOR + ELEMENT], ".", ["1", "1"])])
        page = _render(ev)
        table = re.search(r"Recall by event size.*?</details>", page, re.S).group(0)
        cells = [re.sub(r"<[^>]+>", "", cell)
                 for cell in re.findall(r'<td class=ids>(.*?)</td>', table)]
        assert cells == ["correctly placed: chrT:1000; missed: chrT:5000"], cells
    print("PASS test_the_size_table_names_the_missed_loci_apart_from_the_recovered_ones")


def test_two_elements_stacked_at_one_position_get_distinct_labels():
    """
    Two elements at one site are two loci at one chrom:pos, so the position alone would
    print the same name twice and identify neither. The repeats are suffixed with the same
    numbering Evaluate names their per-locus JSON with, so a label points at one file.
    """
    samples = ["S0"]
    with tempfile.TemporaryDirectory() as d:
        ev = _quiet(_evaluate, d, samples,
                    [_ins(1000, ["1"]), _ins(1000, ["1"]), _ins(5000, ["1"])],
                    samples,
                    [(1000, "p1", ANCHOR, [ANCHOR + ELEMENT], ".", ["1"])])
        labels = [ref.text for ref in report.locus_labels(ev.loci)]
        assert len(set(labels)) == len(labels), labels
        assert labels[:2] == ["chrT:1000", "chrT:1000 #2"], labels
        # the suffix matches the one on that locus's own JSON file
        assert ev.loci[1]["locus_file"].endswith("chrT_1000-2.json"), ev.loci[1]
        page = _render(ev)
        assert "chrT:1000 #2" in page
    print("PASS test_two_elements_stacked_at_one_position_get_distinct_labels")


def test_a_stratification_into_one_group_draws_no_chart():
    """
    One family means one row of bars, which only restates the headline metrics the tiles
    already carry. The table still reports it.
    """
    samples = ["S0", "S1"]
    with tempfile.TemporaryDirectory() as d:
        ev = _quiet(_evaluate, d, samples,
                    [_ins(1000, ["1", "1"]), _ins(5000, ["1", "0"])],
                    samples,
                    [(1000, "p1", ANCHOR, [ANCHOR + ELEMENT], ".", ["1", "1"])])
        assert len(ev.summary["by_family"]) == 1, ev.summary["by_family"]
        assert report.stratum_bars(ev.summary["by_family"],
                                   [("recovery rate", "var(--series-1)",
                                     lambda s: s["recovery_rate"])]) is None
        page = _render(ev)
        assert "By TE family" in page
        assert "TY1-FULL" in page
    print("PASS test_a_stratification_into_one_group_draws_no_chart")


# ---- degrading honestly ------------------------------------------------------


def test_carrier_f1_is_not_drawn_when_no_genome_could_be_paired():
    """
    A BED prediction carries no samples, so the confusion matrix is empty and
    ``calculate_metrics`` reports an F1 of 0 for it. A bar drawn from that would read as
    "every carrier was called wrong" rather than "nothing was compared".
    """
    samples = ["S0", "S1"]
    with tempfile.TemporaryDirectory() as d:
        ev = _quiet(_evaluate, d, samples,
                    [_ins(1000, ["1", "1"]), _ins(5000, ["1", "0"]),
                     _ins(9000, ["0", "1"])],
                    [], [("chrT", 1000, 1001, "b1"), ("chrT", 5000, 5001, "b2")],
                    predType="BED")
        assert ev.meta["sample_pairs"] == [], ev.meta["sample_pairs"]
        assert ev.summary["overall"]["carriers"]["f1"] == 0, ev.summary["overall"]
        page = _render(ev)
        assert "carrier F1</li>" not in page
        assert "genotype concordance</li>" not in page
        # nor a tile claiming a genotype accuracy that was never measured
        assert "Carrier F1</div>" not in page
        # the figure it *can* draw is still there
        assert "Breakpoint resolution" in page
    print("PASS test_carrier_f1_is_not_drawn_when_no_genome_could_be_paired")


def test_a_run_that_recovered_nothing_renders_without_dividing_by_zero():
    samples = ["S0"]
    with tempfile.TemporaryDirectory() as d:
        ev = _quiet(_evaluate, d, samples, [_ins(1000, ["1"])], samples,
                    [(900000, "far", ANCHOR, [ANCHOR + ELEMENT], ".", ["1"])])
        assert ev.summary["overall"]["n_detected"] == 0, ev.summary["overall"]
        page = _render(ev)
        assert "No call was placed on a simulated locus" in page
    print("PASS test_a_run_that_recovered_nothing_renders_without_dividing_by_zero")


# ---- the side-by-side viewer -------------------------------------------------


def test_the_viewer_links_a_matched_pair_from_both_sides():
    """
    The pairing is the thing the viewer exists to show, and it has to be followable in
    either direction: from a simulated locus to the call that found it, and from a call
    back to the locus it was credited against.
    """
    samples = ["S0"]
    with tempfile.TemporaryDirectory() as d:
        ev = _quiet(_evaluate, d, samples, [_ins(1000, ["1"]), _ins(5000, ["1"])], samples,
                    [(5002, "p1", ANCHOR, [ANCHOR + ELEMENT], ".", ["1"])])
        page = _render(ev)
        # the prediction's only record was paired with the second simulated locus
        assert 'id="t1" data-match="p0"' in page, _viewer_row(page, "t1")
        assert 'id="p0" data-match="t1"' in page, _viewer_row(page, "p0")
        # and the locus nothing was paired with points at a sentence, not at a record
        assert 'id="t0" data-note=' in page, _viewer_row(page, "t0")
        assert "this locus was missed" in page
    print("PASS test_the_viewer_links_a_matched_pair_from_both_sides")


def test_the_viewer_says_which_side_has_nothing_to_show():
    """
    "There is no record there" and "the viewer did not do anything" look identical in a
    pane that simply fails to move, so an unpaired record carries the sentence the other
    pane puts up instead of a jump.
    """
    samples = ["S0"]
    with tempfile.TemporaryDirectory() as d:
        ev = _quiet(_evaluate, d, samples, [_ins(1000, ["1"])], samples,
                    [(900000, "spurious", ANCHOR, [ANCHOR + ELEMENT], ".", ["1"])])
        page = _render(ev)
        assert "No prediction was paired with chrT:1000 — this locus was missed." in page
        assert "No simulated locus was paired with the call at chrT:900000." in page
        assert 'data-match=' not in page, "neither record has a partner to jump to"
    print("PASS test_the_viewer_says_which_side_has_nothing_to_show")


def test_a_displaced_call_is_linked_but_marked_as_displaced():
    """
    A call anchored past the solo LTR did find its locus, so the viewer must take you to
    it -- but it is not a correct call, and a row that reads the same as a clean match
    would lose the only distinction Evaluate drew.
    """
    excision = (2000, "TY1-FULL#LTR/Copia_2", ANCHOR,
                [ANCHOR + ELEMENT, ANCHOR + SOLO_LTR],
                "TYPE=INS;EVENTTYPE=INS,EXC;MEPRESENT=1,1;LTRLEN=336", ["2"])
    samples = ["S0"]
    with tempfile.TemporaryDirectory() as d:
        ev = _quiet(_evaluate, d, samples, [excision], samples,
                    [(2336, "p1", ANCHOR, [ANCHOR + SOLO_LTR], ".", ["1"])])
        page = _render(ev)
        assert 'class="vcf-row displ" id="t0" data-match="p0"' in page, _viewer_row(page, "t0")
        assert 'class="vcf-row displ" id="p0" data-match="t0"' in page, _viewer_row(page, "p0")
    print("PASS test_a_displaced_call_is_linked_but_marked_as_displaced")


def test_a_filtered_record_still_gets_a_row_saying_why():
    """
    Dropping the records --TEtype removed would show a smaller simulation than the one that
    ran. They stay, marked unscored, with the filter that removed them named.
    """
    samples = ["S0"]
    other = "TY3-FULL#LTR/Gypsy_9INDEL"
    with tempfile.TemporaryDirectory() as d:
        ev = _quiet(_evaluate, d, samples,
                    [_ins(1000, ["1"]), _ins(5000, ["1"], varID=other)], samples,
                    [(1000, "p1", ANCHOR, [ANCHOR + ELEMENT], ".", ["1"])],
                    TEtype="TY1-FULL")
        assert ev.summary["overall"]["n_loci"] == 1, "only one family was scored"
        page = _render(ev)
        assert len(ev.truth_rows) == 2, ev.truth_rows
        assert 'class="vcf-row skip" id="t1"' in page, _viewer_row(page, "t1")
        assert "dropped by --TEtype TY1-FULL" in page
    print("PASS test_a_filtered_record_still_gets_a_row_saying_why")


def test_the_viewer_shortens_allele_sequences_but_not_the_record():
    """
    A 6 kb ALT is 300x the rest of the line; inlining every one of them whole would put
    megabytes of sequence in a page that has to open off a mounted filesystem. The rest of
    the record -- INFO especially -- is what a reader is here for and is left alone.
    """
    samples = ["S0"]
    with tempfile.TemporaryDirectory() as d:
        ev = _quiet(_evaluate, d, samples, [_ins(1000, ["1"])], samples,
                    [(1000, "p1", ANCHOR, [ANCHOR + ELEMENT], ".", ["1"])])
        page = _render(ev)
        assert ELEMENT not in page, "the whole element must not be inlined"
        assert f"({len(ANCHOR + ELEMENT):,} bp)" in page, "its length is shown instead"
        assert "TYPE=INS;EVENTTYPE=INS" in page, "INFO is untouched"
    print("PASS test_the_viewer_shortens_allele_sequences_but_not_the_record")


def test_every_row_leads_somewhere_and_no_link_dangles():
    """
    A row links to its partner's DOM id, so a partner that never got a row of its own is a
    click that silently does nothing. Every row must therefore carry either a link to a
    record that exists or the sentence to show in its place.
    """
    excision = (2000, "TY1-FULL#LTR/Copia_2", ANCHOR,
                [ANCHOR + ELEMENT, ANCHOR + SOLO_LTR],
                "TYPE=INS;EVENTTYPE=INS,EXC;MEPRESENT=1,1;LTRLEN=336", ["2"])
    samples = ["S0"]
    with tempfile.TemporaryDirectory() as d:
        ev = _quiet(_evaluate, d, samples,
                    [_ins(1000, ["1"]), excision, _ins(9000, ["1"])], samples,
                    [(1000, "p1", ANCHOR, [ANCHOR + ELEMENT], ".", ["1"]),
                     (2336, "p2", ANCHOR, [ANCHOR + SOLO_LTR], ".", ["1"]),
                     (700000, "p3", ANCHOR, [ANCHOR + ELEMENT], ".", ["1"]),
                     (800000, "p4", ANCHOR, [ANCHOR + ELEMENT], ".", ["0"])])
        page = _render(ev)
        rows = re.findall(r'<div class="vcf-row \w+" id="(\w+)"([^>]*)>', page)
        assert rows, page
        ids = {dom_id for dom_id, _rest in rows}
        for dom_id, rest in rows:
            link = re.search(r'data-match="(\w+)"', rest)
            if link:
                assert link.group(1) in ids, f"{dom_id} points at missing {link.group(1)}"
            else:
                assert 'data-note="' in rest, f"{dom_id} leads nowhere and says nothing"
        # the all-reference record is present but unscored, so it says so rather than
        # pretending to be a spurious call
        assert page.count('class="vcf-row skip"') == 1, page
    print("PASS test_every_row_leads_somewhere_and_no_link_dangles")


def test_the_viewer_is_a_tab_reached_from_the_file_names():
    """
    The obvious thing to click when you want to see a file is its name, so that is what
    opens the pane holding it.
    """
    samples = ["S0"]
    with tempfile.TemporaryDirectory() as d:
        ev = _quiet(_evaluate, d, samples, [_ins(1000, ["1"])], samples,
                    [(1000, "p1", ANCHOR, [ANCHOR + ELEMENT], ".", ["1"])])
        page = _render(ev)
        assert '<nav class="tabs" role="tablist" id="tabs" hidden>' in page
        assert 'id="panel-report"' in page and 'id="panel-vcf"' in page
        assert 'data-pane="truth" href="#pane-truth"' in page or \
               'href="#pane-truth" data-pane="truth"' in page, page
        assert 'href="#pane-pred" data-pane="pred"' in page, page
        # the viewer is inside its own panel, not trailing the report
        panel = page.index('id="panel-vcf"')
        assert page.index("Truth and prediction, side by side") > panel
        assert page.index("Breakpoint resolution") < panel
    print("PASS test_the_viewer_is_a_tab_reached_from_the_file_names")


def test_a_locus_in_a_table_links_to_its_record():
    """
    A name in a table and a row in the viewer are the same locus, so the name is the way
    to the record. The link is a real anchor as well, so a page whose tabs never came up
    still lands on the row.
    """
    samples = ["S0"]
    with tempfile.TemporaryDirectory() as d:
        ev = _quiet(_evaluate, d, samples,
                    [_ins(1000, ["1"]), _ins(5000, ["1"])], samples,
                    [(1000, "p1", ANCHOR, [ANCHOR + ELEMENT], ".", ["1"])])
        page = _render(ev)
        link = ('<a class="locus-link" href="#t0" data-locus="t0">chrT:1000</a>')
        assert link in page, page[page.index("breakpoint error (bp)"):][:600]
        # every locus link points at a viewer row that exists
        rows = {dom_id for dom_id, _rest
                in re.findall(r'<div class="vcf-row \w+" id="(\w+)"([^>]*)>', page)}
        for target in set(re.findall(r'data-locus="(\w+)"', page)):
            assert target in rows, target
    print("PASS test_a_locus_in_a_table_links_to_its_record")


def test_locus_names_stay_plain_text_when_there_is_no_viewer():
    """The report is still renderable on its own -- from a summary JSON, say -- and then
    there is nothing for a name to link to."""
    samples = ["S0"]
    with tempfile.TemporaryDirectory() as d:
        ev = _quiet(_evaluate, d, samples, [_ins(1000, ["1"])], samples,
                    [(1000, "p1", ANCHOR, [ANCHOR + ELEMENT], ".", ["1"])])
        page = report.render(ev.meta, ev.summary, ev.loci)
        assert "chrT:1000" in page
        assert 'class="locus-link"' not in page, "no record to link a locus to"
        assert 'class="file-link"' not in page, "no pane to link a file name to"
        assert 'id="panel-vcf"' not in page and 'role="tablist"' not in page
        # and nothing anywhere points at an anchor the page does not hold
        for target in set(re.findall(r'href="#([\w-]+)"', page)):
            assert f'id="{target}"' in page, target
    print("PASS test_locus_names_stay_plain_text_when_there_is_no_viewer")


def test_a_bed_prediction_still_gets_a_pane():
    samples = ["S0"]
    with tempfile.TemporaryDirectory() as d:
        ev = _quiet(_evaluate, d, samples, [_ins(1000, ["1"]), _ins(5000, ["1"])],
                    [], [("chrT", 1000, 1001, "b1")], predType="BED")
        page = _render(ev)
        assert 'id="pane-pred"' in page and 'id="pane-truth"' in page
        assert 'id="p0" data-match="t0"' in page, _viewer_row(page, "p0")
    print("PASS test_a_bed_prediction_still_gets_a_pane")


# ---- the page ----------------------------------------------------------------


def test_the_page_is_self_contained():
    """
    These reports are read wherever the simulation ran -- over a mounted filesystem, or
    copied off a cluster -- so anything fetched at read time is a blank box. Nothing may
    reference an external host or a sibling file.
    """
    samples = ["S0", "S1"]
    with tempfile.TemporaryDirectory() as d:
        ev = _quiet(_evaluate, d, samples,
                    [_ins(1000, ["1", "1"]), _ins(5000, ["1", "0"])],
                    samples,
                    [(1000, "p1", ANCHOR, [ANCHOR + ELEMENT], ".", ["1", "1"])])
        page = _render(ev)
        for pattern in (r"<script[^>]+src=", r"<link\b", r"@import",
                        r"https?://", r"<img\b"):
            assert not re.search(pattern, page), pattern
    print("PASS test_the_page_is_self_contained")


def test_the_page_names_the_headline_numbers():
    samples = ["S0", "S1"]
    with tempfile.TemporaryDirectory() as d:
        ev = _quiet(_evaluate, d, samples,
                    [_ins(1000, ["1", "1"]), _ins(5000, ["1", "0"]),
                     _ins(9000, ["0", "1"])],
                    samples,
                    [(1000, "p1", ANCHOR, [ANCHOR + ELEMENT], ".", ["1", "1"]),
                     (5000, "p2", ANCHOR, [ANCHOR + ELEMENT], ".", ["1", "0"])])
        page = _render(ev)
        assert "66.7%" in page, "the locus recall is the hero figure"
        assert "Breakpoint resolution" in page
        assert "Allele length accuracy" in page
        assert "Recall by event size" in page
        assert "--max_dist 100" in page
        assert "S0" in page, "the paired genomes are named"
    print("PASS test_the_page_names_the_headline_numbers")


def test_evaluate_writes_the_report_beside_its_json_unless_told_not_to():
    samples = ["S0"]
    with tempfile.TemporaryDirectory() as d:
        ev = _quiet(_evaluate, d, samples, [_ins(1000, ["1"])], samples,
                    [(1000, "p1", ANCHOR, [ANCHOR + ELEMENT], ".", ["1"])])
        assert ev.html_file == os.path.join(d, "out.html"), ev.html_file
        assert os.path.isfile(ev.html_file)
    with tempfile.TemporaryDirectory() as d:
        ev = _quiet(_evaluate, d, samples, [_ins(1000, ["1"])], samples,
                    [(1000, "p1", ANCHOR, [ANCHOR + ELEMENT], ".", ["1"])],
                    no_html=True)
        assert ev.html_file is None, ev.html_file
        assert not os.path.exists(os.path.join(d, "out.html"))
    print("PASS test_evaluate_writes_the_report_beside_its_json_unless_told_not_to")


def test_a_locus_id_with_markup_in_it_is_escaped():
    """A TE family name is copied out of the VCF, so it is not trusted as markup."""
    samples = ["S0"]
    nasty = "<script>alert(1)</script>#LTR/Copia_1"
    with tempfile.TemporaryDirectory() as d:
        ev = _quiet(_evaluate, d, samples, [_ins(1000, ["1"], varID=nasty)], samples,
                    [(1000, "p1", ANCHOR, [ANCHOR + ELEMENT], ".", ["1"])])
        page = _render(ev)
        assert "<script>alert(1)</script>" not in page
        assert "&lt;script&gt;alert(1)&lt;/script&gt;" in page
    print("PASS test_a_locus_id_with_markup_in_it_is_escaped")


if __name__ == "__main__":
    test_an_exact_call_gets_a_bar_of_its_own()
    test_the_error_axis_stays_symmetric_so_a_bias_cannot_hide()
    test_every_value_lands_in_exactly_one_bin()
    test_a_handful_of_element_sizes_keep_their_own_bars()
    test_many_distinct_sizes_are_binned_to_round_widths()
    test_the_resolution_readout_is_cumulative_and_stops_at_max_dist()
    test_the_breakpoint_figure_holds_one_bar_per_placed_call()
    test_a_displaced_call_is_kept_out_of_the_breakpoint_figure()
    test_the_outcome_segments_partition_the_simulated_distribution()
    test_the_size_table_names_the_missed_loci_apart_from_the_recovered_ones()
    test_two_elements_stacked_at_one_position_get_distinct_labels()
    test_a_stratification_into_one_group_draws_no_chart()
    test_carrier_f1_is_not_drawn_when_no_genome_could_be_paired()
    test_a_run_that_recovered_nothing_renders_without_dividing_by_zero()
    test_the_viewer_links_a_matched_pair_from_both_sides()
    test_the_viewer_says_which_side_has_nothing_to_show()
    test_a_displaced_call_is_linked_but_marked_as_displaced()
    test_a_filtered_record_still_gets_a_row_saying_why()
    test_the_viewer_shortens_allele_sequences_but_not_the_record()
    test_every_row_leads_somewhere_and_no_link_dangles()
    test_the_viewer_is_a_tab_reached_from_the_file_names()
    test_a_locus_in_a_table_links_to_its_record()
    test_locus_names_stay_plain_text_when_there_is_no_viewer()
    test_a_bed_prediction_still_gets_a_pane()
    test_the_page_is_self_contained()
    test_the_page_names_the_headline_numbers()
    test_evaluate_writes_the_report_beside_its_json_unless_told_not_to()
    test_a_locus_id_with_markup_in_it_is_escaped()
    print("\nAll report tests passed.")
