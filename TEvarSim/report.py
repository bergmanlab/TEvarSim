'''
``tevarsim Evaluate`` -- the HTML report.

``Evaluate`` prints a fixed-width summary and writes the same numbers to JSON. Both are
faithful and neither shows a *distribution*: "breakpoint offset: mean +0.07 bp, SD 0.60 bp"
cannot distinguish a caller that is exact on every locus but 30bp out on two of them from
one that scatters a bp either side of the truth on all of them. Those are different
problems, and the shape that separates them is a histogram of the per-locus errors -- the
figure the Tangram paper (Wu et al. 2014, BMC Genomics 15:795, Figure 2) reports breakpoint
resolution with, and Figure 1's overlay of the missed events on the simulated size
distribution for detection.

So this module renders the per-locus records ``Evaluate`` already computes as one
self-contained HTML page:

  * the signed breakpoint error of every correctly placed call, as a histogram, plus the
    cumulative "exact / within N bp" readout,
  * the signed allele length error, the same way,
  * the simulated event sizes, with each size's loci stacked by what became of them, so a
    sensitivity floor at a particular element size is visible rather than averaged away,
  * and the existing strata (event type, family, size, allele frequency, carriers) as
    grouped bars beside the tables the text report prints.

Everything is generated here: the charts are inline SVG built from Python, so the page
needs no plotting library at run time and no network at read time. That matters because
these reports are read wherever the simulation ran -- over a mounted filesystem, or copied
off a cluster -- and a chart that fetches a script from a CDN is a blank box there. One
file, opened in any browser, offline.

The palette is the validated default categorical/ordinal set: it clears the colour-vision
separation, chroma and contrast gates in both light and dark mode. Do not substitute
colours here by eye -- an outcome stack in green/amber/red, the obvious choice, fails
deuteranopia separation outright (dE 4.1 against a floor of 8).
'''
import html
import math
import os
from collections import OrderedDict

# The breakdowns, in the order the text report prints them. Shared with ``Evaluate`` so the
# two renderings of the same summary cannot drift apart.
STRATUM_TITLES = (
    ("By event type", "by_event_type"),
    ("By event history", "by_event_history"),
    ("By TE family", "by_family"),
    ("By TE superfamily", "by_superfamily"),
    ("By event size", "by_size"),
    ("By allele frequency", "by_allele_frequency"),
    ("By number of carrier genomes", "by_carrier_count"),
)

# How close a call has to land to count at each step of the resolution readout. Trimmed to
# --max_dist when the report is written, since nothing further out could have matched.
RESOLUTION_STEPS = (0, 1, 2, 5, 10, 15, 25, 50, 100)

# Bin widths a histogram is allowed to use, so the axis reads in round numbers whatever the
# spread of the data. 1bp first: at single-nucleotide resolution one bar per base is the
# whole point of the figure.
_BIN_LADDER = (1, 2, 5, 10, 25, 50, 100, 250, 500, 1000, 2500, 5000, 10000, 25000, 50000)

# Chart geometry, in SVG user units. The page scales the SVG to its column, so these are
# proportions rather than pixels; they are sized so an 11pt axis label is legible at the
# width a report is normally read at.
_CHART_W = 720
_PLOT_H = 220
_PAD_L, _PAD_R, _PAD_T = 56, 14, 14
_AXIS_H = 34          # room under the plot for the tick labels and the axis title
_BAR_MAX = 24         # bars are capped rather than filling their band; the rest is air
_GAP = 2              # the surface gap that separates touching marks


# ---- formatting --------------------------------------------------------------


def _esc(value):
    return html.escape(str(value), quote=True)


def _pct(value, digits=1):
    return "n/a" if value is None else f"{value * 100:.{digits}f}%"


def _signed(value):
    '''
    A signed number with its sign always shown, so the direction of a bias reads at a
    glance. "n/a" rather than "-" for a missing value, which would look like a minus sign,
    and a bare "0" for no error at all, which has no direction to point in.
    '''
    if value is None:
        return "n/a"
    return "0" if value == 0 else f"{value:+,g}"


def _num(value):
    if value is None:
        return "-"
    if isinstance(value, float):
        return f"{value:,g}"
    return f"{value:,}"


def locus_label(locus):
    '''
    How a locus is named in the table behind a figure: ``chrom:pos``, where it was
    simulated.

    Not ``INFO/ID``. A truth VCF straight out of ``Simulate`` gives every insertion of one
    family the same ID -- 25 of the 27 loci in a yeast lineage read ``TY1-FULL#LTR/Copia_``
    -- so a list of those would name nothing. The position is unique to the locus, is what
    ``locus_filenames`` names the per-locus JSON after, and is what the rest of the report
    already identifies a locus by, so a label here leads straight to that locus's file.
    '''
    return f"{locus['chrom']}:{locus['pos']}"


class Html(str):
    '''
    A string that is already markup. ``_table`` inserts these as they are and escapes
    everything else, so a cell can only carry markup by having been built as markup.
    '''


class LocusRef:
    '''
    One locus as it appears in a table: what to call it, and the record it stands for.

    ``anchor`` is the DOM id of that locus's row in the side-by-side viewer, so a name in a
    table is also the way to the record behind it. None when there is no viewer to point
    at, in which case the name renders as plain text.
    '''

    __slots__ = ("text", "anchor")

    def __init__(self, text, anchor=None):
        self.text = text
        self.anchor = anchor

    def __str__(self):
        return self.text

    def html(self):
        if not self.anchor:
            return _esc(self.text)
        return (f'<a class="locus-link" href="#{self.anchor}" '
                f'data-locus="{self.anchor}">{_esc(self.text)}</a>')


def locus_labels(loci, locus_rows=()):
    '''
    Name every locus for the tables, disambiguating the ones that share a position.

    Two elements stacked at one site are two loci at one ``chrom:pos``, so the position
    alone would print the same name two or three times in a row and name none of them. The
    repeats take a ``#2``, ``#3`` suffix in the order they appear -- the same rule, over
    the same ordering, that ``locus_filenames`` uses to name their per-locus JSON ``-2``,
    ``-3``, so a label in a table points at exactly one file on disk.

    ``locus_rows`` is which record of the truth file each locus came from, in the same
    order, which is what lets a name in a table link to that record in the viewer.

    Returns a list of LocusRef parallel to ``loci``.
    '''
    seen, labels = {}, []
    for index, locus in enumerate(loci):
        base = locus_label(locus)
        seen[base] = seen.get(base, 0) + 1
        text = base if seen[base] == 1 else f"{base} #{seen[base]}"
        row = locus_rows[index] if index < len(locus_rows) else None
        labels.append(LocusRef(text, None if row is None else f"t{row}"))
    return labels


def _refs_html(refs):
    '''A run of locus names as links, comma separated.'''
    return Html(", ".join(ref.html() for ref in refs))


def _spread_text(spread, unit="bp"):
    '''Render one of ``Evaluate``'s mean/SD blocks as "+0.07 +/- 0.6 bp".'''
    if not spread or not spread.get("n"):
        return "n/a"
    return f"{_signed(spread['mean'])} ± {_num(spread['sd'])} {unit}"


# ---- binning -----------------------------------------------------------------


def _bin_width(span, max_bins):
    '''The narrowest round bin width that keeps a span inside ``max_bins`` bars.'''
    if span <= 0:
        return _BIN_LADDER[0]
    needed = span / max(1, max_bins - 1)
    for width in _BIN_LADDER:
        if width >= needed:
            return width
    return _BIN_LADDER[-1]


def signed_bins(values, max_bins=41):
    '''
    Bin a signed error into bars centred on zero, symmetrically.

    Centred rather than edge-aligned so that an exact call -- error 0, the value the whole
    figure is about -- gets a bar of its own instead of being split across two or lumped in
    with the calls that were 1bp out. Symmetric about zero so a bias reads as an off-centre
    mass: a caller that is only ever late looks quite different from one that jitters, and
    a plot cropped to the data it happens to hold would hide that.

    Returns (centres, width, counts).
    '''
    span = max((abs(v) for v in values), default=0)
    width = _bin_width(2 * span, max_bins)
    reach = max(1, int(span // width) + 1)
    centres = [k * width for k in range(-reach, reach + 1)]
    counts = [0] * len(centres)
    for value in values:
        counts[int(round(value / width)) + reach] += 1
    return centres, width, counts


def value_bins(values, max_bins=24, max_distinct=14):
    '''
    Lay an axis out over a set of values, binning them only when there are enough distinct
    ones to need it.

    A simulation drawing from a handful of elements produces a handful of sizes -- a solo
    LTR, a full-length element, a stacked pair -- and binning those into round widths
    scatters two real classes over a mostly empty axis. Below ``max_distinct`` values each
    one therefore keeps its own bar; above it they are binned to round widths, with the
    empty bins in between kept so the shape of the distribution is not compressed.

    Returns (centres, width, key), where ``key(value)`` names the centre a value belongs to
    and ``width`` is 0 when the values are their own bins.
    '''
    distinct = sorted(set(values))
    if not distinct:
        return [], 0, (lambda v: v)
    if len(distinct) <= max_distinct:
        return distinct, 0, (lambda v: v)
    lo, hi = distinct[0], distinct[-1]
    width = _bin_width(hi - lo, max_bins)
    first, last = int(round(lo / width)), int(round(hi / width))
    centres = [k * width for k in range(first, last + 1)]
    return centres, width, (lambda v: int(round(v / width)) * width)


def _nice_step(raw):
    '''Round a raw axis step up to 1, 2, 2.5 or 5 times a power of ten.'''
    if raw <= 0:
        return 1
    exponent = math.floor(math.log10(raw))
    base = 10.0 ** exponent
    for multiple in (1, 2, 2.5, 5, 10):
        if raw <= multiple * base:
            step = multiple * base
            break
    return int(step) if step >= 1 else step


def count_ticks(top, target=4):
    '''Tick values for a count axis running 0..top, at a round step.'''
    if top <= 0:
        return [0, 1], 1
    step = max(1, _nice_step(top / max(1, target)))
    ticks, value = [], 0
    while value < top + step:
        ticks.append(value)
        value += step
    return ticks, ticks[-1]


def _label_stride(labels, band):
    '''
    How many bars to skip between x-axis labels so they do not collide. Measured from the
    widest label rather than a fixed rule, because "-5,587" needs four times the room "0"
    does and a stride that works for one overprints the other.
    '''
    widest = max((len(str(label)) for label in labels), default=1) * 6.4 + 8
    return max(1, int(math.ceil(widest / max(band, 1))))


# ---- SVG primitives ----------------------------------------------------------


def _r(value):
    '''Round a coordinate; SVG does not need more than two decimals and they cost bytes.'''
    return f"{value:.2f}".rstrip("0").rstrip(".")


def _column(x, y, w, h, radius=4):
    '''
    A column with rounded data-end and square base: the mark grows from the baseline, so
    only the end that carries the value is rounded.
    '''
    radius = min(radius, w / 2.0, h)
    if radius <= 0:
        return f'<path d="M{_r(x)} {_r(y + h)}h{_r(w)}v{_r(-h)}h{_r(-w)}Z"/>'
    return (f'<path d="M{_r(x)} {_r(y + h)}V{_r(y + radius)}'
            f'A{_r(radius)} {_r(radius)} 0 0 1 {_r(x + radius)} {_r(y)}'
            f'H{_r(x + w - radius)}'
            f'A{_r(radius)} {_r(radius)} 0 0 1 {_r(x + w)} {_r(y + radius)}'
            f'V{_r(y + h)}Z"/>')


def _bar(x, y, w, h, radius=4):
    '''A horizontal bar: rounded at the tip, square where it leaves the baseline.'''
    radius = min(radius, h / 2.0, w)
    if radius <= 0:
        return f'<path d="M{_r(x)} {_r(y)}h{_r(w)}v{_r(h)}h{_r(-w)}Z"/>'
    return (f'<path d="M{_r(x)} {_r(y)}H{_r(x + w - radius)}'
            f'A{_r(radius)} {_r(radius)} 0 0 1 {_r(x + w)} {_r(y + radius)}'
            f'V{_r(y + h - radius)}'
            f'A{_r(radius)} {_r(radius)} 0 0 1 {_r(x + w - radius)} {_r(y + h)}'
            f'H{_r(x)}Z"/>')


def _mark(path, fill, tip):
    '''Wrap a mark in the group that carries its colour and its hover text.'''
    return (f'<g fill="{fill}"><title>{_esc(tip)}</title>{path}</g>')


def _svg_open(width, height, label):
    return (f'<svg class="chart" viewBox="0 0 {width} {height}" width="100%" '
            f'preserveAspectRatio="xMidYMin meet" role="img" '
            f'aria-label="{_esc(label)}">')


def _y_axis(ticks, top, plot_top, plot_h, plot_l, plot_r):
    '''Recessive hairline gridlines and their labels, one per tick.'''
    out = []
    for tick in ticks:
        y = plot_top + plot_h - (tick / top) * plot_h if top else plot_top + plot_h
        out.append(f'<line class="grid" x1="{plot_l}" y1="{_r(y)}" '
                   f'x2="{plot_r}" y2="{_r(y)}"/>')
        out.append(f'<text class="tick" x="{plot_l - 8}" y="{_r(y + 3.5)}" '
                   f'text-anchor="end">{_num(tick)}</text>')
    return out


# ---- charts ------------------------------------------------------------------


def error_histogram(items, x_title, y_title, fill="var(--series-1)",
                    unit="bp", max_bins=41):
    '''
    A histogram of one signed error, as columns centred on zero.

    ``items`` is a list of (value, locus label) pairs. The labels are carried through to
    the table so a bar can be read back to the loci that made it -- "which two loci are the
    ones 30bp out?" is the question a histogram provokes and cannot answer on its own.

    Returns (svg, rows) -- the figure, and the bin / count / member rows behind it so the
    page can offer the same numbers as a table. A chart whose values are only reachable by
    pointing at it is not readable in print or by a screen reader, so every figure here
    ships its table.
    '''
    if not items:
        return None, []
    values = [value for value, _label in items]
    centres, width, counts = signed_bins(values, max_bins)
    members = [[] for _ in centres]
    for value, label in items:
        members[int(round(value / width)) + (len(centres) - 1) // 2].append(label)
    top_count = max(counts)
    ticks, top = count_ticks(top_count)

    plot_l, plot_r = _PAD_L, _CHART_W - _PAD_R
    plot_w, plot_top = plot_r - plot_l, _PAD_T
    height = _PAD_T + _PLOT_H + _AXIS_H

    band = plot_w / len(centres)
    bar_w = min(_BAR_MAX, max(1.0, band - _GAP))

    out = [_svg_open(_CHART_W, height, f"{x_title}: {len(values)} loci")]
    out += _y_axis(ticks, top, plot_top, _PLOT_H, plot_l, plot_r)

    zero_index = centres.index(0) if 0 in centres else len(centres) // 2
    stride = _label_stride([_signed(c) for c in centres], band)
    for index, (centre, count) in enumerate(zip(centres, counts)):
        x = plot_l + index * band + (band - bar_w) / 2.0
        if count:
            h = (count / top) * _PLOT_H
            label = (f"{_signed(centre)} {unit}" if width == 1 else
                     f"{_signed(centre - width / 2.0)} to {_signed(centre + width / 2.0)} {unit}")
            out.append(_mark(_column(x, plot_top + _PLOT_H - h, bar_w, h), fill,
                             f"{label}: {count} loci"))
        if (index - zero_index) % stride == 0:
            out.append(f'<text class="tick" x="{_r(plot_l + index * band + band / 2.0)}" '
                       f'y="{plot_top + _PLOT_H + 15}" text-anchor="middle">'
                       f'{_esc(_signed(centre))}</text>')
    out.append(f'<line class="axis" x1="{plot_l}" y1="{plot_top + _PLOT_H}" '
               f'x2="{plot_r}" y2="{plot_top + _PLOT_H}"/>')
    out.append(f'<text class="axis-title" x="{_r((plot_l + plot_r) / 2.0)}" '
               f'y="{plot_top + _PLOT_H + 30}" text-anchor="middle">{_esc(x_title)}</text>')
    out.append(f'<text class="axis-title" transform="translate(13 {_r(plot_top + _PLOT_H / 2.0)}) '
               f'rotate(-90)" text-anchor="middle">{_esc(y_title)}</text>')
    out.append("</svg>")

    # The members are handed back as they came in, not joined: how a locus is written down
    # is the page's business, and the page turns each one into a link to its record.
    rows = [[_signed(c) if width == 1 else
             f"{_signed(c - width / 2.0)} .. {_signed(c + width / 2.0)}", n, m]
            for c, n, m in zip(centres, counts, members) if n]
    return "\n".join(out), rows


def outcome_columns(centres, width, series, x_title, y_title):
    '''
    One column per simulated event size, split into what became of the loci at that size.

    Stacked rather than overlaid: the column's full height is every locus simulated at that
    size -- the distribution the Tangram figure draws as its reference trace -- and the
    segments are a partition of it, so "how many were missed *here*" is read directly
    instead of by subtracting one trace from another.

    ``series`` is a list of (label, colour, {centre: [locus label, ...]}) -- the members
    rather than a tally, so the table under the figure can name the loci in each segment.
    '''
    totals = {centre: sum(len(members.get(centre, ())) for _l, _c, members in series)
              for centre in centres}
    top_count = max(totals.values(), default=0)
    if not top_count:
        return None, []
    ticks, top = count_ticks(top_count)

    plot_l, plot_r = _PAD_L, _CHART_W - _PAD_R
    plot_w, plot_top = plot_r - plot_l, _PAD_T
    height = _PAD_T + _PLOT_H + _AXIS_H

    band = plot_w / len(centres)
    bar_w = min(_BAR_MAX, max(1.0, band - _GAP))

    out = [_svg_open(_CHART_W, height, f"{x_title}, by outcome")]
    out += _y_axis(ticks, top, plot_top, _PLOT_H, plot_l, plot_r)

    labels = [_num(c) for c in centres]
    stride = _label_stride(labels, band)
    for index, centre in enumerate(centres):
        x = plot_l + index * band + (band - bar_w) / 2.0
        cursor = plot_top + _PLOT_H
        # Only the segment on top of the stack carries the column's data-end, so it is the
        # only one that is rounded; an interior segment rounded the same way would read as
        # the top of a column that is not there.
        stacked = [entry for entry in series if entry[2].get(centre)]
        for position, (label, colour, members) in enumerate(stacked):
            count = len(members[centre])
            h = (count / top) * _PLOT_H
            # The surface gap goes inside the segment, so the stack still sums to the
            # column height while touching segments stay separated by background rather
            # than by a stroke drawn around them.
            drawn = max(1.0, h - _GAP)
            radius = 4 if position == len(stacked) - 1 else 0
            out.append(_mark(_column(x, cursor - drawn, bar_w, drawn, radius), colour,
                             f"{_num(centre)} bp, {label}: {count} loci"))
            cursor -= h
        if index % stride == 0:
            out.append(f'<text class="tick" x="{_r(plot_l + index * band + band / 2.0)}" '
                       f'y="{plot_top + _PLOT_H + 15}" text-anchor="middle">'
                       f'{_esc(labels[index])}</text>')
    out.append(f'<line class="axis" x1="{plot_l}" y1="{plot_top + _PLOT_H}" '
               f'x2="{plot_r}" y2="{plot_top + _PLOT_H}"/>')
    out.append(f'<text class="axis-title" x="{_r((plot_l + plot_r) / 2.0)}" '
               f'y="{plot_top + _PLOT_H + 30}" text-anchor="middle">{_esc(x_title)}</text>')
    out.append(f'<text class="axis-title" transform="translate(13 {_r(plot_top + _PLOT_H / 2.0)}) '
               f'rotate(-90)" text-anchor="middle">{_esc(y_title)}</text>')
    out.append("</svg>")

    rows = []
    for centre in centres:
        if not totals[centre]:
            continue
        row = [_num(centre) if not width else
               f"{_num(centre - width // 2)} .. {_num(centre + width // 2)}", totals[centre]]
        row += [len(members.get(centre, ())) for _l, _c, members in series]
        # The members are kept split by segment, not flattened: a bin holding both
        # recovered and missed loci is exactly the case this figure exists to show, and a
        # bare list of everything at that size cannot say which of them were missed.
        row.append([(label, members[centre]) for label, _c, members in series
                    if members.get(centre)])
        rows.append(row)
    return "\n".join(out), rows


def stratum_bars(strata, series):
    '''
    One group of horizontal bars per stratum, one bar per metric.

    Horizontal because the labels are words -- family names, "0.25-0.5", "nested INS" --
    and words set along a vertical axis stay readable at any length, where a column chart
    would either rotate them or truncate them. Every bar is a rate on the same 0-100%
    axis, which is the only reason three metrics can share one plot at all.

    ``series`` is a list of (label, colour, getter). Metrics that are unavailable for the
    whole run -- genotype concordance with no genomes paired -- are dropped, rather than
    drawn as a row of empty tracks.
    '''
    live = [s for s in series
            if any(s[2](stats) is not None for stats in strata.values())]
    # A stratification into one group has nothing to compare against: the bars would
    # restate the headline metrics the tiles already carry. The table still shows it.
    if len(strata) < 2 or not live:
        return None

    bar_h, bar_gap, group_gap = 11, _GAP, 16
    group_h = len(live) * (bar_h + bar_gap) - bar_gap + group_gap
    label_w = 168
    plot_l = label_w + 10
    plot_r = _CHART_W - _PAD_R - 34          # room for the value at the tip of a full bar
    plot_w = plot_r - plot_l
    plot_top = _PAD_T + 4
    height = plot_top + len(strata) * group_h + 20

    out = [_svg_open(_CHART_W, height, "Recall and accuracy by stratum")]
    for fraction in (0, 0.25, 0.5, 0.75, 1.0):
        x = plot_l + fraction * plot_w
        out.append(f'<line class="grid" x1="{_r(x)}" y1="{plot_top - 4}" '
                   f'x2="{_r(x)}" y2="{_r(height - 18)}"/>')
        out.append(f'<text class="tick" x="{_r(x)}" y="{_r(height - 5)}" '
                   f'text-anchor="middle">{int(fraction * 100)}%</text>')

    for row, (label, stats) in enumerate(strata.items()):
        top = plot_top + row * group_h
        centre = top + (len(live) * (bar_h + bar_gap) - bar_gap) / 2.0
        out.append(f'<text class="row-label" x="{label_w}" y="{_r(centre - 1)}" '
                   f'text-anchor="end">{_esc(label)}</text>')
        out.append(f'<text class="row-note" x="{label_w}" y="{_r(centre + 11)}" '
                   f'text-anchor="end">n={stats["n_loci"]}</text>')
        for index, (name, colour, getter) in enumerate(live):
            value = getter(stats)
            y = top + index * (bar_h + bar_gap)
            if value is None:
                continue
            w = max(1.0, value * plot_w)
            out.append(_mark(_bar(plot_l, y, w, bar_h), colour,
                             f"{label} — {name}: {_pct(value)}"))
    out.append("</svg>")
    return "\n".join(out)


def _legend(entries):
    '''
    The identity channel for two or more series. A single series needs none: there is only
    one colour on the plot and the title already says what it is, so a one-swatch box just
    restates the heading.
    '''
    if len(entries) < 2:
        return ""
    swatches = "".join(
        f'<li><span class="swatch" style="background:{colour}"></span>{_esc(label)}</li>'
        for label, colour in entries)
    return f'<ul class="legend">{swatches}</ul>'


# ---- tables ------------------------------------------------------------------


def _table(headers, rows, sortable=False, numeric_from=1, wrap_cols=()):
    '''
    A table. ``numeric_from`` is the first column to right-align and sort numerically --
    every table here leads with a label and continues with numbers. ``wrap_cols`` names the
    columns that hold free text instead, which are left-aligned and allowed to wrap: a list
    of every locus in a bar is as long as the bar is tall, and holding it on one line would
    push the counts off the side of the page.
    '''
    if not rows:
        return '<p class="empty">(nothing to show)</p>'

    def cell(tag, index, value):
        css = ("ids" if index in wrap_cols else
               "num" if index >= numeric_from else "")
        body = value if isinstance(value, Html) else _esc(value)
        return f'<{tag}{f" class={css}" if css else ""}>{body}</{tag}>'

    head = "".join(cell("th", i, h) for i, h in enumerate(headers))
    body = []
    for row in rows:
        cells = "".join(cell("td", i, c) for i, c in enumerate(row))
        body.append(f"<tr>{cells}</tr>")
    css = "data sortable" if sortable else "data"
    return (f'<table class="{css}"><thead><tr>{head}</tr></thead>'
            f'<tbody>{"".join(body)}</tbody></table>')


def _details(summary_text, body):
    return (f'<details class="table-view"><summary>{_esc(summary_text)}</summary>'
            f'{body}</details>')


def _recall(stats):
    '''Locus recall; an Evaluate JSON from before it was named that has only recovery_rate.'''
    return stats.get("recall", stats.get("recovery_rate"))


def _f1(stats, cls):
    '''An allele class's F1, from the counts block, or the old carrier F1 where there is none.'''
    counts = stats.get("counts")
    if counts:
        return counts[cls]["f1"]
    return stats["carriers"]["f1"] if cls == "carrier" else None


COUNT_HEADERS = ["level", "total", "called", "TP", "FP", "FN", "recall", "precision", "F1"]


def counts_rows(counts, unit):
    rows = []
    for key, label in (("loci", "loci"), ("carrier", f"carrier {unit}s"),
                       ("noncarrier", f"noncarrier {unit}s"), ("all", f"all {unit}s")):
        c = counts[key]
        rows.append([label, c["total"], c["genotyped"], c["tp"], _num(c["fp"]), c["fn"],
                     _num(c["recall"]), _num(c["precision"]), _num(c["f1"])])
    return rows


def _counts_section(summary, meta):
    counts = summary["overall"].get("counts")
    if not counts or not meta.get("sample_pairs"):
        return ""
    unit = meta.get("genotype_label", "genotype")
    return (
        '<section class="card"><h2>By locus and by ' + _esc(unit) + '</h2>'
        f'<p class="sub">Over {len(meta["sample_pairs"])} paired genomes. Each allele class is '
        'the positive class in turn; "all" sums the two.</p>'
        + _table(COUNT_HEADERS, counts_rows(counts, unit))
        + '<p class="note">A locus is a TP when a call was placed on it and an FN when none was; a '
        'call on no simulated locus is a locus FP. A ' + _esc(unit) + ' is a TP when genotyped as '
        'its own class, an FN when genotyped as the other class or never genotyped (missed or '
        'displaced locus, or a "." genotype), and an FP of the class it was wrongly genotyped as; '
        'one never genotyped is an FP of neither. The carriers of a call on no simulated locus are '
        'carrier FPs. "called" is the loci recalled and the ' + _esc(unit) + 's that got a '
        'genotype.</p></section>')


def stratum_rows(strata, unit):
    '''The breakdown table, with the same columns the text report prints.'''
    headers = ["stratum", unit, "found", "displaced", "recall", "allele ok",
               "mean offset", "mean length err", "carrier F1", "noncarrier F1", "all F1",
               "genotype"]
    rows = []
    for label, stats in strata.items():
        rows.append([
            label,
            stats["n_loci"],
            stats["n_recovered"],
            stats["n_displaced"],
            _pct(_recall(stats)),
            _pct(stats["allele_concordance_rate"]),
            _signed(stats["breakpoint_offset_bp"]["mean"]),
            _signed(stats["allele_length_error_bp"]["mean"]),
            _num(_f1(stats, "carrier")),
            _num(_f1(stats, "noncarrier")),
            _num(_f1(stats, "all")),
            _pct(stats["genotypes"]["concordance"]),
        ])
    return headers, rows


def resolution_rows(offsets, max_dist):
    '''
    The cumulative breakpoint-resolution readout: of the calls that were placed, how many
    landed exactly on the simulated breakpoint, and how many within N bp of it. This is the
    figure that is quoted in prose ("65% exact, over 99% within 15bp"), and it is the one
    number a mean and SD cannot be turned back into.
    '''
    if not offsets:
        return []
    steps = [s for s in RESOLUTION_STEPS if s < max_dist] + [max_dist]
    rows = []
    for step in steps:
        within = sum(1 for offset in offsets if abs(offset) <= step)
        label = "exact (0 bp)" if step == 0 else f"within {step} bp"
        rows.append([label, within, _pct(within / len(offsets))])
    return rows


# ---- page --------------------------------------------------------------------


CSS = '''
:root {
  color-scheme: light;
  --surface:        #fcfcfb;
  --plane:          #f9f9f7;
  --ink:            #0b0b0b;
  --ink-2:          #52514e;
  --ink-muted:      #898781;
  --grid:           #e1e0d9;
  --axis:           #c3c2b7;
  --border:         rgba(11, 11, 11, 0.10);
  --series-1:       #2a78d6;
  --series-2:       #eb6834;
  --series-3:       #1baf7a;
  --step-near:      #86b6ef;
  --step-mid:       #3987e5;
  --step-far:       #184f95;
}
@media (prefers-color-scheme: dark) {
  :root:not([data-theme="light"]) {
    color-scheme: dark;
    --surface:      #1a1a19;
    --plane:        #0d0d0d;
    --ink:          #ffffff;
    --ink-2:        #c3c2b7;
    --ink-muted:    #898781;
    --grid:         #2c2c2a;
    --axis:         #383835;
    --border:       rgba(255, 255, 255, 0.10);
    --series-1:     #3987e5;
    --series-2:     #d95926;
    --series-3:     #199e70;
    /* The same ramp, stepped for the dark surface: the outcome that matters most
       stays the step furthest from the background in either mode. */
    --step-near:    #184f95;
    --step-mid:     #3987e5;
    --step-far:     #86b6ef;
  }
}
* { box-sizing: border-box; }
body {
  margin: 0;
  padding: 32px 24px 64px;
  background: var(--plane);
  color: var(--ink);
  font: 15px/1.55 system-ui, -apple-system, "Segoe UI", sans-serif;
}
main { max-width: 980px; margin: 0 auto; }
h1 { font-size: 22px; font-weight: 600; margin: 0 0 4px; letter-spacing: -0.01em; }
h2 { font-size: 17px; font-weight: 600; margin: 0 0 4px; }
h3 { font-size: 14px; font-weight: 600; margin: 22px 0 8px; color: var(--ink-2); }
p  { margin: 0 0 10px; }
a  { color: var(--series-1); }
.sub { color: var(--ink-2); font-size: 13.5px; }
.note { color: var(--ink-muted); font-size: 12.5px; margin: 8px 0 0; }
code { font-family: ui-monospace, SFMono-Regular, Menlo, monospace; font-size: 12.5px; }

header { margin-bottom: 22px; }
.meta { display: flex; flex-wrap: wrap; gap: 6px; margin-top: 10px; }
.meta span {
  border: 1px solid var(--border); border-radius: 999px;
  padding: 2px 10px; font-size: 12px; color: var(--ink-2); background: var(--surface);
}

.card {
  background: var(--surface); border: 1px solid var(--border);
  border-radius: 10px; padding: 20px 22px; margin-bottom: 18px;
}
.hero { display: flex; align-items: baseline; gap: 16px; flex-wrap: wrap; }
.hero .figure { font-size: 52px; font-weight: 600; line-height: 1; letter-spacing: -0.02em; }
.hero .caption { color: var(--ink-2); font-size: 14px; }

.tiles { display: grid; grid-template-columns: repeat(auto-fit, minmax(168px, 1fr)); gap: 1px;
         background: var(--border); border: 1px solid var(--border);
         border-radius: 10px; overflow: hidden; margin-bottom: 18px; }
.tile { background: var(--surface); padding: 14px 16px; }
.tile .label { font-size: 12px; color: var(--ink-2); }
.tile .value { font-size: 22px; font-weight: 600; margin-top: 2px; }
.tile .aside { font-size: 12px; color: var(--ink-muted); margin-top: 1px; }

.chart { display: block; margin: 6px 0 2px; overflow: visible; }
.chart .grid { stroke: var(--grid); stroke-width: 1; }
.chart .axis { stroke: var(--axis); stroke-width: 1; }
.chart .tick { fill: var(--ink-muted); font-size: 11px; font-variant-numeric: tabular-nums; }
.chart .axis-title { fill: var(--ink-2); font-size: 12px; }
.chart .row-label { fill: var(--ink); font-size: 12px; }
.chart .row-note { fill: var(--ink-muted); font-size: 11px; font-variant-numeric: tabular-nums; }
.chart g:hover path { opacity: 0.82; }

.legend { display: flex; flex-wrap: wrap; gap: 16px; list-style: none;
          margin: 10px 0 0; padding: 0; font-size: 12.5px; color: var(--ink-2); }
.legend li { display: flex; align-items: center; gap: 7px; }
.swatch { width: 11px; height: 11px; border-radius: 3px; display: inline-block; }

.split { display: grid; grid-template-columns: minmax(0, 1fr) 224px; gap: 24px; align-items: start; }
@media (max-width: 760px) { .split { grid-template-columns: minmax(0, 1fr); } }

table.data { border-collapse: collapse; width: 100%; font-size: 12.5px; }
table.data th, table.data td {
  text-align: left; padding: 5px 10px 5px 0; border-bottom: 1px solid var(--border);
  white-space: nowrap;
}
table.data th { color: var(--ink-2); font-weight: 600; }
table.data td.num, table.data th.num { text-align: right; font-variant-numeric: tabular-nums; }
/* The member list of a bar is as long as the bar is tall, so it wraps inside its cell
   rather than pushing the counts off the side of the page. */
table.data td.ids, table.data th.ids {
  white-space: normal; word-break: break-word; min-width: 220px;
  font-family: ui-monospace, SFMono-Regular, Menlo, monospace;
  font-size: 11.5px; color: var(--ink-2);
}
table.data.sortable th { cursor: pointer; user-select: none; }
table.data.sortable th:hover { color: var(--ink); }
table.data.sortable th[aria-sort]::after { content: " \\2191"; font-size: 10px; }
table.data.sortable th[aria-sort="descending"]::after { content: " \\2193"; }
.scroll { overflow-x: auto; }

details.table-view { margin-top: 12px; }
details.table-view > summary {
  cursor: pointer; font-size: 12.5px; color: var(--ink-2); padding: 4px 0;
}
details.table-view > summary:hover { color: var(--ink); }
details.table-view[open] > summary { margin-bottom: 6px; }
.empty { color: var(--ink-muted); font-size: 12.5px; }

footer { color: var(--ink-muted); font-size: 12px; margin-top: 28px; }

.tabs { display: flex; gap: 4px; margin: 0 0 18px; border-bottom: 1px solid var(--border); }
.tabs[hidden] { display: none; }
.tabs button {
  appearance: none; background: none; border: 0; border-bottom: 2px solid transparent;
  margin-bottom: -1px; padding: 8px 14px; cursor: pointer;
  font: inherit; font-size: 13.5px; color: var(--ink-2);
}
.tabs button:hover { color: var(--ink); }
.tabs button[aria-selected="true"] {
  color: var(--ink); font-weight: 600; border-bottom-color: var(--series-1);
}
.file-link, .locus-link { color: var(--series-1); text-decoration: none; }
.file-link:hover, .locus-link:hover { text-decoration: underline; }
.locus-link { white-space: nowrap; }

/* ---- the side-by-side VCF viewer ---- */
.viewer { display: grid; grid-template-columns: 1fr 1fr; gap: 14px; margin-top: 4px; }
@media (max-width: 760px) { .viewer { grid-template-columns: minmax(0, 1fr); } }
/* min-width:0 or a grid column refuses to shrink below its widest record and the two
   panes push the page sideways instead of scrolling inside themselves. */
.pane-wrap { position: relative; min-width: 0; }
.pane-head { font-size: 12px; color: var(--ink-2); margin-bottom: 5px; }
.pane {
  height: 440px; overflow: auto; overscroll-behavior: contain;
  border: 1px solid var(--border); border-radius: 8px; background: var(--surface);
  /* Positioned so a row's offsetTop is measured from the pane's own content box, which is
     the coordinate its scrollTop is in. Left static, offsetParent would be .pane-wrap and
     every jump would land one pane-heading too low. */
  position: relative;
}
.vcf-row {
  display: flex; gap: 8px; align-items: baseline; width: max-content; min-width: 100%;
  padding: 2px 8px; cursor: pointer; border-left: 3px solid transparent;
  font-family: ui-monospace, SFMono-Regular, Menlo, monospace; font-size: 11px;
  white-space: pre; color: var(--ink-2);
}
.vcf-row:hover { background: var(--plane); }
.vcf-row:focus-visible { outline: 2px solid var(--series-1); outline-offset: -2px; }
.vcf-row .tag {
  flex: none; width: 38px; font-size: 10px; letter-spacing: 0.03em;
  color: var(--ink-muted); text-transform: uppercase;
}
.vcf-row.match { border-left-color: var(--series-1); }
.vcf-row.displ { border-left-color: var(--series-2); }
.vcf-row.alone { border-left-color: var(--series-3); }
.vcf-row.skip  { opacity: 0.55; }
/* The picked pair is marked on both sides at once -- that pairing is the thing the
   viewer exists to show, and a highlight on only the clicked side does not show it. */
.vcf-row.picked { background: var(--plane); color: var(--ink); font-weight: 600; }
.bubble {
  position: absolute; left: 50%; top: 50%; transform: translate(-50%, -50%);
  max-width: 84%; padding: 9px 13px; border-radius: 8px; z-index: 2;
  background: var(--ink); color: var(--surface); font-size: 12px; line-height: 1.4;
  box-shadow: 0 3px 14px rgba(0, 0, 0, 0.28);
}
.bubble[hidden] { display: none; }
'''

# Column sorting, so a family table of any length can be read by the column that matters.
# Written out rather than pulled from a library: the page has to work with no network.
JS = '''
document.querySelectorAll("table.sortable").forEach(function (table) {
  var body = table.tBodies[0];
  var rows = Array.prototype.slice.call(body.rows);
  Array.prototype.forEach.call(table.tHead.rows[0].cells, function (th, column) {
    th.addEventListener("click", function () {
      var descending = th.getAttribute("aria-sort") !== "descending";
      Array.prototype.forEach.call(table.tHead.rows[0].cells, function (other) {
        other.removeAttribute("aria-sort");
      });
      th.setAttribute("aria-sort", descending ? "descending" : "ascending");
      var read = function (row) {
        var text = row.cells[column].textContent.trim().replace(/[,%]/g, "");
        var number = parseFloat(text);
        return isNaN(number) ? text.toLowerCase() : number;
      };
      rows.slice().sort(function (a, b) {
        var x = read(a), y = read(b);
        if (x === y) { return 0; }
        return (x < y ? -1 : 1) * (descending ? -1 : 1);
      }).forEach(function (row) { body.appendChild(row); });
    });
  });
});

// The side-by-side viewer, the tab that holds it, and the links into it. Two panes are
// kept together two ways: a click follows the pairing Evaluate worked out, and a scroll
// falls back to relative position, since two files of different lengths have no
// record-for-record correspondence to follow.
(function () {
  var viewer = document.getElementById("vcf-viewer");
  if (!viewer) { return; }
  var panes = { truth: document.getElementById("pane-truth"),
                pred: document.getElementById("pane-pred") };
  var bubbles = { truth: document.getElementById("bubble-truth"),
                  pred: document.getElementById("bubble-pred") };
  if (!panes.truth || !panes.pred) { return; }

  // ---- tabs ----
  var tabs = document.getElementById("tabs");
  var panels = { report: document.getElementById("panel-report"),
                 vcf: document.getElementById("panel-vcf") };
  var buttons = { report: document.getElementById("tab-report"),
                  vcf: document.getElementById("tab-vcf") };

  function showTab(name) {
    if (!tabs) { return; }
    Object.keys(panels).forEach(function (key) {
      if (panels[key]) { panels[key].hidden = key !== name; }
      if (buttons[key]) { buttons[key].setAttribute("aria-selected", key === name); }
    });
  }

  if (tabs && panels.report && panels.vcf) {
    tabs.hidden = false;
    showTab("report");
    Object.keys(buttons).forEach(function (name) {
      if (buttons[name]) {
        buttons[name].addEventListener("click", function () { showTab(name); });
      }
    });
  }

  // Scrolling one pane scrolls the other, whose own scroll event would scroll the first
  // back -- so whichever pane the user is driving owns the sync until it goes quiet, and
  // the events its movement induces elsewhere are ignored.
  var owner = null, release = null;
  // A jump moves both panes deliberately, and neither move may be answered by the sync:
  // the mate has been put where it belongs and proportional scrolling would drag it off.
  var muted = false, unmute = null;
  function claim(pane) {
    owner = pane;
    clearTimeout(release);
    release = setTimeout(function () { owner = null; }, 150);
  }
  function quiet(move) {
    muted = true;
    try { move(); } finally {
      clearTimeout(unmute);
      // Scroll events arrive after the assignment, so the mute outlives this turn.
      unmute = setTimeout(function () { muted = false; }, 150);
    }
  }

  function proportional(from, to) {
    var fromMax = from.scrollHeight - from.clientHeight;
    var toMax = to.scrollHeight - to.clientHeight;
    // A pane shorter than its box has nowhere to scroll to; 0/0 would be NaN.
    if (fromMax <= 0 || toMax <= 0) { return; }
    to.scrollTop = (from.scrollTop / fromMax) * toMax;
  }

  function hideBubbles() {
    Object.keys(bubbles).forEach(function (side) {
      if (bubbles[side]) { bubbles[side].hidden = true; }
    });
  }

  function clearPicked() {
    Array.prototype.forEach.call(viewer.querySelectorAll(".vcf-row.picked"),
      function (row) { row.classList.remove("picked"); });
  }

  function reveal(pane, row) {
    pane.scrollTop = row.offsetTop - (pane.clientHeight - row.offsetHeight) / 2;
  }

  // `arriving` is set when the pick came from a link rather than from a click on the row
  // itself, in which case the row is not on screen yet and its own pane has to move too.
  function pick(side, row, arriving) {
    var pane = panes[side];
    var otherSide = side === "truth" ? "pred" : "truth";
    var other = panes[otherSide];
    hideBubbles();
    clearPicked();
    row.classList.add("picked");
    var mate = row.getAttribute("data-match");
    mate = mate ? document.getElementById(mate) : null;
    if (mate) { mate.classList.add("picked"); }
    quiet(function () {
      if (arriving) { reveal(pane, row); }
      if (mate) {
        reveal(other, mate);
      } else {
        // Nothing to jump to. Go where the record would be if it were there, and say why
        // it is not -- an unmoved pane cannot be told from a viewer that did nothing.
        proportional(pane, other);
      }
    });
    if (mate) { return; }
    var bubble = bubbles[otherSide];
    if (bubble) {
      bubble.textContent = row.getAttribute("data-note") || "No matching record.";
      bubble.hidden = false;
    }
  }

  Object.keys(panes).forEach(function (side) {
    var pane = panes[side];
    var other = panes[side === "truth" ? "pred" : "truth"];
    pane.addEventListener("scroll", function () {
      if (muted || (owner && owner !== pane)) { return; }
      claim(pane);
      // The bubble points at a place in the pane, so it stops being true as soon as the
      // pane moves. The picked pair is kept: scrolling away to read around a match and
      // back is the normal thing to do, and losing the highlight for it would be hostile.
      hideBubbles();
      proportional(pane, other);
    });
    pane.addEventListener("click", function (event) {
      var row = event.target;
      while (row && row !== pane && !row.classList.contains("vcf-row")) {
        row = row.parentNode;
      }
      if (row && row.classList && row.classList.contains("vcf-row")) { pick(side, row); }
    });
    // Same from the keyboard: the rows are focusable, so they have to be operable.
    pane.addEventListener("keydown", function (event) {
      if (event.key !== "Enter" && event.key !== " ") { return; }
      var row = event.target;
      if (row && row.classList && row.classList.contains("vcf-row")) {
        event.preventDefault();
        pick(side, row);
      }
    });
  });

  // A file name opens its pane; a locus name in any table opens the record behind it.
  // Both are real anchors, so without this handler they still jump to the right element
  // on a page whose tabs were never switched on.
  document.addEventListener("click", function (event) {
    var link = event.target;
    while (link && link !== document && link.tagName !== "A") { link = link.parentNode; }
    if (!link || link.tagName !== "A") { return; }
    var locus = link.getAttribute("data-locus");
    var side = link.getAttribute("data-pane");
    if (!locus && !side) { return; }
    event.preventDefault();
    showTab("vcf");
    if (locus) {
      var row = document.getElementById(locus);
      if (row) { pick("truth", row, true); }
    }
  });
})();
'''


def _tile(label, value, aside=None):
    extra = f'<div class="aside">{_esc(aside)}</div>' if aside else ""
    return (f'<div class="tile"><div class="label">{_esc(label)}</div>'
            f'<div class="value">{_esc(value)}</div>{extra}</div>')


def _file_name(path, side, linked):
    '''
    A file name, as the way into the pane holding it when there is one to go to.

    Only linked when the viewer was built: a name pointing at a pane that is not on the
    page is a link that does nothing, which is worse than plain text.
    '''
    name = _esc(path)
    if not linked:
        return f"<code>{name}</code>"
    return (f'<code><a class="file-link" href="#pane-{side}" '
            f'data-pane="{side}">{name}</a></code>')


def _header(meta, has_viewer=False):
    pairs = meta.get("sample_pairs") or []
    shown = ", ".join(f"{t}={p}" if t != p else t for t, p in pairs[:6])
    if len(pairs) > 6:
        shown += f", … (+{len(pairs) - 6})"
    lines = [
        '<header>',
        '<h1>tevarsim Evaluate</h1>',
        f'<p class="sub">{meta["n_loci"]} simulated loci over '
        f'{meta["n_truth_samples"]} genomes, against {meta["n_pred_records"]} '
        f'{meta["predType"]} records over {meta["n_pred_samples"]} samples</p>',
        # The file names are the way into the viewer: the obvious thing to click when you
        # want to see a file is its name. They are plain <a href="#pane-...">, so with no
        # JavaScript running they still jump to the pane further down the page.
        f'<p class="sub">truth {_file_name(meta["truth"], "truth", has_viewer)}'
        f' &nbsp;&middot;&nbsp; '
        f'prediction {_file_name(meta["pred"], "pred", has_viewer)}</p>',
    ]
    if pairs:
        lines.append(f'<p class="note">{len(pairs)} genome(s) paired by '
                     f'{_esc(meta["pairing"])}: {_esc(shown)}</p>')
    else:
        lines.append('<p class="note">No genomes could be paired, so carrier and genotype '
                     'statistics are unavailable. Pass --sample_map to pair them.</p>')
    chips = [f"--max_dist {meta['max_dist']}", f"--gt_len_tol {meta['gt_len_tol']}"]
    if meta.get("nHap", 1) > 1:
        chips.append(f"--nHap {meta['nHap']}")
    if meta.get("filters"):
        chips.append(meta["filters"])
    lines.append('<div class="meta">'
                 + "".join(f"<span>{_esc(chip)}</span>" for chip in chips)
                 + "</div>")
    if meta.get("filters"):
        lines.append('<p class="note">The truth is filtered but the prediction is not, so '
                     'unmatched predictions include events the filter removed.</p>')
    lines.append("</header>")
    return "\n".join(lines)


def _headline(summary, meta):
    overall = summary["overall"]
    predictions = summary["predictions"]
    carriers = overall["carriers"]
    genotypes = overall["genotypes"]
    label = meta.get("genotype_label", "genotype")

    hero = (
        '<section class="card hero">'
        f'<div class="figure">{_esc(_pct(_recall(overall)))}</div>'
        f'<div class="caption">locus recall over the {overall["n_loci"]} simulated loci'
        f' &mdash; {overall["n_recovered"]} found, '
        f'{overall["n_loci"] - overall["n_recovered"]} missed</div>'
        '</section>'
    )

    tiles = [
        _tile("Correctly placed", _num(overall["n_detected"]),
              f'{overall["n_displaced"]} displaced past the solo LTR'
              if overall["n_displaced"] else "no displaced calls"),
        _tile("Allele concordant", _pct(overall["allele_concordance_rate"]),
              f'{overall["n_allele_concordant"]} of {overall["n_detected"]} placed'),
        _tile("Locus precision", _pct(predictions["precision"]),
              f'{predictions["unmatched"]} of {predictions["total"]} calls unmatched'),
        _tile("Breakpoint offset", _spread_text(overall["breakpoint_offset_bp"]),
              f'mean ± SD over n={overall["breakpoint_offset_bp"]["n"]}'),
        _tile("Allele length error", _spread_text(overall["allele_length_error_bp"]),
              f'mean ± SD over n={overall["allele_length_error_bp"]["n"]}'),
    ]
    if meta.get("sample_pairs"):
        c = overall.get("counts", {}).get("carrier") or carriers
        tiles.append(_tile("Carrier F1", _num(c["f1"]),
                           f'TP {c["tp"]}, FP {c["fp"]}, FN {c["fn"]}'))
        tiles.append(_tile(f"{label.capitalize()} concordance",
                           _pct(genotypes["concordance"]),
                           f'{genotypes["concordant"]} of {genotypes["compared"]} compared'))
    return hero + f'<section class="tiles">{"".join(tiles)}</section>'


def _breakpoint_section(loci, labels, meta):
    placed = [(locus["match"]["pos_offset"], label)
              for locus, label in zip(loci, labels) if locus["detected"]]
    offsets = [offset for offset, _label in placed]
    displaced = sum(1 for locus in loci if locus.get("displaced"))
    # Said in both branches: a reader looking at an empty or thinned-out figure is owed the
    # reason it is empty, and "every recovery was displaced" is a result, not a blank.
    note = ("A displaced call — one anchored past the surviving solo LTR rather than at "
            "its start — earns detection credit but no breakpoint credit, so the "
            f"{displaced} of them are not plotted here."
            if displaced else
            "Negative is short of the simulated breakpoint, positive is past it.")
    svg, rows = error_histogram(
        placed,
        "breakpoint error (called position − simulated position, bp)",
        "loci")
    if svg is None:
        return ('<section class="card"><h2>Breakpoint resolution</h2>'
                '<p class="empty">No call was placed on a simulated locus.</p>'
                f'<p class="note">{_esc(note)}</p></section>')
    resolution = resolution_rows(offsets, meta["max_dist"])
    return (
        '<section class="card">'
        '<h2>Breakpoint resolution</h2>'
        f'<p class="sub">Signed error of every correctly placed call, one bar per bin. '
        f'n={len(offsets)}.</p>'
        '<div class="split">'
        f'<div>{svg}<p class="note">{_esc(note)}</p></div>'
        f'<div><h3>Cumulative</h3>'
        f'{_table(["placed within", "loci", "share"], resolution)}</div>'
        '</div>'
        + _details("Show the values behind this figure",
                   _table(["breakpoint error (bp)", "count", "loci"],
                          [[a, b, _refs_html(m)] for a, b, m in rows], wrap_cols={2}))
        + '</section>'
    )


def _length_section(loci, labels):
    errors = [(locus["match"]["length_error"], label)
              for locus, label in zip(loci, labels)
              if locus["detected"] and locus["match"]["length_error"] is not None]
    svg, rows = error_histogram(
        errors,
        "allele length error (called length − simulated length, bp)", "loci")
    if svg is None:
        return ""
    return (
        '<section class="card">'
        '<h2>Allele length accuracy</h2>'
        f'<p class="sub">How far the called allele is from the simulated element, for the '
        f'{len(errors)} placed calls that carry a comparable sequence.</p>'
        f'{svg}'
        '<p class="note">Negative is short of the simulated element, positive is past it. '
        'A well-placed call that truncated the element shows up here and not in the '
        'breakpoint figure.</p>'
        + _details("Show the values behind this figure",
                   _table(["allele length error (bp)", "count", "loci"],
                          [[a, b, _refs_html(m)] for a, b, m in rows], wrap_cols={2}))
        + '</section>'
    )


def _segment_html(segments):
    '''
    A bin's loci, named by the segment they fall in. A single-outcome bin has nothing to
    tell apart, so it is written as a plain run of names.
    '''
    if len(segments) == 1:
        return _refs_html(segments[0][1])
    return Html("; ".join(f"{_esc(label)}: {_refs_html(refs)}"
                          for label, refs in segments))


def _size_section(loci, labels):
    '''
    The size distribution of every simulated locus, stacked by outcome -- the shape that
    exposes a sensitivity floor at a particular element size, which no aggregate can.
    '''
    sized = [(locus, label) for locus, label in zip(loci, labels)
             if locus.get("size_bp") is not None]
    if not sized:
        return ""
    centres, width, key = value_bins([locus["size_bp"] for locus, _label in sized])
    placed, displaced, missed = {}, {}, {}
    for locus, label in sized:
        centre = key(locus["size_bp"])
        bucket = (placed if locus["detected"] else
                  displaced if locus.get("displaced") else missed)
        bucket.setdefault(centre, []).append(label)
    series = [
        ("correctly placed", "var(--step-near)", placed),
        ("displaced", "var(--step-mid)", displaced),
        ("missed", "var(--step-far)", missed),
    ]
    live = [entry for entry in series if entry[2]]
    svg, rows = outcome_columns(
        centres, width, live,
        "simulated event size (net length change, bp)", "loci")
    if svg is None:
        return ""
    unsized = len(loci) - len(sized)
    note = ("Sizes keep their sign, so an insertion reads positive and an excision or "
            "deletion negative; the bar heights are the simulated distribution and the "
            "segments partition it.")
    if width:
        note += (f" Sizes are binned to {width} bp, so a bar is labelled with the centre "
                 "of its bin; the table gives the range.")
    if unsized:
        note += (f" {unsized} locus/loci carry only symbolic alleles and have no length "
                 "to place on this axis.")
    return (
        '<section class="card">'
        '<h2>Recall by event size</h2>'
        f'<p class="sub">Every simulated locus, at the size it was simulated at, split by '
        f'what the prediction made of it.</p>'
        f'{svg}'
        + _legend([(label, colour) for label, colour, _counts in live])
        + f'<p class="note">{_esc(note)}</p>'
        + _details("Show the values behind this figure",
                   _table(["event size (bp)", "count"]
                          + [label for label, _c, _n in live] + ["loci"],
                          [row[:-1] + [_segment_html(row[-1])] for row in rows],
                          wrap_cols={len(live) + 2}))
        + '</section>'
    )


def _strata_section(summary, meta):
    metrics = [("locus recall", "var(--series-1)", _recall)]
    # With no genomes paired there are no carrier counts to score, and
    # ``calculate_metrics`` reports an F1 of 0 for the empty confusion matrix -- a bar
    # drawn from that would read as "the caller got every carrier wrong" rather than
    # "nothing was compared". The tables still carry the raw zeros, as the text report does.
    if meta.get("sample_pairs"):
        metrics.append(("carrier F1", "var(--series-2)", lambda s: _f1(s, "carrier")))
        metrics.append(("genotype concordance", "var(--series-3)",
                        lambda s: s["genotypes"]["concordance"]))
    out = []
    for title, key in STRATUM_TITLES:
        strata = summary.get(key)
        if not strata:
            continue
        unit = "events" if key == "by_event_type" else "loci"
        headers, rows = stratum_rows(strata, unit)
        svg = stratum_bars(strata, metrics)
        live = [(label, colour) for label, colour, getter in metrics
                if any(getter(stats) is not None for stats in strata.values())]
        block = [f'<section class="card"><h2>{_esc(title)}</h2>']
        if key == "by_event_type":
            block.append('<p class="sub">A locus counts under each event class in its '
                         'history, so these rows count events, not loci, and total more '
                         'than the loci.</p>')
        if svg:
            block.append(svg)
            block.append(_legend(live))
        block.append(_details("Show the full table",
                              f'<div class="scroll">'
                              f'{_table(headers, rows, sortable=True)}</div>'))
        block.append("</section>")
        out.append("".join(block))
    return "\n".join(out)


# A file with more records than this is shown truncated rather than inlined whole: the page
# has to stay a file a browser will open. The cut is always stated on the page -- a viewer
# that silently stops part way through reads as a shorter VCF than the one that was scored.
MAX_VIEWER_ROWS = 5000

# REF and ALT hold whole elements. A 6 kb ALT per record is 300x the rest of the line and
# would put megabytes of sequence in the page, so an allele past this length is shown as its
# first bases and its length -- the same thing any VCF viewer does, and the whole sequence
# is in the VCF the pane names.
ALLELE_PREVIEW_BP = 24


def _elide_allele(allele):
    if len(allele) <= ALLELE_PREVIEW_BP + 12:
        return allele
    return f"{allele[:ALLELE_PREVIEW_BP]}…({len(allele):,} bp)"


def _vcf_text(text):
    '''One record as it stands in the file, with its allele sequences shortened.'''
    fields = text.split("\t")
    if len(fields) >= 8:                      # a VCF line; a BED one is left alone
        for column in (3, 4):                 # REF, ALT
            fields[column] = ",".join(_elide_allele(a) for a in fields[column].split(","))
    return "  ".join(fields)


def _viewer_rows(rows, side, mate_side, mate_of, missing_note):
    '''
    Render one pane's records.

    ``mate_of(scored)`` returns the record on the other side this one is paired with, or
    None. Every row carries either the DOM id of its mate or the sentence to show when
    there is none, so the page needs no lookup table beside the markup.
    '''
    out, states = [], {"match": 0, "displ": 0, "alone": 0, "skip": 0}
    for row in rows[:MAX_VIEWER_ROWS]:
        where = f"{row.chrom}:{row.pos}"
        if row.scored is None:
            state, mate, note = "skip", None, f"{where} was not scored — {row.skipped}."
        else:
            mate, displaced = mate_of(row.scored)
            state = "displ" if displaced else "match" if mate is not None else "alone"
            note = None if mate is not None else missing_note(where)
        states[state] += 1
        anchor = (f' data-match="{mate_side}{mate.index}"' if mate is not None
                  else f' data-note="{_esc(note)}"')
        out.append(f'<div class="vcf-row {state}" id="{side}{row.index}"'
                   f'{anchor} tabindex="0">'
                   f'<span class="tag">{state}</span>'
                   f'<span class="rec">{_esc(_vcf_text(row.text))}</span></div>')
    return "\n".join(out), states


def _viewer_section(truth_rows, pred_rows, meta):
    '''
    Truth and prediction side by side, each file whole, with the pairing wired up.

    Every aggregate above this reduces the two files to a number. This is the step before
    that: the records themselves, with a click on either side taking you to the record it
    was paired with. Where a record has no partner the other pane says so rather than
    leaving you to work out whether you missed it -- "there is nothing there" and "I cannot
    find it" look identical in a viewer that only scrolls.
    '''
    if not truth_rows and not pred_rows:
        return ""

    def truth_mate(locus):
        if locus.match is not None:
            return locus.match, False
        return locus.displaced, locus.displaced is not None

    def pred_mate(record):
        if record.matched_for is not None:
            return record.matched_for, False
        return record.displaced_for, record.displaced_for is not None

    truth_html, truth_states = _viewer_rows(
        truth_rows, "t", "p", truth_mate,
        lambda where: f"No prediction was paired with {where} — this locus was missed.")
    pred_html, pred_states = _viewer_rows(
        pred_rows, "p", "t", pred_mate,
        lambda where: f"No simulated locus was paired with the call at {where}.")

    def pane(side, title, path, rows, html, states):
        note = ""
        if len(rows) > MAX_VIEWER_ROWS:
            note = (f'<p class="note">Showing the first {MAX_VIEWER_ROWS:,} of '
                    f'{len(rows):,} records; the rest are in the VCF.</p>')
        return (
            f'<div class="pane-wrap">'
            f'<div class="pane-head"><strong>{_esc(title)}</strong> '
            f'<code>{_esc(os.path.basename(path))}</code> · {len(rows):,} records</div>'
            f'<div class="pane" id="pane-{side}">{html}</div>'
            f'<div class="bubble" id="bubble-{side}" hidden></div>'
            f'{note}</div>')

    legend = _legend([
        ("match — paired with the other side", "var(--series-1)"),
        ("displ — found, one solo LTR downstream", "var(--series-2)"),
        ("alone — nothing paired with it", "var(--series-3)"),
        ("skip — not scored", "var(--ink-muted)"),
    ])
    return (
        '<section class="card">'
        '<h2>Truth and prediction, side by side</h2>'
        '<p class="sub">Click a record to jump the other pane to the one it was paired '
        'with. Scrolling either pane carries the other to the same relative position, so '
        'two files of different lengths stay roughly aligned.</p>'
        f'<div class="viewer" id="vcf-viewer">'
        + pane("truth", "truth", meta["truth"], truth_rows, truth_html, truth_states)
        + pane("pred", "prediction", meta["pred"], pred_rows, pred_html, pred_states)
        + '</div>'
        + legend
        + f'<p class="note">Allele sequences longer than {ALLELE_PREVIEW_BP} bp are shown '
        'as their first bases and their length.</p>'
        '</section>'
    )


def _unsupported_section(summary):
    unsupported = summary.get("unsupported_events") or []
    if not unsupported:
        return ""
    rows = [[f'{item["chrom"]}:{item["pos"]}', item.get("id", "."),
             "; ".join(item["declared"])] for item in unsupported]
    return (
        '<section class="card">'
        f'<h2>Loci declaring an event no allele bears out ({len(unsupported)})</h2>'
        '<p class="sub">Not counted under that event type: the allele does not change '
        'length the way the event it is labelled with would.</p>'
        f'<div class="scroll">{_table(["locus", "id", "declared"], rows, numeric_from=3)}</div>'
        '</section>'
    )


def _tabs(has_viewer):
    '''
    The tab strip, hidden until the script that drives it runs.

    Marked hidden in the markup and unhidden by the script rather than the other way
    round: with no JavaScript the page stays one long document, where every anchor still
    resolves and nothing is stranded behind a control that cannot be operated.
    '''
    if not has_viewer:
        return ""
    return (
        '<nav class="tabs" role="tablist" id="tabs" hidden>'
        '<button type="button" role="tab" id="tab-report" aria-controls="panel-report" '
        'aria-selected="true">Report</button>'
        '<button type="button" role="tab" id="tab-vcf" aria-controls="panel-vcf" '
        'aria-selected="false">VCF viewer</button>'
        '</nav>'
    )


def render(meta, summary, loci, locus_dir=None, summary_file=None,
           truth_rows=(), pred_rows=(), locus_rows=()):
    '''Build the whole page as one string.'''
    labels = locus_labels(loci, locus_rows)
    viewer = _viewer_section(truth_rows, pred_rows, meta)
    report = "\n".join(part for part in (
        _headline(summary, meta),
        _counts_section(summary, meta),
        _breakpoint_section(loci, labels, meta),
        _length_section(loci, labels),
        _size_section(loci, labels),
        _strata_section(summary, meta),
        _unsupported_section(summary),
    ) if part)
    body = "\n".join(part for part in (
        _header(meta, bool(viewer)),
        _tabs(bool(viewer)),
        f'<div id="panel-report" role="tabpanel" aria-labelledby="tab-report">'
        f'{report}</div>',
        (f'<div id="panel-vcf" role="tabpanel" aria-labelledby="tab-vcf">{viewer}</div>'
         if viewer else ""),
    ) if part)
    trail = []
    if summary_file:
        trail.append(f"summary: <code>{_esc(os.path.basename(summary_file))}</code>")
    if locus_dir:
        trail.append(f"per-locus JSON: <code>{_esc(os.path.basename(locus_dir))}/</code>")
    footer = ("<footer>Written by tevarsim Evaluate. "
              + " &middot; ".join(trail) + "</footer>") if trail else ""
    return (
        "<!DOCTYPE html>\n"
        '<html lang="en">\n<head>\n<meta charset="utf-8">\n'
        '<meta name="viewport" content="width=device-width, initial-scale=1">\n'
        f"<title>tevarsim Evaluate — {_esc(os.path.basename(meta['pred']))}</title>\n"
        f"<style>{CSS}</style>\n</head>\n<body>\n<main>\n{body}\n{footer}\n</main>\n"
        f"<script>{JS}</script>\n</body>\n</html>\n"
    )


def write_report(path, meta, summary, loci, locus_dir=None, summary_file=None,
                 truth_rows=(), pred_rows=(), locus_rows=()):
    '''Render the report and write it to ``path``.'''
    with open(path, "w") as fo:
        fo.write(render(meta, summary, loci, locus_dir, summary_file,
                        truth_rows, pred_rows, locus_rows))
    return path
