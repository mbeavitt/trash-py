"""HOR threshold sweep — an interactive self-similarity dot-plot you can scrub
across divergence thresholds with a slider.

The one-shot :mod:`hor_plot` renders a single ``--threshold`` value. This module
re-runs the HOR scan (``_ext.find_hors``) at every integer threshold from 1% up
to a maximum, then bakes all those frames into one self-contained HTML page:
Plotly (loaded from CDN) draws the same magma-on-dark dot-plot as
:mod:`hor_plot`, and a slider swaps between the pre-rendered frames.

Why a re-scan per threshold rather than filtering one big scan: the threshold is
*inside* the scanner — a HOR block is extended only while consecutive repeat
pairs stay under the SNV cutoff, and it closes the moment one exceeds it. So a
lower threshold yields shorter, more-fragmented blocks, not merely a subset of a
higher-threshold run. Each frame therefore needs its own scan. The expensive
step (MAFFT) is done once upstream; these scans reuse that alignment and are
cheap, and identical ``threshold_SNV`` values (the ``floor`` collapses several
percentages onto the same integer) are scanned once and cached.

Self-comparison (method 1) only — that is where the dot-plot lives.
"""
from __future__ import annotations

import json
import math
import time
from dataclasses import dataclass
from pathlib import Path

from . import _log as log
from . import _ext

# Match hor_plot's "stained-glass" palette so the interactive view reads the same.
_BG = "#0d0d12"
_INK = "#c8c8d0"
_GRID = "#2a2a33"

# Matplotlib's "magma" colormap as explicit Plotly stops (dark purple -> magenta
# -> orange -> light yellow). Plotly.js has no built-in colorscale named "magma"
# — passing the string silently falls back to a default scale — so we ship the
# actual ramp to match hor_plot / plot_hor_table pixel-for-hue. Used with
# reversescale=true so conserved HORs (low SNV) glow yellow and divergent ones
# recede to purple, exactly like the static plots' magma(1 - t).
_MAGMA = [
    [0.0, "rgb(0,0,4)"],
    [0.0667, "rgb(11,9,36)"],
    [0.1333, "rgb(32,17,75)"],
    [0.2, "rgb(59,15,112)"],
    [0.2667, "rgb(87,21,126)"],
    [0.3333, "rgb(114,31,129)"],
    [0.4, "rgb(140,41,129)"],
    [0.4667, "rgb(168,50,125)"],
    [0.5333, "rgb(196,60,117)"],
    [0.6, "rgb(222,73,104)"],
    [0.6667, "rgb(241,96,93)"],
    [0.7333, "rgb(250,127,94)"],
    [0.8, "rgb(254,159,109)"],
    [0.8667, "rgb(254,191,132)"],
    [0.9333, "rgb(253,222,160)"],
    [1.0, "rgb(252,253,191)"],
]

# Points kept per threshold frame (before the y=x mirror doubles them). Caps the
# HTML size on huge arrays; a strided sample, like hor_plot's plot_cap. 30k keeps
# the dot-plot dense/faithful; on a ~1M-HOR centromeric array a 50-frame sweep is
# then ~34 MB — fine to open locally, though above some upload limits.
SWEEP_PLOT_CAP = 30_000
_PLOTLY_CDN = "https://cdn.plot.ly/plotly-2.35.2.min.js"


def _median(values) -> float:
    s = sorted(values)
    n = len(s)
    if not n:
        return 0.0
    mid = n // 2
    return float(s[mid]) if n % 2 else (s[mid - 1] + s[mid]) / 2.0


@dataclass
class SweepData:
    """The full result of a threshold sweep, in one structured object.

    Everything downstream (the 2D slider, the NPZ dump, the 3D stack, any custom
    replot) is derived from this. HOR records are stored once per *distinct*
    ``threshold_SNV`` — the ``floor`` collapses several percentages onto the same
    cutoff — with ``pct_snv`` mapping each threshold percent back to its scan.
    Records keep the raw unit indices (``sa..eb``, 1-based into the repeats), not
    bp: compact, and every bp/size/divergence column re-derives from the repeat
    coordinates with the exact upstream formulas."""
    starts: "object"          # np.int64[n_rep]  repeat start (bp)
    ends: "object"            # np.int64[n_rep]  repeat end (bp)
    widths: "object"          # np.int32[n_rep]  monomer width (bp)
    median_w: float
    min_len: int
    chrA: str
    hor_class: str
    thresholds: list          # [1..max_threshold]
    pct_snv: "object"         # np.int32[T]  threshold_SNV per percent
    scans: dict               # threshold_SNV -> np.int32[(n, 6)] : sa,ea,sb,eb,dir,tv
    unit_val: float = 1_000_000.0   # bp per plot unit (Mbp)

    def counts(self):
        import numpy as np
        return np.array([len(self.scans[int(s)]) for s in self.pct_snv], dtype=np.int64)


def scan_sweep(aligned: Path, repeats: list, max_threshold: int = 30,
               min_len: int = 3, workers: int = 1) -> SweepData:
    """Run the HOR scan at every threshold 1..max_threshold (%) on an existing
    alignment and return a :class:`SweepData`. Distinct ``threshold_SNV`` values
    are scanned once. With ``workers > 1``, independent cutoffs run concurrently;
    each native scan releases the GIL."""
    import numpy as np

    n_rep = len(repeats)
    starts = np.fromiter((r.start for r in repeats), np.int64, n_rep)
    ends = np.fromiter((r.end for r in repeats), np.int64, n_rep)
    widths = np.fromiter((r.width for r in repeats), np.int32, n_rep)
    median_w = _median([r.width for r in repeats])

    thresholds = list(range(1, max(1, max_threshold) + 1))
    pct_snv = np.array([math.floor(p * median_w / 100) for p in thresholds], np.int32)

    # Stream each distinct scan's raw int32 rows to a temp file and read them back
    # as an array, rather than materialising `_ext.find_hors`'s Python list of
    # tuples — at high thresholds a single scan can emit tens of millions of HORs,
    # where the list (≈120 B/HOR) would dwarf the int32 array (24 B/HOR).
    aligned = Path(aligned)
    distinct = list(dict.fromkeys(pct_snv.tolist()))
    scans: dict[int, "object"] = {}
    pcts_by_snv: dict[int, list[int]] = {}
    for pct, snv in zip(thresholds, pct_snv.tolist()):
        pcts_by_snv.setdefault(snv, []).append(pct)
    scan_number = {snv: i for i, snv in enumerate(distinct, start=1)}

    def scan_one(snv: int):
        pcts = pcts_by_snv[snv]
        pct_label = (f"{pcts[0]}%" if len(pcts) == 1
                     else f"{pcts[0]}–{pcts[-1]}%")
        label = f"HOR sweep scan {scan_number[snv]}/{len(distinct)}: {pct_label}, SNV cutoff {snv}"
        log.detail(f"{label} — starting")
        started = time.perf_counter()
        raw = aligned.with_suffix(f".sweep.{snv}.raw")
        try:
            n = _ext.find_hors_stream(
                str(aligned), str(raw), 1, int(snv), min_len, 1
            )
            rows = (np.fromfile(raw, dtype=np.int32).reshape(-1, 6) if n
                    else np.empty((0, 6), np.int32))
            log.detail(
                f"{label} — {n:,} HORs, {log.format_elapsed(time.perf_counter() - started)}"
            )
            return snv, rows
        finally:
            raw.unlink(missing_ok=True)

    if workers > 1 and len(distinct) > 1:
        from concurrent.futures import ThreadPoolExecutor
        with ThreadPoolExecutor(max_workers=min(workers, len(distinct))) as pool:
            for snv, rows in pool.map(scan_one, distinct):
                scans[snv] = rows
    else:
        for snv in distinct:
            key, rows = scan_one(snv)
            scans[key] = rows
    return SweepData(starts=starts, ends=ends, widths=widths, median_w=median_w,
                     min_len=min_len, chrA=chrA_of(repeats), hor_class="",
                     thresholds=thresholds, pct_snv=pct_snv, scans=scans)


def chrA_of(repeats: list) -> str:
    return repeats[0].seqID if repeats else ""


def _derive(sw: SweepData, snv: int):
    """Return (xa, yb, snv_per_kbp) in plot units for one scan's HOR records."""
    import numpy as np
    rows = sw.scans[int(snv)]
    if len(rows) == 0:
        return np.empty(0), np.empty(0), np.empty(0)
    sa, ea, sb, eb, tv = (rows[:, 0].astype(np.int64), rows[:, 1].astype(np.int64),
                          rows[:, 2].astype(np.int64), rows[:, 3].astype(np.int64),
                          rows[:, 5].astype(np.int64))
    sabp = sw.starts[sa - 1]; eabp = sw.ends[ea - 1]
    sbbp = sw.starts[sb - 1]; ebbp = sw.ends[eb - 1]
    ba = eabp - sabp + 1
    xa = (sabp + eabp) / 2.0 / sw.unit_val
    yb = (sbbp + ebbp) / 2.0 / sw.unit_val
    snv_kbp = 1000.0 * tv / ba
    return xa, yb, snv_kbp


def run_hor_sweep(aligned: Path, repeats: list, out_html: Path, chrA: str,
                  hor_class: str, max_threshold: int = 30, min_len: int = 3,
                  plot_cap: int = SWEEP_PLOT_CAP, dump: bool = True,
                  make_3d: bool = True, workers: int = 1) -> Path | None:
    """Scan the aligned repeats across thresholds 1..max_threshold and write the
    outputs: the interactive 2D slider dot-plot (``out_html``), and — when
    enabled — the structured NPZ dump (``*.npz``) and the 3D stacked view
    (``*_3d.html``) alongside it.

    Returns the 2D HTML path, or ``None`` if no HORs were found at any threshold."""
    sw = scan_sweep(aligned, repeats, max_threshold=max_threshold, min_len=min_len,
                    workers=workers)
    sw.chrA = chrA
    sw.hor_class = hor_class
    if int(sw.counts().sum()) == 0:
        return None

    def write_phase(label: str, path: Path, writer) -> None:
        log.detail(f"{label}...")
        started = time.perf_counter()
        writer()
        size_mb = path.stat().st_size / 1_000_000
        log.detail(
            f"wrote {path.name} ({size_mb:,.1f} MB), "
            f"{log.format_elapsed(time.perf_counter() - started)}"
        )

    out_html = Path(out_html)
    write_phase(
        "building 2D sweep HTML", out_html,
        lambda: out_html.write_text(_render_html(_build_2d_payload(sw, plot_cap))),
    )
    if dump:
        npz = out_html.with_suffix(".npz")
        write_phase("writing compressed sweep NPZ", npz,
                    lambda: write_sweep_npz(sw, npz))
    if make_3d:
        three_d = out_html.with_name(out_html.stem + "_3d.html")
        write_phase(
            "building 3D sweep HTML", three_d,
            lambda: three_d.write_text(_render_3d_html(_build_3d_payload(sw))),
        )
    return out_html


def _build_2d_payload(sw: SweepData, plot_cap: int) -> dict:
    import numpy as np
    frames: list[dict] = []
    lo = hi = None
    for pct, snv in zip(sw.thresholds, sw.pct_snv.tolist()):
        xa, yb, snv_kbp = _derive(sw, snv)
        n_total = len(xa)
        if plot_cap and n_total > plot_cap:                   # strided subsample
            step = n_total // plot_cap or 1
            xa, yb, snv_kbp = xa[::step], yb[::step], snv_kbp[::step]
        if len(xa):                                           # shared square window
            flo = float(min(xa.min(), yb.min())); fhi = float(max(xa.max(), yb.max()))
            lo = flo if lo is None else min(lo, flo)
            hi = fhi if hi is None else max(hi, fhi)
        frames.append({
            "pct": pct, "snv": int(snv), "n": int(n_total),
            # Per-frame colour ceiling, exactly as the static plot_hors: vmax =
            # threshold%*10 — spans the full magma range at each threshold.
            "vmax": max(pct * 10.0, 1e-9),
            "xa": [round(v, 4) for v in xa.tolist()],
            "yb": [round(v, 4) for v in yb.tolist()],
            "c": [round(v, 2) for v in snv_kbp.tolist()],
        })
    if lo is None or hi is None or hi <= lo:
        lo, hi = 0.0, 1.0
    pad = 0.02 * (hi - lo) or 1.0
    return {
        "frames": frames, "lo": lo - pad, "hi": hi + pad,
        "chrA": sw.chrA, "cls": sw.hor_class,
        "bg": _BG, "ink": _INK, "grid": _GRID, "cs": _MAGMA,
    }


def _render_html(payload: dict) -> str:
    data_json = json.dumps(payload, separators=(",", ":"))
    title = f"HOR sweep — {payload['cls']}, {payload['chrA']}"
    return _HTML_TEMPLATE.replace("__TITLE__", title) \
                         .replace("__CDN__", _PLOTLY_CDN) \
                         .replace("__BG__", payload["bg"]) \
                         .replace("__DATA__", data_json)


# The page is intentionally dependency-light: one <script> from the Plotly CDN
# plus the baked-in frame data. Frames are stored un-mirrored (half the bytes);
# the y=x mirror that makes the self-comparison symmetric is applied in JS.
_HTML_TEMPLATE = """<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>__TITLE__</title>
<script src="__CDN__" charset="utf-8"></script>
<style>
  html, body { margin: 0; height: 100%; background: __BG__; }
  #plot { width: 100vw; height: 100vh; }
  #err { color: #e0a0a0; font-family: sans-serif; padding: 1rem; }
</style>
</head>
<body>
<div id="plot"></div>
<div id="err"></div>
<script>
const D = __DATA__;

function mirror(f) {
  // Symmetric self-dot-plot: one dot per HOR (block A on x, block B on y) plus
  // its reflection across y = x. Divergence colour is carried on both halves.
  return { x: f.xa.concat(f.yb), y: f.yb.concat(f.xa), marker: { color: f.c.concat(f.c) } };
}

function countAnno(f) {
  return [{
    xref: 'paper', yref: 'paper', x: 0.02, y: 0.98, xanchor: 'left', yanchor: 'top',
    showarrow: false, align: 'left',
    text: f.n.toLocaleString() + ' HORs &nbsp;·&nbsp; SNV cutoff ' + f.snv,
    font: { color: D.ink, size: 12 }, bgcolor: 'rgba(13,13,18,0.6)'
  }];
}

try {
  // Pre-mirror every frame once (block A/B swap across y = x).
  const M = D.frames.map(mirror);
  const f0 = D.frames[0], m0 = M[0];
  const trace = {
    type: 'scattergl', mode: 'markers', x: m0.x, y: m0.y,
    hovertemplate: 'A %{x:.3f} Mbp<br>B %{y:.3f} Mbp<br>%{marker.color:.1f} SNV/kbp<extra></extra>',
    marker: {
      size: 4, opacity: 0.5, color: m0.marker.color,
      colorscale: D.cs, reversescale: true, cmin: 0, cmax: f0.vmax,
      colorbar: {
        title: { text: 'SNV per kbp', font: { color: D.ink } },
        tickfont: { color: D.ink }, outlinecolor: D.grid, thickness: 14
      }
    }
  };

  // Drive the slider with method:'update' (a direct restyle+relayout of the one
  // trace) rather than animate/frames: scattergl traces are known to drop their
  // points under frame animation, but a plain update redraws them reliably. Each
  // step also rescales the colour ceiling to that threshold's vmax (= pct*10) so
  // the full magma range — up to purple — is used at every threshold.
  const steps = D.frames.map((f, i) => ({
    method: 'update', label: String(f.pct),
    args: [
      { x: [M[i].x], y: [M[i].y], 'marker.color': [M[i].marker.color], 'marker.cmax': [f.vmax] },
      { annotations: countAnno(f) },
      [0]
    ]
  }));

  const layout = {
    paper_bgcolor: D.bg, plot_bgcolor: D.bg, font: { color: D.ink },
    title: { text: 'HOR self-similarity — ' + D.cls + ', ' + D.chrA, font: { color: D.ink } },
    margin: { l: 70, r: 20, t: 60, b: 110 },
    xaxis: { range: [D.lo, D.hi], title: { text: 'HOR block A position — ' + D.chrA + ' (Mbp)' },
             gridcolor: D.grid, zeroline: false, color: D.ink, constrain: 'domain' },
    yaxis: { range: [D.lo, D.hi], title: { text: 'HOR block B position — ' + D.chrA + ' (Mbp)' },
             gridcolor: D.grid, zeroline: false, color: D.ink, scaleanchor: 'x', scaleratio: 1 },
    shapes: [{ type: 'line', x0: D.lo, y0: D.lo, x1: D.hi, y1: D.hi, line: { color: D.grid, width: 1 } }],
    annotations: countAnno(f0),
    sliders: [{
      active: 0, x: 0.5, xanchor: 'center', y: 0, yanchor: 'top', len: 0.9, pad: { t: 20, b: 10 },
      currentvalue: { prefix: 'threshold = ', suffix: ' %', font: { color: D.ink, size: 14 } },
      font: { color: D.ink }, bgcolor: D.grid, activebgcolor: D.ink,
      tickcolor: D.grid, steps: steps
    }]
  };

  Plotly.newPlot('plot', [trace], layout, { responsive: true, displaylogo: false });
} catch (e) {
  document.getElementById('err').textContent =
    'Failed to render (is the Plotly CDN reachable?): ' + e;
}
</script>
</body>
</html>
"""


# --------------------------------------------------------------------------
# Structured dump (NPZ) + loader
# --------------------------------------------------------------------------

def write_sweep_npz(sw: SweepData, path: Path) -> Path:
    """Persist the whole sweep to a compressed ``.npz`` (numpy-native, no extra
    deps). HOR records for the distinct scans are concatenated with an offsets
    array; ``load_sweep`` reconstitutes per-threshold views. Round-trips exactly
    — every derived column recomputes from these arrays."""
    import numpy as np

    snv_values = sorted(set(int(s) for s in sw.pct_snv.tolist()))
    blocks = [sw.scans[s] for s in snv_values]                # each (n, 6) int32
    offsets = np.zeros(len(snv_values) + 1, np.int64)
    offsets[1:] = np.cumsum([len(b) for b in blocks])
    rec = (np.concatenate(blocks, axis=0) if blocks
           else np.empty((0, 6), np.int32))                   # sa,ea,sb,eb,dir,tv

    path = Path(path)
    np.savez_compressed(
        path,
        rep_start=sw.starts, rep_end=sw.ends, rep_width=sw.widths,
        pct=np.array(sw.thresholds, np.int16),
        pct_snv=sw.pct_snv.astype(np.int32),
        snv_values=np.array(snv_values, np.int32),
        scan_offsets=offsets,
        sa=rec[:, 0], ea=rec[:, 1], sb=rec[:, 2], eb=rec[:, 3],
        direction=rec[:, 4].astype(np.int8), total_variant=rec[:, 5],
        median_w=np.float64(sw.median_w), min_len=np.int32(sw.min_len),
        unit_val=np.float64(sw.unit_val),
        chrA=np.array(sw.chrA), hor_class=np.array(sw.hor_class),
    )
    # np.savez appends .npz if absent — return the real path either way.
    return path if path.suffix == ".npz" else path.with_suffix(".npz")


class SweepDump:
    """Loaded view over a sweep ``.npz``. Hands you numpy arrays (or a pandas
    DataFrame) per threshold plus the whole-sweep summary, for arbitrary replots.

    >>> s = load_sweep("HORs_sweep_178_1_CP116280.1.npz")
    >>> s.thresholds                 # [1, 2, ..., 50]
    >>> rec = s.records(25)          # dict of arrays for threshold 25%
    >>> rec["snv_per_kbp"], rec["block_A_centre_mbp"]
    >>> s.summary()                  # per-threshold count / divergence / block size
    >>> s.dotplot(25); s.count_curve()   # quick matplotlib recipes
    """

    def __init__(self, npz: dict):
        import numpy as np
        self._np = np
        self.rep_start = npz["rep_start"]; self.rep_end = npz["rep_end"]
        self.rep_width = npz["rep_width"]
        self.thresholds = npz["pct"].astype(int).tolist()
        self.pct_snv = npz["pct_snv"]
        self.snv_values = npz["snv_values"]
        self._off = npz["scan_offsets"]
        self._sa = npz["sa"]; self._ea = npz["ea"]; self._sb = npz["sb"]
        self._eb = npz["eb"]; self._dir = npz["direction"]; self._tv = npz["total_variant"]
        self.median_w = float(npz["median_w"]); self.min_len = int(npz["min_len"])
        self.unit_val = float(npz["unit_val"])
        self.chrA = str(npz["chrA"]); self.hor_class = str(npz["hor_class"])

    def snv_for(self, pct: int) -> int:
        return int(self.pct_snv[self.thresholds.index(pct)])

    def _slice(self, pct: int):
        snv = self.snv_for(pct)
        i = int(self._np.searchsorted(self.snv_values, snv))
        a, b = int(self._off[i]), int(self._off[i + 1])
        return a, b

    def records(self, pct: int) -> dict:
        """All HORs at threshold ``pct`` (%), with derived columns, as a dict of
        numpy arrays. Positions are in bp; centres also given in plot units (Mbp
        by default). ``snv_per_kbp`` and block sizes use the exact upstream
        formulas (SNV/kbp divides by block A size, matching the HOR table)."""
        np = self._np
        a, b = self._slice(pct)
        sa = self._sa[a:b].astype(np.int64); ea = self._ea[a:b].astype(np.int64)
        sb = self._sb[a:b].astype(np.int64); eb = self._eb[a:b].astype(np.int64)
        tv = self._tv[a:b].astype(np.int64)
        sabp = self.rep_start[sa - 1]; eabp = self.rep_end[ea - 1]
        sbbp = self.rep_start[sb - 1]; ebbp = self.rep_end[eb - 1]
        ba = eabp - sabp + 1; bb = ebbp - sbbp + 1
        u = self.unit_val
        return {
            "start_A_bp": sabp, "end_A_bp": eabp, "start_B_bp": sbbp, "end_B_bp": ebbp,
            "block_size_units": ea - sa + 1,
            "block_A_size_bp": ba, "block_B_size_bp": bb,
            "direction": self._dir[a:b],           # 1 = tandem/parallel, 2 = inverted
            "total_variant": tv,
            "snv_per_kbp": 1000.0 * tv / ba,
            "block_A_centre_mbp": (sabp + eabp) / 2.0 / u,
            "block_B_centre_mbp": (sbbp + ebbp) / 2.0 / u,
        }

    def frame(self, pct: int):
        """`records(pct)` as a pandas DataFrame (requires pandas)."""
        import pandas as pd
        return pd.DataFrame(self.records(pct))

    def counts(self):
        return self._np.array([self._off[int(self._np.searchsorted(self.snv_values,
                               self.snv_for(p))) + 1]
                               - self._off[int(self._np.searchsorted(self.snv_values,
                               self.snv_for(p)))] for p in self.thresholds],
                              dtype=self._np.int64)

    def summary(self):
        """Per-threshold summary table (dict of arrays; DataFrame if pandas is
        present): HOR count, mean/median SNV/kbp, mean block size in units."""
        np = self._np
        pct = np.array(self.thresholds)
        cnt, msnv, medsnv, mblk = [], [], [], []
        for p in self.thresholds:
            r = self.records(p); s = r["snv_per_kbp"]
            cnt.append(len(s))
            msnv.append(float(s.mean()) if len(s) else 0.0)
            medsnv.append(float(np.median(s)) if len(s) else 0.0)
            mblk.append(float(r["block_size_units"].mean()) if len(s) else 0.0)
        out = {"threshold_pct": pct, "threshold_snv": self.pct_snv.copy(),
               "n_hors": np.array(cnt), "mean_snv_per_kbp": np.array(msnv),
               "median_snv_per_kbp": np.array(medsnv), "mean_block_units": np.array(mblk)}
        try:
            import pandas as pd
            return pd.DataFrame(out)
        except Exception:
            return out

    # -- quick matplotlib recipes (best-effort; import lazily) ---------------
    def dotplot(self, pct: int, ax=None):
        """Static self-similarity dot-plot at one threshold, magma-on-dark, the
        same look as ``hor_plot`` — a convenience for replotting from the dump."""
        import numpy as np
        import matplotlib.pyplot as plt
        r = self.records(pct)
        xa, yb, snv = r["block_A_centre_mbp"], r["block_B_centre_mbp"], r["snv_per_kbp"]
        xs = np.concatenate([xa, yb]); ys = np.concatenate([yb, xa])
        vmax = max(pct * 10.0, 1e-9)
        shade = np.clip(np.concatenate([snv, snv]) / vmax, 0, 1)
        colors = plt.get_cmap("magma")(1.0 - shade); colors[:, 3] = 0.5
        if ax is None:
            _, ax = plt.subplots(figsize=(8, 8), dpi=160)
        ax.set_facecolor(_BG)
        ax.scatter(xs, ys, s=3, c=colors, linewidths=0)
        ax.set_aspect("equal")
        ax.set_title(f"{self.hor_class}, {self.chrA} — threshold {pct}%", color=_INK)
        ax.set_xlabel("block A (Mbp)"); ax.set_ylabel("block B (Mbp)")
        return ax

    def count_curve(self, ax=None):
        """HOR count vs threshold — the sweep's headline curve."""
        import matplotlib.pyplot as plt
        if ax is None:
            _, ax = plt.subplots(figsize=(7, 4), dpi=160)
        ax.plot(self.thresholds, self.counts(), marker="o", ms=3, color="#e0578d")
        ax.set_xlabel("threshold (%)"); ax.set_ylabel("HORs identified")
        ax.set_title(f"HOR count vs threshold — {self.hor_class}, {self.chrA}")
        return ax


def load_sweep(path: Path) -> SweepDump:
    """Load a sweep ``.npz`` written by :func:`write_sweep_npz`."""
    import numpy as np
    with np.load(Path(path)) as z:
        return SweepDump({k: z[k] for k in z.files})


# --------------------------------------------------------------------------
# 3D stack — the sweep as a physical structure (one dot-plot layer per threshold)
# --------------------------------------------------------------------------

# Points kept per threshold *layer* in the 3D view (before mirroring). 3D scatter
# is far heavier than 2D, so this is much smaller than SWEEP_PLOT_CAP. Kept lean
# because the markers are opaque (for correct depth) — too many would wall off
# the interior of the stack.
SWEEP_3D_CAP = 1_400


def _build_3d_payload(sw: SweepData, cap: int = SWEEP_3D_CAP) -> dict:
    """Stack every threshold as a z-layer. z = pct and the z-axis is reversed, so
    the highest threshold sits at the bottom of the pile and t=1 rests on top —
    the layers read as slices through one 3D structure. Colour is divergence
    normalised to each layer's own ceiling (snv / (pct*10)), so the full magma
    range shows at every depth."""
    import numpy as np
    X, Y, Z, C = [], [], [], []
    lo = hi = None
    for pct, snv in zip(sw.thresholds, sw.pct_snv.tolist()):
        xa, yb, snv_kbp = _derive(sw, snv)
        n = len(xa)
        if n == 0:
            continue
        if cap and n > cap:
            step = n // cap or 1
            xa, yb, snv_kbp = xa[::step], yb[::step], snv_kbp[::step]
        shade = np.clip(snv_kbp / max(pct * 10.0, 1e-9), 0.0, 1.0)
        # mirror across the A/B diagonal so each layer is the symmetric dot-plot
        xs = np.concatenate([xa, yb]); ys = np.concatenate([yb, xa])
        cs = np.concatenate([shade, shade])
        X.append(xs); Y.append(ys)
        Z.append(np.full(len(xs), float(pct)))
        C.append(cs)
        flo = float(min(xs.min(), ys.min())); fhi = float(max(xs.max(), ys.max()))
        lo = flo if lo is None else min(lo, flo)
        hi = fhi if hi is None else max(hi, fhi)
    if not X:
        lo, hi = 0.0, 1.0
        x = y = z = c = []
    else:
        x = np.concatenate(X); y = np.concatenate(Y)
        z = np.concatenate(Z); c = np.concatenate(C)
        x = [round(v, 4) for v in x.tolist()]
        y = [round(v, 4) for v in y.tolist()]
        z = [int(v) for v in z.tolist()]
        c = [round(v, 3) for v in c.tolist()]
    if lo is None or hi is None or hi <= lo:
        lo, hi = 0.0, 1.0
    pad = 0.02 * (hi - lo) or 1.0
    return {
        "x": x, "y": y, "z": z, "c": c,
        "lo": lo - pad, "hi": hi + pad,
        "zmin": min(sw.thresholds), "zmax": max(sw.thresholds),
        "chrA": sw.chrA, "cls": sw.hor_class,
        "bg": _BG, "ink": _INK, "grid": _GRID, "cs": _MAGMA,
    }


def _render_3d_html(payload: dict) -> str:
    data_json = json.dumps(payload, separators=(",", ":"))
    title = f"HOR sweep 3D — {payload['cls']}, {payload['chrA']}"
    return _HTML_3D_TEMPLATE.replace("__TITLE__", title) \
                            .replace("__CDN__", _PLOTLY_CDN) \
                            .replace("__BG__", payload["bg"]) \
                            .replace("__DATA__", data_json)


_HTML_3D_TEMPLATE = """<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>__TITLE__</title>
<script src="__CDN__" charset="utf-8"></script>
<style>
  html, body { margin: 0; height: 100%; background: __BG__; }
  #plot { width: 100vw; height: 100vh; }
  #err { color: #e0a0a0; font-family: sans-serif; padding: 1rem; }
</style>
</head>
<body>
<div id="plot"></div>
<div id="err"></div>
<script>
const D = __DATA__;
try {
  const trace = {
    type: 'scatter3d', mode: 'markers', x: D.x, y: D.y, z: D.z,
    hovertemplate: 'A %{x:.3f} Mbp<br>B %{y:.3f} Mbp<br>threshold %{z}%<extra></extra>',
    marker: {
      // Fully opaque on purpose: plotly scatter3d skips depth-writing for
      // opacity < 1, so transparent points draw in insertion order and the far
      // side bleeds over the near side. Opaque markers restore real occlusion.
      size: 2.2, opacity: 1.0, color: D.c,
      colorscale: D.cs, reversescale: true, cmin: 0, cmax: 1,
      colorbar: {
        title: { text: 'divergence<br>(fraction of<br>threshold)', font: { color: D.ink } },
        tickfont: { color: D.ink }, outlinecolor: D.grid, thickness: 14
      }
    }
  };
  const ax = (t) => ({ title: { text: t }, color: D.ink, gridcolor: D.grid,
                       zerolinecolor: D.grid, backgroundcolor: D.bg, showbackground: true });
  const layout = {
    paper_bgcolor: D.bg, font: { color: D.ink },
    title: { text: 'HOR sweep stack — ' + D.cls + ', ' + D.chrA
                   + '  (t=' + D.zmax + ' bottom → t=1 top)', font: { color: D.ink } },
    margin: { l: 0, r: 0, t: 50, b: 0 },
    scene: {
      xaxis: ax('block A (Mbp)'), yaxis: ax('block B (Mbp)'),
      // reversed z so the highest threshold layer sits at the bottom of the pile.
      zaxis: Object.assign(ax('threshold (%)'), { autorange: 'reversed' }),
      aspectmode: 'manual', aspectratio: { x: 1, y: 1, z: 0.55 },
      camera: { eye: { x: 1.5, y: 1.5, z: 0.6 } }
    }
  };
  Plotly.newPlot('plot', [trace], layout, { responsive: true, displaylogo: false });
} catch (e) {
  document.getElementById('err').textContent =
    'Failed to render (is the Plotly CDN reachable?): ' + e;
}
</script>
</body>
</html>
"""
