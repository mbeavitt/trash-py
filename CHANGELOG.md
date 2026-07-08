# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [2.3.1] - 2026-07-08

### Changed

- HOR detection now validates the requested sequence IDs against the repeats
  table up front and **aborts (running nothing) if any are missing**, printing a
  truncated list (e.g. `not running HOR identification for 3 ID(s) not in the
  repeats table: Chr3, Chr4, Chr5`) — so requested sequences can't be silently
  skipped.

## [2.3.0] - 2026-07-08

### Added

- **`scripts/plot_hor_table.py`** — a standalone script that renders a
  chromosome's HOR dot-plot directly from its `HORs_<class>_<seq>.csv` table, so
  a plot can be regenerated or tweaked without re-running detection.
  `--label` sets the display name (e.g. `Chr3`).

### Changed

- **HOR plot reworked into a per-HOR dot-plot.** One translucent dot per HOR at
  the centres of its paired blocks (block A on x, block B on y), mirrored across
  y=x — replacing the block-to-block line segments. Tandem HORs land near the
  diagonal, inversions on the anti-diagonal, coloured by divergence (magma) on a
  dark ground. Axes are now labelled "HOR block A/B position", and the render is
  2× the previous resolution.

## [2.2.0] - 2026-07-08

### Added

- **HORT.R drop-in compatibility.** `trash-py hor` now accepts the full upstream
  `HORT.R` getopt surface — `-r/--repeats` (the table can be a flag as well as
  the positional arg), `-m/--method`, `-s/--saveR`, `-p/--plot_simple`,
  `-A/--chrA`, `-B/--chrB`, and the underscore long-names (`--output_folder`,
  `--hor_threshold`, …) — so replacing `Rscript HORT.R` with `trash-py hor` in an
  existing script works unchanged.
- **Multi-threaded MAFFT.** MAFFT is ~65–70% of HOR runtime and was
  single-threaded; it now runs `--thread`-parallel (`-T/--threads` for
  `trash-py hor`, default all cores; wired to `-p` in a pipeline run). Verified
  byte-identical to the single-threaded alignment, so the HOR tables are
  unchanged — it's a free ~4× on the alignment step (≈40s → ≈10s at 8 threads).

### Changed

- **HOR plot is now a "stained-glass" self-similarity dot-plot**
  (`HORs_dotplot_<class>_<seq>.png`): each HOR block is a translucent diagonal
  segment on a chromosome-vs-chromosome plane (parallel along the diagonal,
  inverted repeats along the anti-diagonal, mirrored across y=x), coloured by
  divergence on a perceptual sequential ramp over a dark ground. Replaces the old
  line plot.

## [2.1.0] - 2026-07-08

### Changed

- **HOR detection now has a dedicated `trash-py hor` subcommand** for running on
  an existing `…_repeats_with_seq.csv` (e.g. to re-run for a different class or
  set of sequences without redoing the pipeline):
  `trash-py hor repeats_with_seq.csv --chr-list Chr1,Chr2 -o out`.
- The in-pipeline "run and done" path is retained: the trigger flag is renamed
  `--chr-list` → **`--hor-chr-list`** (and the other in-pipeline HOR options are
  now `--hor-`prefixed: `--hor-ChrA`/`--hor-ChrB`/`--hor-class`/…) so they read
  clearly alongside the pipeline flags. `trash-py -f g.fasta -o out --hor-chr-list Chr1`
  runs the pipeline and then HOR detection on its output.

### Added

- **Automatic class selection.** `--class`/`--hor-class` is now optional; when
  omitted, the most abundant repeat class on the target sequence(s) is selected
  (and logged). This makes HOR detection usable on a new organism before you
  know which satellite is the major one.

## [2.0.0] - 2026-07-07

### Added

- **HOR (higher-order repeat) detection** — a port of the upstream `HORT.R`
  module, run as an optional stage after the main pipeline. Enable it with
  `--chr-list seq1,seq2,...` (self-comparison for several sequences — the common
  case) or `--ChrA seq` (single sequence), optionally with `--ChrB seq` for a
  cross-region comparison, plus `-c/--class`, `--hor-threshold` (default 25) and
  `--hor-min-len` (default 3).
- The HOR scan is a native reimplementation (`trash_py._ext.find_hors`) of the
  legacy `HOR.V3.3` binary, reverse-engineered from its disassembly. It
  reproduces the reference HOR tables essentially **byte-for-byte** across all
  five *A. thaliana* chromosomes (up to 1.45M rows each) — including the upstream
  quirks catalogued in `docs/HOR_source_bugs.md` (edge-close SNV under-count, SNV
  carry-over, the `SNV_per_kbp` block-A-twice typo, and more). The only
  differences are R's scientific-notation contractions (e.g. `2e+05`), which we
  deliberately keep as full integers, and a handful (~5 in ~3.9M rows) of
  last-digit float-rounding differences in `SNV_per_kbp` where R's and Python's
  formatters round a tie in opposite directions.
- Per-sequence outputs `HORs_<class>_<seq>.csv`,
  `repeats_with_hors_<class>_<seq>.csv`, and a matplotlib dot-plot
  `HORs_lines_<class>_<seq>.png`. The R `pretty()`/`hist()` binning used for the
  `start.adjusted` column is reproduced exactly.
- Optional `plot` extra (`pip install trash-py[plot]`) for the HOR plots; HOR
  detection additionally requires MAFFT on `PATH`.

## [1.2.0] - 2026-05-25

### Added

- **Programmatic annotation API.** New `trash_py.reads` module exposes
  `annotate()` and `annotate_sequences()` — thin wrappers around the
  pipeline that accept in-memory sequence strings and return typed
  dataclasses (`RepeatUnit`, `RepeatArray`, `SequenceAnnotation`)
  instead of writing CSV files. Useful for annotating individual long
  reads or insertion sequences inside larger Python workflows without
  going through temporary FASTA files. Input sequences shorter than
  `window_size` are rejected up front with a clear `ValueError`.
- **Built-in organism templates.** New `trash_py.templates` module
  ships `ARABIDOPSIS_CEN178`, the 177 bp Arabidopsis CEN178
  centromeric satellite consensus. Passing a template to `annotate()`
  anchors the repeat units to a fixed canonical rotation so breakpoint
  positions are directly comparable across reads.
- Top-level package now re-exports `annotate`, `annotate_sequences`,
  `RepeatUnit`, `RepeatArray`, `SequenceAnnotation`, and
  `ARABIDOPSIS_CEN178` for ergonomic imports.

[1.2.0]: https://github.com/mbeavitt/trash-py/releases/tag/v1.2.0

## [1.1.1] - 2026-05-21

### Added

- `-V` / `--version` flag — prints the installed version and exits.

[1.1.1]: https://github.com/mbeavitt/trash-py/releases/tag/v1.1.1

## [1.1.0] - 2026-05-21

### Added

- **Parallel processing.** The new `-p` / `--processes` flag runs the
  array-identification and repeat-mapping stages across multiple worker
  processes (default `1` = serial). This supersedes the v1.0.0 advice to
  parallelise externally with GNU Parallel — `trash-py` now handles
  per-sequence parallelism internally.

### Fixed

- The startup log header now reports the correct version. Previously
  `__version__` was pinned at `0.1.0` regardless of release.

[1.1.0]: https://github.com/mbeavitt/trash-py/releases/tag/v1.1.0

## [1.0.0] - 2026-04-28

Initial stable release of `trash-py` — a Python port of the
[TRASH](https://github.com/vlothec/TRASH_2) program (originally written in R)
for tandem-repeat array identification. Algorithmic credit for the
underlying approach belongs to Piotr Włodzimierz and the original TRASH
authors; this is an independent rewrite that preserves the structure and
logic of the upstream pipeline.

### What's different from upstream TRASH

- A more flexible Python/C runner-and-libraries framework: reusable
  repeat-annotation functions are exposed for incorporation into
  workflows beyond just running the bundled CLI.
- Hotter functions have been ported to C to maximise performance. On
  smaller, less repetitive genomes (e.g. Arabidopsis, human) the
  bottleneck is now nhmmer rather than the TRASH pipeline itself, which
  takes a fraction of the original time to complete.

### CLI differences from upstream TRASH

- `-q` is now available to silence logs.
- `-p` (multiprocessing) is **not** currently supported. This may be
  implemented in future. In the meantime, it's recommended to split
  input `.fasta` files and parallelise by chromosome/sequence using an
  external tool like GNU Parallel, and merge later.

[1.0.0]: https://github.com/mbeavitt/trash-py/releases/tag/v1.0.0
