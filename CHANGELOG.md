# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [2.6.0] - 2026-07-22

### Added

- **Read the input fasta from stdin with `-f -`.** Resolves to `/dev/stdin`, so
  trash-py can sit in a pipeline (`samtools faidx ref.fa chr1 | trash-py -f -
  -o out`). Gzip is now detected by magic bytes instead of the `.gz` suffix,
  which a pipe doesn't have, so `cat x.fa.gz | trash-py -f -` works as well;
  on-disk `.gz` inputs are unaffected. POSIX only — `/dev/stdin` does not exist
  on Windows.

- **`-n` / `--name` sets the output filename prefix.** Defaults to the input
  filename as before (which is `stdin` when reading from a pipe). Mainly useful
  for piped runs, where the default prefix carries no information and repeated
  runs would otherwise overwrite each other.

### Changed

- **The HOR divergence threshold now defaults to 4%, down from 25%.** This
  affects `--hor-threshold` on the main pipeline, `-t/--threshold` on
  `trash-py hor`, and the colour scale in `scripts/plot_hor_table.py`. The
  default is now a stringent HOR definition rather than upstream HORT.R's
  permissive one.

  **This changes output for anyone who does not pass the flag explicitly.** The
  threshold is a percentage of the median repeat width, so for a 178 bp monomer
  it moves `threshold_SNV` from 44 to 7; on the test sample that takes the HOR
  count from 77 to 58. Pass `-t 25` (or `--hor-threshold 25`) to restore the
  previous behaviour.

  Golden-file and upstream-parity tests now pin `hor_threshold=25` explicitly,
  since the reference outputs they compare against were generated at that
  threshold.

## [2.5.1] - 2026-07-16

### Fixed

- **More informative HOR abort errors.** When HOR identification is aborted because requested sequence IDs are missing from the repeats table, the tool now prints an explicit `ERROR:` message stating that it is skipping all HOR identification, specifies the repeat class (species of HOR) being trialled, and lists the missing IDs.

## [2.5.0] - 2026-07-09

### Changed

- **CLI clarity pass — no flags removed, everything still parses.** The
  higher-order-repeat surface was confusing for new users; this cleans it up
  without breaking the HORT.R drag-and-drop contract.
  - The `hor` subcommand is now **listed under `trash-py --help`** (previously it
    was dispatched invisibly), and `trash-py hor` with no arguments prints its
    help instead of a terse one-liner.
  - `trash-py hor --help` is **decluttered**: it shows a compact canonical set and
    hides the HORT.R-compat aliases (`-r`, `--output_folder`, `-m/--method`,
    `-s/--saveR`, `-p/--plot-simple`, lowercase `--chrA`/`--chrB`, underscore
    long-names). All of them still parse; they're just no longer in the listing.
  - New canonical spellings **`-t/--threshold`** and **`-l/--min-len`** (the older
    `--hor-threshold`/`--hor-min-len`/underscore forms remain as hidden aliases).
  - The main pipeline's HOR options keep their `--hor-*` prefix, and the help now
    frames them explicitly as **a separate second stage** that runs after array
    identification — so a `--hor-*` option is never mistaken for a setting on the
    tandem-repeat step.
  - Redundant/duplicate aliases (e.g. `--chrA` vs `--ChrA`) now emit a one-line
    **deprecation nudge** on stderr but continue to work. Genuine HORT.R flags stay
    silent.
  - Error messages from the shared HOR runner now name the flags for the context
    they were reached from (`--ChrA`/`--ChrB` for the subcommand,
    `--hor-ChrA`/`--hor-ChrB` for the pipeline).

## [2.4.0] - 2026-07-08

### Performance

- **HOR detection now streams its output** instead of materialising the whole
  result in memory. The native scan writes raw HOR rows to disk
  (`_ext.find_hors_stream`) and the Python layer annotates + formats them in
  chunks, so **peak memory is O(N repeats), not O(#HORs)**. A single homogeneous
  array can produce tens of millions of HORs (10 GB+ tables) — previously that
  needed tens of GB of RAM (and could OOM); now it doesn't. Output is
  byte-identical; `hors_formed_count` is computed via an O(N) difference array
  and plots use a strided subsample.
- The per-pair SNV comparison **early-exits** once the divergence threshold is
  exceeded — byte-identical (the exact count is only ever used for pairs that
  pass), and faster on the O(N²) scan (helps most when many pairs are divergent).

## [2.3.2] - 2026-07-08

### Changed

- HOR **auto-class selection** now ranks classes by total repeat bp (Σ monomer
  width ≈ genomic coverage) instead of raw monomer count. A short, high-count
  microsatellite (e.g. a 13 bp repeat) no longer outranks the genome's actual
  major satellite. The log line reports coverage, e.g. `auto-selected class
  '43_28' (0.76 Mb of repeats, 18146 monomers, …); next by coverage: …`.
  (A. thaliana still auto-selects `178_1`.)

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
[2.5.1]: https://github.com/mbeavitt/trash-py/releases/tag/v2.5.1
[2.6.0]: https://github.com/mbeavitt/trash-py/releases/tag/v2.6.0
