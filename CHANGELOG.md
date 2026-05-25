# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

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
