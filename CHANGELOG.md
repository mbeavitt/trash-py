# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

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
