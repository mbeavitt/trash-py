# trash-py

Tandem-repeat array identifier — a Python port of the TRASH program.

## Origin and acknowledgements

trash-py is a Python re-implementation of the TRASH program written by
Piotr Włodzimierz (<pwlodzimierz@ibb.waw.pl>, Institute of Biochemistry
and Biophysics, Polish Academy of Sciences). The upstream repository lives
at <https://github.com/vlothec/TRASH_2>. All algorithmic credit for the
underlying approach belongs to the original author; this port is an
independent rewrite (no code is byte-identical with the upstream) that
retains the structure and logic of the upstream pipeline.

The upstream MIT license is reproduced in
[`LICENSES/TRASH-UPSTREAM-MIT.txt`](LICENSES/TRASH-UPSTREAM-MIT.txt) to
satisfy its notice-preservation clause.

## What's different?
trash-py aims to build on the substantial work done in developing the original
TRASH repeat annotation pipeline by adopting a more flexible (and performant)
Python/C runner/libs framework. A library of reusable functions is exposed which
can be incorporated into diverse repeat-annotation related workflows beyond
simply running the TRASH pipeline, and the hotter functions have been ported to
C to maximise performance.

Currently, on smaller less repetitive genomes (e.g. Arabidopsis, Human genome)
the bottleneck is nhmmer rather than the the TRASH pipeline itself, which takes
a fraction of the original time to complete.

The stderr logs have also been changed to reflect the new internals and to give
trash-py a little of its own personality.

## Installation
### conda

The best way to install the tool is using conda/mamba/micromamba:

```
conda install -c bioconda trash-py
```
### source
If instead you'd like to build/install from source, please install [nhmmer](http://hmmer.org/) and 
[Clustal Omega](https://bioconda.github.io/recipes/clustalo/README.html)
and ensure they are available on the PATH. HOR detection additionally needs
[MAFFT](https://mafft.cbrc.jp/alignment/software/) on the PATH (and `matplotlib`
for its plots — `pip install trash-py[plot]`). Additionally, please ensure you have
a suitable C compiler installed (gcc, clang) for the C extensions. 
There is no need to separately compile these, python should recognise the instructions for compilation in setup.py.

trash-py can then be installed by:

```
git clone https://github.com/mbeavitt/trash-py
cd trash-py
pip install .
```

## Usage

```
trash-py -f input.fasta -o output_dir
```

Currently the CLI aims to mirror the one in the original TRASH tool as closely
as possible, to present a drag-and-drop replacement. trash-py adds two options
over upstream: `-q` to silence logs, and `-p` to run the array-identification
and repeat-mapping stages across multiple worker processes.

### HOR detection

Higher-order-repeat (HOR) detection — the port of the upstream `HORT.R` module —
can be run two ways.

1. Inline HOR detection, as part of the main pipeline — add `--hor-chr-list` to a
normal run and HORs are detected on the repeat table the pipeline just produced, using the most abundant class of repeats (e.g. 178_1 in Arabidopsis):

```
trash-py -f genome.fasta -o out --hor-chr-list Chr1,Chr2,Chr3
```

2. Standalone, to re-run HOR finding on an existing
`<fasta>_repeats_with_seq.csv` without redoing the whole pipeline — handy for targeting a different class or set of sequences. `-c/--class` picks the class explicitly — any class, e.g. `178_1`, `178_2`, or, as below, the 113 bp `113_15` satellite:

```
# a non-primary class across several sequences
trash-py hor out/genome.fasta_repeats_with_seq.csv -c 113_15 --chr-list Chr1,Chr2 -o out

# a single sequence
trash-py hor out/genome.fasta_repeats_with_seq.csv -c 113_15 --ChrA Chr1 -o out

# cross-region comparison (region A vs region B)
trash-py hor out/genome.fasta_repeats_with_seq.csv -c 113_15 --ChrA Chr1 --ChrB Chr2 -o out
```

You don't have to know the satellite class up front: **the class is optional and
defaults to the most abundant class on the target sequence(s)** — usually the
centromeric satellite — with the choice logged. Override with `--hor-class`
(pipeline) / `-c/--class` (subcommand), or run it again for another class.

For each sequence this writes `HORs_<class>_<seq>.csv` (the HOR table),
`repeats_with_hors_<class>_<seq>.csv` (per-repeat annotation), and a
`HORs_lines_<class>_<seq>.png` dot-plot. Quirks preserved from the original tool are
documented in [`docs/HOR_source_bugs.md`](docs/HOR_source_bugs.md).

## Benchmarks

![trash-py runtime and parallel speedup vs. process count](docs/images/runtime_plot.png)

Benchmarked on the *A. thaliana* GCA_028009825.2 assembly (144 MB) on a
32-core host. Increasing `-p` cuts wall-clock runtime from ~230 s at `-p 1`
to ~59 s at `-p 16`. Speedup is sub-linear and flattens beyond 16 processes,
so `-p 16` is a reasonable default on this class of host.

## How to cite

If you use `trash-py` in academic work, please cite the original TRASH
publication:

Wlodzimierz, P., Hong, M., & Henderson, I. R. (2023). TRASH: tandem
repeat annotation and structural hierarchy. *Bioinformatics*, 39(5),
btad308.

BibTeX:

```bibtex
@article{wlodzimierz2023trash,
  title={TRASH: tandem repeat annotation and structural hierarchy},
  author={Wlodzimierz, Piotr and Hong, Michael and Henderson, Ian R},
  journal={Bioinformatics},
  volume={39},
  number={5},
  pages={btad308},
  year={2023},
  publisher={Oxford University Press}
}
```

## License

`trash-py` is released under the MIT License — see [`LICENSE`](LICENSE).

## Command-line reference

```
$ trash-py --help
usage: trash-py [-h] [-V] -f FASTA -o OUTPUT [-m MAX_REP_SIZE]
                [-i MIN_REP_SIZE] [-t TEMPLATES] [-q] [-p PROCESSES]

TRASH — tandem-repeat array identifier (Python)

options:
  -h, --help            show this help message and exit
  -V, --version         show program's version number and exit
  -f, --fasta FASTA     input fasta
  -o, --output OUTPUT   output directory
  -m, --max-rep-size MAX_REP_SIZE
  -i, --min-rep-size MIN_REP_SIZE
  -t, --templates TEMPLATES
                        optional template fasta — assigns class names from
                        headers
  -q, --quiet           suppress progress output
  -p, --processes PROCESSES
                        parallel worker processes for the array-identification
                        and repeat-mapping stages (default 1 = serial)
```
