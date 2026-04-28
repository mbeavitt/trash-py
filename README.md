# trash-py

Tandem-repeat array identifier — a Python port of the TRASH program written in R.

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

## Installation
Please install [nhmmer](http://hmmer.org/) and [Clustal Omega](https://bioconda.github.io/recipes/clustalo/README.html)
and ensure they are available on the PATH. Additionally, please ensure you have
a suitable C compiler installed (gcc, clang).

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
as possible, to present a drag-and-drop replacement.

The exceptions are that -q is now available to silence logs, and
currently __-p for multiprocessing is not supported__. This may be implemented in
future - currently however, it's recommended to split input .fasta files and
parallelise by chromosome/sequence using an external tool like GNU Parallel, and
merge later.

## Benchmarks


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
