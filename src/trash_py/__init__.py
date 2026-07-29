"""trash-py — Python port of TRASH.

A tandem-repeat array identifier for assembled genomes. See `cli.py`
for the command-line entry point, `pipeline.run_pipeline` to drive
the pipeline programmatically, or `reads.annotate` / `reads.annotate_sequences`
for in-memory annotation of individual sequences (e.g. long reads).
"""

__version__ = "2.7.1"

from .reads import (  # noqa: E402
    RepeatArray,
    RepeatUnit,
    SequenceAnnotation,
    annotate,
    annotate_sequences,
)
from .templates import ARABIDOPSIS_CEN178  # noqa: E402

__all__ = [
    "annotate",
    "annotate_sequences",
    "RepeatArray",
    "RepeatUnit",
    "SequenceAnnotation",
    "ARABIDOPSIS_CEN178",
]
