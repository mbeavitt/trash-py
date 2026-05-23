"""Built-in sequence templates for TRASH classification.

Templates are ``(name, sequence)`` tuples that can be passed to
:func:`trash_py.annotate` (or :func:`trash_py.annotate_sequences`) via the
``templates`` parameter.  Supplying a template anchors the canonical rotation
used by ``shift_and_compare`` so that all arrays are oriented consistently
relative to the reference sequence, making cross-read comparisons valid.

Example::

    from trash_py import annotate
    from trash_py.templates import ARABIDOPSIS_CEN178

    result = annotate(read_seq, name="read1", templates=[ARABIDOPSIS_CEN178])
    for unit in result.repeats:
        print(unit.start, unit.end, unit.class_)
"""

# Arabidopsis thaliana CEN178 centromeric satellite consensus sequence (177 bp).
# Source: CEN178 canonical monomer — Arabidopsis thaliana centromeric repeat.
# Used as template to ensure consistent canonical rotation when annotating
# centromeric satellite arrays across multiple reads.
ARABIDOPSIS_CEN178: tuple[str, str] = (
    "CEN178",
    "AGTATAAGAACTTAAACCGCAACCGATCTTAAAAGCCTAAGTAGTGTTTCCTTGTTAGAA"
    "GACACAAAGCCAAAGACTCATATGGACTTTGGCTACACCATGAAAGCTTTGAGAAGCAAG"
    "AAGAAGGTTGGTTAGTGTTTTGGAGTCGAATATGACTTGATGTCATGTGTATGATTG",
)
