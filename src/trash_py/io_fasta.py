"""Port of `src/in_out/read_fasta_and_list.R`."""
from __future__ import annotations

from pathlib import Path


def read_fasta_and_list(path: Path) -> list[tuple[str, str]]:
    """Return [(name, sequence_lowercase), ...].

    Matches `ape::read.dna(..., as.character=TRUE, as.matrix=FALSE)`:
    sequences are lowercased; headers are split on whitespace and only
    the first token is kept. Dashes and ambiguity codes are preserved
    as-is.
    """
    records: list[tuple[str, str]] = []
    name: str | None = None
    chunks: list[str] = []
    with Path(path).open() as f:
        for line in f:
            line = line.rstrip("\r\n")
            if line.startswith(">"):
                if name is not None:
                    records.append((name, "".join(chunks).lower()))
                name = line[1:].split()[0] if line[1:].strip() else ""
                chunks = []
            else:
                chunks.append(line)
    if name is not None:
        records.append((name, "".join(chunks).lower()))
    return records
