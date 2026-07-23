"""Port of `src/in_out/read_fasta_and_list.R`."""
from __future__ import annotations

import gzip
import io
from pathlib import Path

GZIP_MAGIC = b"\x1f\x8b"


def _open_text(path: Path) -> io.TextIOBase:
    """Open `path` for text reading, decompressing gzip transparently.

    Detection is by magic bytes rather than by suffix so that pipes work
    (`-f -` resolves to /dev/stdin, which has no informative name). `peek`
    inspects the buffer without consuming it, which is what makes this
    safe on a non-seekable stream.
    """
    raw = open(path, "rb")
    if raw.peek(len(GZIP_MAGIC))[:len(GZIP_MAGIC)] == GZIP_MAGIC:
        return io.TextIOWrapper(gzip.open(raw, "rb"))
    return io.TextIOWrapper(raw)


def read_fasta_and_list(path: Path) -> list[tuple[str, str]]:
    """Return [(name, sequence_lowercase), ...].

    Matches `ape::read.dna(..., as.character=TRUE, as.matrix=FALSE)`:
    sequences are lowercased; headers are split on whitespace and only
    the first token is kept. Dashes and ambiguity codes are preserved
    as-is. Gzipped inputs are decompressed transparently, detected by
    magic bytes rather than suffix so that pipes (`-f -`) work too.
    """
    path = Path(path)
    records: list[tuple[str, str]] = []
    name: str | None = None
    chunks: list[str] = []
    with _open_text(path) as f:
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
