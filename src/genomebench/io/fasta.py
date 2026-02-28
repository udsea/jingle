from collections.abc import Iterator
from pathlib import Path


def iter_fasta(path: Path) -> Iterator[tuple[str, str]]:
    """Stream a FASTA file and yield (header, sequence).

    The header is the text after '>' (without the '>').
    Sequence lines are concatenated across multiline blocks.
    Blank lines are ignored.

    Args:
        path: Path to a FASTA file.

    Yields:
        Tuples of (header, sequence).
    """
    header: str | None = None
    seq_chunks: list[str] = []

    with path.open("r", encoding="utf-8", errors="ignore") as file:
        for raw in file:
            line = raw.strip()
            if not line:
                continue

            if line.startswith(">"):
                if header is not None:
                    yield (header, "".join(seq_chunks))
                header = line[1:].strip()
                seq_chunks = []
            else:
                if header is None:
                    raise ValueError(
                        f"FASTA file {path} has sequence data before header"
                    )
                seq_chunks.append(line)

        if header is not None:
            yield (header, "".join(seq_chunks))
