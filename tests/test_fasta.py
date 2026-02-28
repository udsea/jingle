from pathlib import Path

from genomebench.io.fasta import iter_fasta


def test_iter_fasta_multiline(tmp_path: Path):
    """Test fucntion for iter_fasta."""
    p = tmp_path / "tiny.fasta"
    p.write_text(">seq1\nACG\nTT\n>seq2 some desc\nGGG\n")

    rows = list(iter_fasta(p))
    assert rows == [
        ("seq1", "ACGTT"),
        ("seq2 some desc", "GGG"),
    ]
