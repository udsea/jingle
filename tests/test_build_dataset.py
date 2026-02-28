from pathlib import Path

from genomebench.data_processing.build_dataset import build_dataset_from_fasta


def test_build_dataset_keeps_valid_record(tmp_path: Path) -> None:
    """Keeps sequences that pass length and ambiguity filters."""
    p = tmp_path / "tiny.fasta"
    # seq length = 6, ambig = 0/6
    p.write_text(">seq1\nACGTAC\n", encoding="utf-8")

    rows = build_dataset_from_fasta(p, taxid=123, min_length=5, max_ambig_frac=0.05)

    assert len(rows) == 1
    r = rows[0]
    assert r["sequence_id"] == "seq1"
    assert r["taxid"] == 123
    assert r["length"] == 6
    assert r["ambig_fraction"] == 0.0
    assert r["sequence"] == "ACGTAC"


def test_build_dataset_filters_short_sequences(tmp_path: Path) -> None:
    """Drops sequences shorter than min_length."""
    p = tmp_path / "tiny.fasta"
    p.write_text(">seq1\nACGT\n", encoding="utf-8")  # length 4

    rows = build_dataset_from_fasta(p, taxid=1, min_length=5, max_ambig_frac=1.0)

    assert rows == []


def test_build_dataset_filters_high_ambiguity(tmp_path: Path) -> None:
    """Drops sequences whose N fraction exceeds max_ambig_frac."""
    p = tmp_path / "tiny.fasta"
    # cleaned stays same; ambig_frac = 3/6 = 0.5
    p.write_text(">seq1\nNNNAAA\n", encoding="utf-8")

    rows = build_dataset_from_fasta(p, taxid=1, min_length=1, max_ambig_frac=0.49)

    assert rows == []


def test_build_dataset_u_to_t_and_iupac_mapping(tmp_path: Path) -> None:
    """Cleans sequences before filtering (U->T, IUPAC->N)."""
    p = tmp_path / "tiny.fasta"
    # U->T, R->N => "ACTN", length 4, ambig 1/4 = 0.25
    p.write_text(">seq1\nACUR\n", encoding="utf-8")

    rows = build_dataset_from_fasta(p, taxid=9, min_length=1, max_ambig_frac=0.25)

    assert len(rows) == 1
    r = rows[0]
    assert r["sequence"] == "ACTN"
    assert r["ambig_fraction"] == 0.25
