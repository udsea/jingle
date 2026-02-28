from pathlib import Path
from typing import Any

from genomebench.data_processing.cleaning import ambig_fraction, clean_seq
from genomebench.io.fasta import iter_fasta


def build_dataset_from_fasta(
    fasta_path: Path,
    taxid: int,
    min_length: int = 5000,
    max_ambig_frac: float = 0.05,
) -> list[dict[str, Any]]:
    """Build cleaned and filtered sequence records from a FASTA file.

    Reads sequences from `fasta_path`, normalizes them
    (uppercase, U->T, ambiguity->N),
    computes ambiguity fraction, and filters records by length and ambiguity.

    Args:
        fasta_path: Path to an input FASTA file containing viral genome sequences.
        taxid: NCBI taxonomy identifier to attach to all returned records.
        min_length: Minimum sequence length (in bases) to keep a record.
        max_ambig_frac: Maximum allowed fraction of ambiguous bases.

    Returns:
        A list of JSON-serializable dict records. Each record includes at least:
        `sequence_id`, `taxid`, `length`, `ambig_fraction`, and `sequence`.
    """
    results: list[dict[str, Any]] = []

    for header, seq in iter_fasta(fasta_path):
        cleaned = clean_seq(seq)
        seq_ambig_frac = ambig_fraction(cleaned)
        if len(cleaned) < min_length:
            continue
        if seq_ambig_frac > max_ambig_frac:
            continue
        results.append(
            {
                "sequence_id": header,
                "taxid": taxid,
                "length": len(cleaned),
                "ambig_fraction": seq_ambig_frac,
                "sequence": cleaned,
            }
        )
    return results
