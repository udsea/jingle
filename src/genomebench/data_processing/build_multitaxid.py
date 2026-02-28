import random
from pathlib import Path
from typing import Any

from genomebench.data_processing.build_dataset import build_dataset_from_fasta


def build_multitaxid_dataset(
    fasta_map: dict[int, Path],
    min_length: int = 5000,
    max_ambig_frac: float = 0.05,
    seed: int = 42,
) -> list[dict[str, Any]]:
    """Build a combined dataset from multiple taxids and FASTA files.

    Args:
        fasta_map: Mapping from NCBI taxid to input FASTA path.
                    Each FASTA is assumed to contain genome sequences
                    for that taxon.
        min_length: Minimum sequence length (in bases) required to keep a record.
        max_ambig_frac: Maximum allowed fraction of ambiguous bases ('N') required to
            keep a record.
        seed: Random seed used to deterministically shuffle the combined dataset.

    Returns:
        A list of JSON-serializable dict records. Each record includes at least:
        `sequence_id`, `taxid`, `length`, `ambig_frac`, and `sequence`.
    """
    all_records: list[dict[str, Any]] = []
    rng = random.Random(seed)

    for taxid, fasta_path in fasta_map.items():
        records = build_dataset_from_fasta(
            fasta_path=fasta_path,
            taxid=taxid,
            min_length=min_length,
            max_ambig_frac=max_ambig_frac,
        )
        for r in records:
            r["label_name"] = taxid
        all_records.extend(records)

    rng.shuffle(all_records)
    return all_records
