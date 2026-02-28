from __future__ import annotations

import random
from typing import Any


def stratified_split(
    records: list[dict[str, Any]],
    train_frac: float = 0.7,
    validate_frac: float = 0.15,
    seed: int = 42,
    label_key: str = "taxid",
) -> tuple[list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]]]:
    """Split records into train/validation/test with label stratification.

    Splits a dataset into train/val/test partitions while preserving the label
    distribution for the provided `label_key` (default: `taxid`). The split is
    deterministic given a fixed seed.

    The test fraction is inferred as:
        test_frac = 1.0 - train_frac - val_frac

    Args:
        records: Dataset records to split. Each record must contain `label_key`.
        train_frac: Fraction of examples to place into the training split.
        validate_frac: Fraction of examples to place into the validation split.
        seed: Random seed for deterministic shuffling within each label group.
        label_key: Field name to stratify on (e.g., "taxid" or "label").

    Returns:
        (train_records, val_records, test_records) with no overlap.

    Raises:
        ValueError: If fractions are invalid or if records are missing `label_key`.
    """
    if not 0.0 < train_frac < 1.0:
        raise ValueError("Train_frac must be in (0,1)")
    if not 0.0 < validate_frac < 1.0:
        raise ValueError("validate_frac must be in (0,1)")
    if 1 - (train_frac + validate_frac) < 0.0:
        raise ValueError("Test_frac must be in (0,1)")

    groups: dict[Any, list[dict[str, Any]]] = {}
    for record in records:
        if label_key not in record:
            raise ValueError(f"Record Missing label_key='{label_key}'")
        groups.setdefault(record[label_key], []).append(record)

    rng = random.Random(seed)
    train: list[dict[str, Any]] = []
    validate: list[dict[str, Any]] = []
    test: list[dict[str, Any]] = []

    for _, items in groups.items():
        items = list(items)
        rng.shuffle(items)

        n = len(items)
        n_train = int(round(train_frac * n))
        n_validate = int(round(validate_frac * n))

        if n_train + n_validate > n:
            n_validate = max(0, n - n_train)

        train.extend(items[:n_train])
        validate.extend(items[n_train : n_train + n_validate])
        test.extend(items[n_train + n_validate :])
    rng.shuffle(train)
    rng.shuffle(validate)
    rng.shuffle(test)

    return train, validate, test
