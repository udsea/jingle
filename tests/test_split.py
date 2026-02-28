from __future__ import annotations

from typing import Any

import pytest

from genomebench.data_processing.split import stratified_split


def _ids(xs: list[dict[str, Any]]) -> set[str]:
    return {r["sequence_id"] for r in xs}


def _count_label(xs: list[dict[str, Any]], label: int) -> int:
    return sum(1 for r in xs if r["taxid"] == label)


def test_stratified_split_no_overlap_and_expected_counts() -> None:
    """Splits are disjoint and preserve per-label counts (stable rounding case)."""
    records = [{"taxid": 1, "sequence_id": f"a{i}"} for i in range(10)] + [
        {"taxid": 2, "sequence_id": f"b{i}"} for i in range(10)
    ]

    train, validate, test = stratified_split(
        records, train_frac=0.7, validate_frac=0.2, seed=123
    )

    assert _ids(train).isdisjoint(_ids(validate))
    assert _ids(train).isdisjoint(_ids(test))
    assert _ids(validate).isdisjoint(_ids(test))

    assert len(train) + len(validate) + len(test) == len(records)

    assert _count_label(train, 1) == 7
    assert _count_label(validate, 1) == 2
    assert _count_label(test, 1) == 1

    assert _count_label(train, 2) == 7
    assert _count_label(validate, 2) == 2
    assert _count_label(test, 2) == 1


def test_stratified_split_deterministic_given_seed() -> None:
    """Same seed -> identical split ordering; different seed -> usually different."""
    records = [{"taxid": 1, "sequence_id": f"x{i}"} for i in range(25)] + [
        {"taxid": 2, "sequence_id": f"y{i}"} for i in range(25)
    ]

    out1 = stratified_split(records, seed=0)
    out2 = stratified_split(records, seed=0)
    assert out1 == out2

    out3 = stratified_split(records, seed=1)
    assert out1 != out3


def test_stratified_split_raises_on_missing_label_key() -> None:
    """Missing label_key in any record should raise ValueError."""
    records = [{"sequence_id": "oops"}]
    with pytest.raises(ValueError):
        stratified_split(records)


def test_stratified_split_raises_on_invalid_fracs_sum_gt_one() -> None:
    """train_frac + validate_frac > 1 should raise ValueError."""
    records = [{"taxid": 1, "sequence_id": "a"}]
    with pytest.raises(ValueError):
        stratified_split(records, train_frac=0.9, validate_frac=0.2)


def test_stratified_split_handles_small_groups_without_crashing() -> None:
    """Tiny label groups should still split without overlap or crashes."""
    records = [
        {"taxid": 1, "sequence_id": "a0"},
        {"taxid": 2, "sequence_id": "b0"},
        {"taxid": 2, "sequence_id": "b1"},
    ]

    train, validate, test = stratified_split(records, seed=7)

    assert len(train) + len(validate) + len(test) == len(records)
    assert _ids(train).isdisjoint(_ids(validate))
    assert _ids(train).isdisjoint(_ids(test))
    assert _ids(validate).isdisjoint(_ids(test))
