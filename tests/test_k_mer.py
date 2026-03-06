"""Tests for the KMerTokenizer."""

import pytest

from genomebench.tokenization.k_mer import KMerTokenizer


def test_init_rejects_non_positive_k() -> None:
    """Tokenizer should reject k < 1."""
    with pytest.raises(ValueError):
        KMerTokenizer(k=0)


def test_vocab_size_for_k2_is_25() -> None:
    """Vocab size should be 5^k over A,C,G,T,N."""
    tok = KMerTokenizer(k=2)
    assert tok.vocab_size == 25


def test_encode_k2_basic_ids() -> None:
    """Encoding for k=2 should use overlapping windows."""
    tok = KMerTokenizer(k=2)
    assert tok.encode("ACGTN") == [1, 7, 13, 19]


def test_encode_normalizes_input() -> None:
    """Encoding should use cleaned sequence (uppercase, U->T, unknown->N)."""
    tok = KMerTokenizer(k=2)
    # "acuX" -> "ACTN", windows: AC, CT, TN
    assert tok.encode("acuX") == [1, 8, 19]


def test_encode_short_sequence_returns_empty() -> None:
    """Encoding should return empty list when len(seq) < k."""
    tok = KMerTokenizer(k=3)
    assert tok.encode("AC") == []


def test_decode_roundtrip_matches_cleaned_sequence() -> None:
    """decode(encode(seq)) should equal cleaned sequence when len >= k."""
    tok = KMerTokenizer(k=3)
    seq = "acguX"
    ids = tok.encode(seq)  # cleaned sequence = "ACGTN"
    assert tok.decode(ids) == "ACGTN"


def test_decode_empty_ids_returns_empty_string() -> None:
    """Decoding empty ids should return empty string."""
    tok = KMerTokenizer(k=4)
    assert tok.decode([]) == ""


def test_decode_raises_on_out_of_range_id() -> None:
    """Decoding invalid ids should raise ValueError."""
    tok = KMerTokenizer(k=2)
    with pytest.raises(ValueError):
        tok.decode([0, 9999])


def test_decode_raises_on_inconsistent_overlap() -> None:
    """Decoding should fail when adjacent k-mers do not overlap."""
    tok = KMerTokenizer(k=3)
    # AAA then CCC do not overlap on two bases.
    with pytest.raises(ValueError):
        tok.decode([0, 31])
