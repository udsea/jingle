"""Tests for the CharTokenizer."""

import pytest

from genomebench.tokenization.char import CharTokenizer


def test_encode_basic_acgtn() -> None:
    """Encoding ACGTN should map to fixed ids."""
    tok = CharTokenizer()
    assert tok.encode("ACGTN") == [0, 1, 2, 3, 4]


def test_encode_normalizes_u_to_t() -> None:
    """Encoding should normalize U->T (RNA to DNA)."""
    tok = CharTokenizer()
    # U should become T during cleaning.
    assert tok.encode("ACGU") == [0, 1, 2, 3]


def test_encode_maps_unknown_to_n() -> None:
    """Unknown characters should map to N."""
    tok = CharTokenizer()
    # X should become N after cleaning.
    assert tok.encode("AXGT") == [0, 4, 2, 3]


def test_decode_basic_roundtrip() -> None:
    """decode(encode(seq)) should equal the cleaned seq."""
    tok = CharTokenizer()
    ids = tok.encode("acguX")  # lower + U + unknown
    # clean_seq would yield "ACGTN"
    assert tok.decode(ids) == "ACGTN"


def test_decode_raises_on_out_of_range_id() -> None:
    """Decoding invalid ids should raise ValueError."""
    tok = CharTokenizer()
    with pytest.raises(ValueError):
        tok.decode([0, 1, 999])


def test_vocab_size_is_five() -> None:
    """Vocab size should be 5 for A,C,G,T,N."""
    tok = CharTokenizer()
    assert tok.vocab_size == 5
