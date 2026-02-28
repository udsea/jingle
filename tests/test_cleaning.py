from genomebench.data_processing.cleaning import ambig_fraction, clean_seq


def test_clean_seq_uppercase_and_u_to_t() -> None:
    """Uppercases sequence and converts RNA base U to DNA base T."""
    assert clean_seq("acgu") == "ACGT"


def test_clean_seq_maps_iupac_to_n() -> None:
    """Maps IUPAC ambiguous bases to 'N'."""
    assert clean_seq("RYSWKMBDHV") == "N" * 10


def test_clean_seq_maps_invalid_to_n() -> None:
    """Maps invalid/non-IUPAC characters to 'N'."""
    assert clean_seq("ACGT-XYZ!") == "ACGTNNNNN"


def test_clean_seq_keeps_valid_bases_and_n() -> None:
    """Leaves valid bases A/C/G/T/N unchanged (aside from normalization)."""
    assert clean_seq("ACGTNacgtn") == "ACGTNACGTN"


def test_ambig_fraction_empty_is_zero() -> None:
    """Returns 0.0 ambiguity fraction for empty sequences."""
    assert ambig_fraction("") == 0.0


def test_ambig_fraction_no_ambiguity() -> None:
    """Returns 0.0 when there are no 'N' bases."""
    assert ambig_fraction("ACGTACGT") == 0.0


def test_ambig_fraction_all_ambiguity() -> None:
    """Returns 1.0 when all bases are 'N'."""
    assert ambig_fraction("NNNN") == 1.0


def test_ambig_fraction_mixed() -> None:
    """Computes correct fraction when sequence has a mix of bases and 'N'."""
    assert ambig_fraction("ACNNG") == 2 / 5
