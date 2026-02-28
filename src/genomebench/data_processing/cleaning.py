VALID = {"A", "C", "T", "G", "N"}
IUPAC_BASES = {"R", "Y", "S", "W", "K", "M", "B", "D", "H", "V"}


def clean_seq(seq: str) -> str:
    """Normalize and clean a nucleotide sequence.

    Converts to uppercase, maps RNA base 'U' to 'T',
    replaces IUPAC ambiguous bases and invalid characters with 'N'.

    Args:
        seq: Raw nucleotide sequence string.

    Returns:
        A cleaned sequence containing only A, C, G, T, N.
    """
    seq = seq.upper()
    seq = seq.replace("U", "T")
    out: list[str] = []
    for chars in seq:
        if chars in VALID:
            out.append(chars)
        elif chars in IUPAC_BASES:
            out.append("N")
        else:
            out.append("N")

    return "".join(out)


def ambig_fraction(seq: str) -> float:
    """Compute fraction of ambiguous bases in a sequence.

    Ambiguity is defined as the proportion of 'N' bases.

    Args:
        seq: Cleaned nucleotide sequence.

    Returns:
        Fraction of ambiguous bases in [0.0, 1.0].
        Returns 0.0 for empty sequences.
    """
    ambig_frac = 0.0
    if len(seq) != 0:
        ambig_frac = seq.count("N") / len(seq)
    else:
        ambig_frac = 0.0

    return ambig_frac
