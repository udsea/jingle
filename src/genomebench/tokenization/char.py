from genomebench.data_processing.cleaning import clean_seq


class CharTokenizer:
    """Character-level tokenizer for DNA sequences."""

    def __init__(self) -> None:
        """Initializes a character-level tokenizer for DNA sequences.

        The vocabulary is fixed to five tokens: A, C, G, T, N.
        Any unknown character encountered during encoding is mapped to 'N'.
        """
        self._token_to_id = {"A": 0, "C": 1, "G": 2, "T": 3, "N": 4}
        self._id_to_token = ["A", "C", "G", "T", "N"]

    def encode(self, seq: str) -> list[int]:
        """Encodes a nucleotide sequence into token IDs.

        The input is normalized via `clean_seq` (uppercase, U->T, non-IUPAC->N).

        Args:
            seq: Input nucleotide sequence (raw or cleaned).

        Returns:
            A list of integer token IDs with the same length as the normalized
            sequence.
        """
        seq = clean_seq(seq)
        encoded: list[int] = []
        for char in seq:
            encoded.append(
                self._token_to_id.get(char, self._token_to_id["N"])
            )  # get gets key if in dict or esle default
        return encoded

    def decode(self, ids: list[int]) -> str:
        """Decodes token IDs back into a nucleotide sequence string.

        Args:
            ids: Token IDs to decode.

        Returns:
            A string sequence over the alphabet {A, C, G, T, N}.

        Raises:
            ValueError: If any ID is outside the valid vocabulary range.
        """
        decoded: list[str] = []
        for idx in ids:
            if 0 <= idx < len(self._id_to_token):
                decoded.append(self._id_to_token[idx])
            else:
                raise ValueError(f"Outside vocabulary range: {idx}")
        return "".join(decoded)

    @property
    def vocab_size(self) -> int:
        """Returns the size of the tokenizer vocabulary."""
        return len(self._id_to_token)
