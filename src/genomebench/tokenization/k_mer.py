from itertools import product

from genomebench.data_processing.cleaning import clean_seq


class KMerTokenizer:
    """K-mer tokenizer for DNA sequences."""

    def __init__(self, k: int = 3) -> None:
        """Initialize a fixed k-mer vocabulary over A, C, G, T, N.

        Args:
            k: K-mer length.

        Raises:
            ValueError: If `k` is less than 1.
        """
        if k < 1:
            raise ValueError(f"`k` must be >= 1, got {k}")

        self._k = k
        bases = ("A", "C", "G", "T", "N")
        self._id_to_token = ["".join(chars) for chars in product(bases, repeat=k)]
        self._token_to_id = {token: idx for idx, token in enumerate(self._id_to_token)}

    def encode(self, seq: str) -> list[int]:
        """Encode a sequence into overlapping k-mer token IDs.

        The input is normalized with `clean_seq` before tokenization.

        Args:
            seq: Input nucleotide sequence (raw or cleaned).

        Returns:
            List of token IDs for windows of length `k` and stride 1.
            Returns an empty list when the normalized sequence is shorter than `k`.
        """
        seq = clean_seq(seq)
        if len(seq) < self._k:
            return []

        encoded: list[int] = []
        for i in range(len(seq) - self._k + 1):
            kmer = seq[i : i + self._k]
            encoded.append(self._token_to_id[kmer])
        return encoded

    def decode(self, ids: list[int]) -> str:
        """Decode k-mer token IDs into a nucleotide sequence.

        Args:
            ids: K-mer token IDs to decode.

        Returns:
            Reconstructed sequence string.

        Raises:
            ValueError: If any token ID is outside the vocabulary range.
            ValueError: If decoded k-mers are not overlap-consistent.
        """
        if not ids:
            return ""

        kmers: list[str] = []
        for idx in ids:
            if 0 <= idx < len(self._id_to_token):
                kmers.append(self._id_to_token[idx])
            else:
                raise ValueError(f"Outside vocabulary range: {idx}")

        if self._k == 1:
            return "".join(kmers)

        decoded = kmers[0]
        for pos, kmer in enumerate(kmers[1:], start=1):
            if decoded[-(self._k - 1) :] != kmer[: self._k - 1]:
                raise ValueError(
                    "Inconsistent k-mer overlap at position "
                    f"{pos}: {decoded[-(self._k - 1) :]} vs {kmer[: self._k - 1]}"
                )
            decoded += kmer[-1]
        return decoded

    @property
    def k(self) -> int:
        """Return k-mer length."""
        return self._k

    @property
    def vocab_size(self) -> int:
        """Return tokenizer vocabulary size."""
        return len(self._id_to_token)
