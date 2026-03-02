import argparse
import logging
import os
import struct
import sys

import torch
import torch.nn as nn
import torch.nn.functional as F
from Bio import SeqIO

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(message)s",
)
log = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Amino acid alphabet (26 symbols)
# ---------------------------------------------------------------------------
_AMINO_STANDARD = "ACDEFGHIKLMNPQRSTVWY"
_AMINO_ALPHABET_SIZE = 26

_AMINO_TO_INT = {}
for _i, _aa in enumerate(_AMINO_STANDARD):
    _AMINO_TO_INT[_aa] = _i

_AMINO_TO_INT["B"] = 20  # Asx (D or N)
_AMINO_TO_INT["Z"] = 21  # Glx (E or Q)
_AMINO_TO_INT["J"] = 22  # Xle (I or L)
_AMINO_TO_INT["X"] = 23  # Unknown
_AMINO_TO_INT["*"] = 24  # Stop
_AMINO_TO_INT["-"] = 25  # Gap
_AMINO_TO_INT["."] = 25  # Gap (alternate)
_AMINO_TO_INT["U"] = _AMINO_TO_INT["C"]  # Selenocysteine -> Cysteine
_AMINO_TO_INT["O"] = _AMINO_TO_INT["K"]  # Pyrrolysine -> Lysine

_AMINO_LOOKUP = [23] * 256
for _aa, _idx in _AMINO_TO_INT.items():
    _AMINO_LOOKUP[ord(_aa)] = _idx
    _AMINO_LOOKUP[ord(_aa.lower())] = _idx
_AMINO_LOOKUP_TENSOR = torch.tensor(_AMINO_LOOKUP, dtype=torch.long)

# Amino ambiguation tables
_AMINO_AMBIG_SPECIFIC = [23] * 26  # default: X
_AMINO_AMBIG_SPECIFIC[2] = 20  # D -> B
_AMINO_AMBIG_SPECIFIC[11] = 20  # N -> B
_AMINO_AMBIG_SPECIFIC[3] = 21  # E -> Z
_AMINO_AMBIG_SPECIFIC[13] = 21  # Q -> Z
_AMINO_AMBIG_SPECIFIC[7] = 22  # I -> J
_AMINO_AMBIG_SPECIFIC[9] = 22  # L -> J
_AMINO_AMBIG_SPECIFIC = torch.tensor(_AMINO_AMBIG_SPECIFIC, dtype=torch.long)

_AMINO_HAS_SPECIFIC = torch.zeros(26, dtype=torch.bool)
for _idx in (2, 3, 7, 9, 11, 13):
    _AMINO_HAS_SPECIFIC[_idx] = True

# Full amino alphabet string for CSV header (indices 0-25)
_AMINO_ALPHABET_STR = _AMINO_STANDARD + "BZJX*-"

# ---------------------------------------------------------------------------
# DNA alphabet (15 symbols: 4 standard + 11 IUPAC ambiguous)
# ---------------------------------------------------------------------------
_DNA_STANDARD = "ACGT"
_DNA_ALPHABET_SIZE = 15

_DNA_TO_INT = {"A": 0, "C": 1, "G": 2, "T": 3, "U": 3}
_DNA_TO_INT["R"] = 4  # A or G
_DNA_TO_INT["Y"] = 5  # C or T
_DNA_TO_INT["S"] = 6  # G or C
_DNA_TO_INT["W"] = 7  # A or T
_DNA_TO_INT["K"] = 8  # G or T
_DNA_TO_INT["M"] = 9  # A or C
_DNA_TO_INT["B"] = 10  # C or G or T
_DNA_TO_INT["D"] = 11  # A or G or T
_DNA_TO_INT["H"] = 12  # A or C or T
_DNA_TO_INT["V"] = 13  # A or C or G
_DNA_TO_INT["N"] = 14  # any

_DNA_LOOKUP = [14] * 256  # default: N
for _ch, _idx in _DNA_TO_INT.items():
    _DNA_LOOKUP[ord(_ch)] = _idx
    _DNA_LOOKUP[ord(_ch.lower())] = _idx
_DNA_LOOKUP_TENSOR = torch.tensor(_DNA_LOOKUP, dtype=torch.long)

# DNA ambiguation: for each standard base, list all IUPAC codes that cover it.
# A(0): R(4),W(7),M(9),D(11),H(12),V(13),N(14)
# C(1): Y(5),S(6),M(9),B(10),H(12),V(13),N(14)
# G(2): R(4),S(6),K(8),B(10),D(11),V(13),N(14)
# T(3): Y(5),W(7),K(8),B(10),D(11),H(12),N(14)
_DNA_AMBIG_OPTIONS = {
    0: [4, 7, 9, 11, 12, 13, 14],
    1: [5, 6, 9, 10, 12, 13, 14],
    2: [4, 6, 8, 10, 11, 13, 14],
    3: [5, 7, 8, 10, 11, 12, 14],
}
# Pad to uniform length for tensor ops (all have 7 options)
_DNA_AMBIG_TABLE = torch.tensor(
    [_DNA_AMBIG_OPTIONS[i] for i in range(4)], dtype=torch.long
)  # [4, 7]

_DNA_ALPHABET_STR = "ACGTRYSW KMBDHVN"
# Clean version: indices 0-14 map to these characters
_DNA_INDEX_TO_CHAR = list("ACGTRYSWKMBDHVN")

# Full DNA alphabet string for CSV header
_DNA_ALPHABET_STR = "ACGTRYSWKMBDHVN"

# ---------------------------------------------------------------------------
# Ambiguation cache
# ---------------------------------------------------------------------------
_amino_ambig_cache = {}
_dna_ambig_cache = {}


def _get_amino_ambig_tables(device):
    if device not in _amino_ambig_cache:
        _amino_ambig_cache[device] = (
            _AMINO_AMBIG_SPECIFIC.to(device),
            _AMINO_HAS_SPECIFIC.to(device),
        )
    return _amino_ambig_cache[device]


def _get_dna_ambig_table(device):
    if device not in _dna_ambig_cache:
        _dna_ambig_cache[device] = _DNA_AMBIG_TABLE.to(device)
    return _dna_ambig_cache[device]


# ---------------------------------------------------------------------------
# Ambiguation functions
# ---------------------------------------------------------------------------


def ambiguate_amino_batch(batch, rate):
    """Randomly replace amino acids with appropriate ambiguous codes."""
    dev = batch.device
    ambig_specific, has_specific = _get_amino_ambig_tables(dev)

    mask = torch.rand(batch.shape, device=dev) < rate
    mask = mask & (batch < 20)

    result = batch.clone()
    masked_vals = batch[mask]

    replacement = ambig_specific[masked_vals]
    use_x = has_specific[masked_vals] & (
        torch.rand(masked_vals.shape, device=dev) < 0.5
    )
    replacement[use_x] = 23

    result[mask] = replacement
    return result


def ambiguate_dna_batch(batch, rate):
    """Randomly replace DNA bases with IUPAC ambiguous codes (uniform choice)."""
    dev = batch.device
    ambig_table = _get_dna_ambig_table(dev)  # [4, 7]

    mask = torch.rand(batch.shape, device=dev) < rate
    mask = mask & (batch < 4)

    result = batch.clone()
    masked_vals = batch[mask]

    # For each masked position, pick a random column from the ambig table
    options = ambig_table[masked_vals]  # [num_masked, 7]
    pick = torch.randint(0, 7, (masked_vals.shape[0],), device=dev)
    replacement = options[torch.arange(masked_vals.shape[0], device=dev), pick]

    result[mask] = replacement
    return result


# ---------------------------------------------------------------------------
# Alphabet config helper
# ---------------------------------------------------------------------------


def get_alphabet_config(alphabet):
    """Return (lookup_tensor, alphabet_size, alphabet_str, ambiguate_fn, index_to_char)."""
    if alphabet == "amino":
        index_to_char = list(_AMINO_STANDARD + "BZJX*-")
        return (
            _AMINO_LOOKUP_TENSOR,
            _AMINO_ALPHABET_SIZE,
            _AMINO_ALPHABET_STR,
            ambiguate_amino_batch,
            index_to_char,
        )
    elif alphabet == "dna":
        return (
            _DNA_LOOKUP_TENSOR,
            _DNA_ALPHABET_SIZE,
            _DNA_ALPHABET_STR,
            ambiguate_dna_batch,
            _DNA_INDEX_TO_CHAR,
        )
    else:
        raise ValueError("Unknown alphabet: %s" % alphabet)


# ---------------------------------------------------------------------------
# Kmer extraction
# ---------------------------------------------------------------------------


def _scan_sequence_lengths(file_path, min_seq_len):
    """First pass: collect (record_index, seq_length) without keeping sequence data."""
    entries = []
    for i, record in enumerate(SeqIO.parse(file_path, "fasta")):
        n = len(record.seq)
        if n >= min_seq_len:
            entries.append((i, n))
    return entries


def encode_sequences(file_path, kmer_size, lookup_tensor, max_centers=0):
    """Encode FASTA sequences into a concatenated tensor with valid center positions.

    Returns (all_encoded, valid_centers) where valid_centers are global indices
    into all_encoded at which a center kmer can be placed such that left/right
    context kmers at up to 3*kmer_size distance stay within the same sequence.

    If max_centers > 0, randomly sample sequences from the file until enough
    valid centers have been collected, keeping memory usage proportional to the
    training budget while representing the full diversity of the input.
    """
    max_dist = 3 * kmer_size
    min_seq_len = max_dist + kmer_size + max_dist

    # -- Lightweight scan to get sequence lengths --
    log.info("Scanning %s for sequence lengths...", file_path)
    entries = _scan_sequence_lengths(file_path, min_seq_len)
    if not entries:
        return torch.empty(0, dtype=torch.long), torch.empty(0, dtype=torch.long)

    total_available = sum(
        n - 2 * max_dist - kmer_size + 1 for _, n in entries
    )
    log.info(
        "Found %d eligible sequences (%d available centers).",
        len(entries),
        total_available,
    )

    # -- Decide which sequences to load --
    if max_centers > 0 and total_available > max_centers:
        num_entries = len(entries)
        perm = torch.randperm(num_entries)
        rec_indices = torch.tensor([e[0] for e in entries], dtype=torch.long)
        seq_lengths = torch.tensor([e[1] for e in entries], dtype=torch.long)
        centers_per_seq = seq_lengths - 2 * max_dist - kmer_size + 1
        shuffled_centers = centers_per_seq[perm]
        cumsum = shuffled_centers.cumsum(0)
        num_needed = int((cumsum >= max_centers).long().argmax().item()) + 1
        selected_indices = set(rec_indices[perm[:num_needed]].tolist())
        collected = int(cumsum[num_needed - 1].item())
        log.info(
            "Sampled %d / %d sequences (%d centers) to meet budget of %d.",
            num_needed,
            num_entries,
            collected,
            max_centers,
        )
    else:
        selected_indices = None  # load all

    # -- Second pass: load only selected sequences --
    encoded_parts = []
    seq_ranges = []
    offset = 0

    for i, record in enumerate(SeqIO.parse(file_path, "fasta")):
        if selected_indices is not None and i not in selected_indices:
            continue
        seq = str(record.seq).upper()
        n = len(seq)
        if n < min_seq_len:
            continue
        ords = torch.tensor([ord(c) for c in seq], dtype=torch.long)
        encoded_parts.append(lookup_tensor[ords])
        seq_ranges.append((offset, n))
        offset += n

    if not encoded_parts:
        return torch.empty(0, dtype=torch.long), torch.empty(0, dtype=torch.long)

    all_encoded = torch.cat(encoded_parts)

    # Valid center positions: center kmer starts at c, occupies [c, c+k-1].
    # Left kmer at max distance starts at c - max_dist, needs c - max_dist >= seq_start.
    # Right kmer at max distance starts at c + max_dist, ends at c + max_dist + k - 1,
    # needs c + max_dist + k - 1 < seq_start + seq_len.
    valid_centers = []
    for start, length in seq_ranges:
        first_valid = start + max_dist
        last_valid = start + length - max_dist - kmer_size
        if first_valid <= last_valid:
            valid_centers.append(torch.arange(first_valid, last_valid + 1))

    if not valid_centers:
        return all_encoded, torch.empty(0, dtype=torch.long)

    return all_encoded, torch.cat(valid_centers)


# ---------------------------------------------------------------------------
# Model
# ---------------------------------------------------------------------------


class KmerNet(nn.Module):
    def __init__(self, dims, k, alphabet_size):
        super().__init__()
        self.emb = nn.Embedding(alphabet_size, dims)
        self.L = nn.Sequential(
            nn.Linear(k * dims, k * dims),
            nn.ELU(),
            nn.Linear(k * dims, k * dims),
            nn.ELU(),
            nn.Linear(k * dims, dims),
        )
        self.T = nn.Linear(dims, dims)

    def forward(self, x):
        return torch.softmax(self.L(torch.flatten(self.emb(x), -2)), dim=-1)


# ---------------------------------------------------------------------------
# Loss
# ---------------------------------------------------------------------------


def negative_skipgram_loss(anchor, targets, negatives):
    """Compute negative skip-gram loss.

    Args:
        anchor: [batch, dims] - the center embedding.
        targets: [num_targets, batch, dims] - positive target embeddings.
        negatives: [num_targets, batch, num_neg, dims] - negative embeddings.
    """
    pos_score = (targets * anchor.unsqueeze(0)).sum(dim=-1)
    neg_score = torch.einsum("tbnd,bd->tbn", negatives, anchor)
    loss = -F.logsigmoid(pos_score).mean() - F.logsigmoid(-neg_score).mean()
    return loss


# ---------------------------------------------------------------------------
# Training
# ---------------------------------------------------------------------------


def train(
    fasta_path,
    alphabet="amino",
    kmer_size=3,
    dims=16,
    num_neg=5,
    batch_size=4096,
    epochs=1,
    lr=1e-3,
    ambig_rate=0.0,
    max_steps=0,
    device=None,
):
    lookup_tensor, alphabet_size, _, ambiguate_fn, _ = get_alphabet_config(alphabet)

    if device is None:
        if torch.cuda.is_available():
            device = "cuda"
        elif torch.backends.mps.is_available():
            device = "mps"
        else:
            device = "cpu"
    device = torch.device(device)
    log.info("Using device: %s", device)

    log.info(
        "Encoding sequences from %s (k=%d, alphabet=%s)...",
        fasta_path,
        kmer_size,
        alphabet,
    )
    max_centers = max_steps * batch_size * epochs if max_steps > 0 else 0
    all_encoded, valid_centers = encode_sequences(
        fasta_path, kmer_size, lookup_tensor, max_centers=max_centers
    )
    log.info("Encoded %d total positions, %d valid centers", all_encoded.shape[0], valid_centers.shape[0])

    all_encoded = all_encoded.to(device)
    valid_centers = valid_centers.to(device)
    num_centers = valid_centers.shape[0]
    num_batches_per_epoch = num_centers // batch_size
    if max_steps > 0 and num_batches_per_epoch > max_steps:
        num_batches_per_epoch = max_steps
    max_dist = 3 * kmer_size
    k_offsets = torch.arange(kmer_size, device=device)

    model = KmerNet(dims, kmer_size, alphabet_size).to(device)
    optimizer = torch.optim.Adam(model.parameters(), lr=lr)

    for epoch in range(epochs):
        total_loss = 0.0
        interval_loss = 0.0

        perm = torch.randperm(num_centers, device=device)

        for step in range(num_batches_per_epoch):
            idx = perm[step * batch_size : (step + 1) * batch_size]
            centers = valid_centers[idx]

            left_dist = torch.randint(1, max_dist + 1, (batch_size,), device=device)
            right_dist = torch.randint(1, max_dist + 1, (batch_size,), device=device)

            center_starts = centers.unsqueeze(1) + k_offsets.unsqueeze(0)
            left_starts = (centers - left_dist).unsqueeze(1) + k_offsets.unsqueeze(0)
            right_starts = (centers + right_dist).unsqueeze(1) + k_offsets.unsqueeze(0)

            center_kmers = all_encoded[center_starts]
            left_kmers = all_encoded[left_starts]
            right_kmers = all_encoded[right_starts]

            if ambig_rate > 0:
                left_kmers = ambiguate_fn(left_kmers, ambig_rate)
                right_kmers = ambiguate_fn(right_kmers, ambig_rate)

            left_emb = model(left_kmers)
            center_emb = model(center_kmers)
            right_emb = model(right_kmers)

            norm_loss = torch.median(center_emb, dim=-1, keepdim=True)[0]
            norm_loss = center_emb[center_emb < norm_loss]
            norm_loss = norm_loss.mean()

            left_emb = model.T(left_emb)
            center_emb = model.T(center_emb)
            right_emb = model.T(right_emb)

            neg_left = left_emb[
                torch.randint(0, batch_size, (batch_size, num_neg), device=device)
            ]
            neg_right = right_emb[
                torch.randint(0, batch_size, (batch_size, num_neg), device=device)
            ]

            targets = torch.stack([left_emb, right_emb])
            negatives = torch.stack([neg_left, neg_right])
            loss = negative_skipgram_loss(center_emb, targets, negatives) + norm_loss

            optimizer.zero_grad(set_to_none=True)
            loss.backward()
            optimizer.step()

            total_loss += loss.item()
            interval_loss += loss.item()

            if (step + 1) % 1000 == 0:
                log.info(
                    "Epoch %d/%d step %d/%d - loss: %.4f",
                    epoch + 1,
                    epochs,
                    step + 1,
                    num_batches_per_epoch,
                    interval_loss / 1000,
                )
                interval_loss = 0.0

        avg_loss = total_loss / max(num_batches_per_epoch, 1)
        log.info(
            "Epoch %d/%d - loss: %.4f (%d batches)",
            epoch + 1,
            epochs,
            avg_loss,
            num_batches_per_epoch,
        )

    log.info("Training complete.")
    return model


# ---------------------------------------------------------------------------
# All-kmer generation
# ---------------------------------------------------------------------------


def generate_all_kmers(kmer_size, alphabet_size):
    """Generate all possible kmers for the given alphabet size and k."""
    cols = []
    for pos in range(kmer_size):
        repeat_each = alphabet_size ** (kmer_size - 1 - pos)
        repeat_block = alphabet_size**pos
        col = torch.arange(alphabet_size, dtype=torch.long).repeat_interleave(
            repeat_each
        )
        cols.append(col.repeat(repeat_block))
    return torch.stack(cols, dim=1)


# ---------------------------------------------------------------------------
# Export
# ---------------------------------------------------------------------------


def _compute_all_embeddings(model, kmer_size, alphabet_size, device):
    """Compute model(kmer) for every possible kmer. Returns (all_kmers, embeddings)."""
    all_kmers = generate_all_kmers(kmer_size, alphabet_size)
    total = all_kmers.shape[0]
    log.info("Computing embeddings for %d kmers...", total)

    model.eval()
    chunk_size = 65536
    embeddings = []
    with torch.no_grad():
        for i in range(0, total, chunk_size):
            chunk = all_kmers[i : i + chunk_size].to(device)
            embeddings.append(model(chunk).cpu())
    return all_kmers, torch.cat(embeddings, dim=0)


def export_csv(model, kmer_size, alphabet, output_path, device):
    """Write CSV with kmer string and model(kmer) embeddings for every possible kmer."""
    _, alphabet_size, alphabet_str, _, index_to_char = get_alphabet_config(alphabet)
    all_kmers, embeddings = _compute_all_embeddings(
        model, kmer_size, alphabet_size, device
    )
    total = all_kmers.shape[0]
    dims = embeddings.shape[1]
    log.info("Writing %s (%d kmers, %d dims)...", output_path, total, dims)

    with open(output_path, "w") as f:
        f.write("# alphabet: %s\n" % alphabet_str)
        header = "kmer," + ",".join("dim_%d" % d for d in range(dims))
        f.write(header + "\n")
        for i in range(total):
            kmer_indices = all_kmers[i]
            kmer_str = "".join(index_to_char[j] for j in kmer_indices.tolist())
            vals = ",".join("%.8f" % v for v in embeddings[i].tolist())
            f.write("%s,%s\n" % (kmer_str, vals))

    log.info("Wrote %s", output_path)


def export_bin(model, kmer_size, alphabet, output_path, device):
    """Write binary kmer embedding file.

    Format:
        4 bytes:  magic "KMER"
        4 bytes:  alphabet_size  (uint32 LE)
        4 bytes:  kmer_size      (uint32 LE)
        4 bytes:  dims           (uint32 LE)
        4 bytes:  alphabet_str_len (uint32 LE)
        N bytes:  alphabet string (alphabet_str_len bytes, no null terminator)
        rest:     float32[alphabet_size**kmer_size][dims] (little-endian)

    Kmers are in lexicographic order over alphabet indices, so the index for a
    kmer [i0, i1, ..., i_{k-1}] is:
        i0 * alphabet_size^(k-1) + i1 * alphabet_size^(k-2) + ... + i_{k-1}
    """
    _, alphabet_size, alphabet_str, _, _ = get_alphabet_config(alphabet)
    _, embeddings = _compute_all_embeddings(model, kmer_size, alphabet_size, device)
    total = embeddings.shape[0]
    dims = embeddings.shape[1]
    log.info("Writing %s (%d kmers, %d dims, binary)...", output_path, total, dims)

    alphabet_bytes = alphabet_str.encode("ascii")
    with open(output_path, "wb") as f:
        f.write(b"KMER")
        f.write(
            struct.pack("<IIII", alphabet_size, kmer_size, dims, len(alphabet_bytes))
        )
        f.write(alphabet_bytes)
        # Write embeddings as contiguous little-endian float32
        f.write(embeddings.numpy().astype("<f4").tobytes())

    log.info("Wrote %s", output_path)


# ---------------------------------------------------------------------------
# Save / Load
# ---------------------------------------------------------------------------


def save_model(path, model, alphabet, kmer_size):
    dims = model.T.out_features
    torch.save(
        {
            "dims": dims,
            "k": kmer_size,
            "alphabet": alphabet,
            "alphabet_size": model.emb.num_embeddings,
            "state_dict": model.state_dict(),
        },
        path,
    )
    log.info("Saved model to %s", path)


def load_model(path, device=None):
    if device is None:
        if torch.cuda.is_available():
            device = "cuda"
        elif torch.backends.mps.is_available():
            device = "mps"
        else:
            device = "cpu"
    device = torch.device(device)

    checkpoint = torch.load(path, map_location=device, weights_only=True)
    dims = checkpoint["dims"]
    k = checkpoint["k"]
    alphabet = checkpoint["alphabet"]
    alphabet_size = checkpoint["alphabet_size"]

    model = KmerNet(dims, k, alphabet_size).to(device)
    model.load_state_dict(checkpoint["state_dict"])
    model.eval()

    log.info(
        "Loaded model from %s (dims=%d, k=%d, alphabet=%s)", path, dims, k, alphabet
    )
    return model, alphabet, k, device


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def parse_args(argv=None):
    p = argparse.ArgumentParser(
        description="Train kmer skip-gram embeddings and export kmer hash table.",
    )

    p.add_argument(
        "fasta_file",
        nargs="?",
        default=None,
        help="Input FASTA file (required unless --load-model is given).",
    )

    g = p.add_argument_group("alphabet")
    g.add_argument(
        "--alphabet",
        choices=["amino", "dna"],
        default="amino",
        help="Sequence alphabet (default: amino).",
    )

    g = p.add_argument_group("training")
    g.add_argument("--kmer-size", type=int, default=3, help="Kmer size (default: 3).")
    g.add_argument(
        "--dims", type=int, default=16, help="Embedding dimensions (default: 16)."
    )
    g.add_argument(
        "--epochs", type=int, default=1, help="Training epochs (default: 2)."
    )
    g.add_argument(
        "--batch-size", type=int, default=4096, help="Batch size (default: 4096)."
    )
    g.add_argument(
        "--lr", type=float, default=1e-3, help="Learning rate (default: 1e-3)."
    )
    g.add_argument(
        "--num-neg",
        type=int,
        default=5,
        help="Negative samples per target (default: 5).",
    )
    g.add_argument(
        "--ambig-rate",
        type=float,
        default=0.005,
        help="Ambiguation rate for context kmers (default: 0.0).",
    )
    g.add_argument(
        "--max-steps",
        type=int,
        default=0,
        help="Max training steps (0 = unlimited). Limits data loaded into memory.",
    )
    g.add_argument(
        "--device",
        choices=["cuda", "mps", "cpu"],
        default=None,
        help="Compute device (default: auto-detect).",
    )

    g = p.add_argument_group("model I/O")
    g.add_argument(
        "--save-model",
        default=None,
        metavar="PATH",
        help="Save trained model to this path.",
    )
    g.add_argument(
        "--load-model",
        default=None,
        metavar="PATH",
        help="Load pretrained model (skip training).",
    )

    g = p.add_argument_group("output")
    g.add_argument(
        "--format",
        choices=["csv", "bin"],
        default="bin",
        help="Output format (default: bin).",
    )
    g.add_argument(
        "--output",
        default=None,
        metavar="PATH",
        help="Output path (default: kmer_hashes.bin or .csv).",
    )

    args = p.parse_args(argv)

    if args.load_model is None and args.fasta_file is None:
        p.error("fasta_file is required unless --load-model is given.")

    return args


def main(argv=None):
    args = parse_args(argv)

    if args.load_model is not None:
        model, alphabet, kmer_size, device = load_model(args.load_model, args.device)
    else:
        model = train(
            fasta_path=args.fasta_file,
            alphabet=args.alphabet,
            kmer_size=args.kmer_size,
            dims=args.dims,
            num_neg=args.num_neg,
            batch_size=args.batch_size,
            epochs=args.epochs,
            lr=args.lr,
            ambig_rate=args.ambig_rate,
            max_steps=args.max_steps,
            device=args.device,
        )
        alphabet = args.alphabet
        kmer_size = args.kmer_size

        if args.device is None:
            if torch.cuda.is_available():
                device = torch.device("cuda")
            elif torch.backends.mps.is_available():
                device = torch.device("mps")
            else:
                device = torch.device("cpu")
        else:
            device = torch.device(args.device)

        if args.save_model is not None:
            save_model(args.save_model, model, alphabet, kmer_size)

    output_path = args.output
    if output_path is None:
        output_path = "kmer_hashes." + args.format

    if args.format == "csv":
        export_csv(model, kmer_size, alphabet, output_path, device)
    else:
        export_bin(model, kmer_size, alphabet, output_path, device)


if __name__ == "__main__":
    main()
