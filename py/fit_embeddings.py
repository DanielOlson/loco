import logging
import os

import torch
import torch.nn as nn
import torch.nn.functional as F
from Bio import SeqIO

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(message)s",
)
log = logging.getLogger(__name__)


# Map all amino acids (standard, ambiguous, nonstandard) to integers 0-26.
# Nonstandard aminos are mapped to their closest standard amino.
# Standard: 0-19, Ambiguous: 20-26
_AMINO_TO_INT = {}

_STANDARD = "ACDEFGHIKLMNPQRSTVWY"
for i, aa in enumerate(_STANDARD):
    _AMINO_TO_INT[aa] = i

# Ambiguous amino acids
_AMINO_TO_INT["B"] = 20  # Asx (D or N)
_AMINO_TO_INT["Z"] = 21  # Glx (E or Q)
_AMINO_TO_INT["J"] = 22  # Xle (I or L)
_AMINO_TO_INT["X"] = 23  # Unknown
_AMINO_TO_INT["*"] = 24  # Stop/translation stop
_AMINO_TO_INT["-"] = 25  # Gap
_AMINO_TO_INT["."] = 25  # Gap (alternate)

# Nonstandard -> closest standard
_AMINO_TO_INT["U"] = _AMINO_TO_INT["C"]  # Selenocysteine -> Cysteine
_AMINO_TO_INT["O"] = _AMINO_TO_INT["K"]  # Pyrrolysine -> Lysine

# Build a lookup table (256 entries) for fast ord-based indexing.
# Default unknown chars map to X (23).
_LOOKUP = [23] * 256
for aa, idx in _AMINO_TO_INT.items():
    _LOOKUP[ord(aa)] = idx
    _LOOKUP[ord(aa.lower())] = idx

_LOOKUP_TENSOR = torch.tensor(_LOOKUP, dtype=torch.long)


def extract_kmer_trios(file_path, kmer_size):
    """Extract all consecutive [kmer1, kmer2, kmer3] from sequences in a FASTA file.

    Args:
        file_path: Path to a FASTA file containing protein sequences.
        kmer_size: Size of each kmer.

    Returns:
        Torch long tensor of shape [num_kmer_trios, 3, kmer_size] where each
        value is an integer-encoded amino acid.
    """
    trio_len = kmer_size * 3
    chunks = []

    for record in SeqIO.parse(file_path, "fasta"):
        seq = str(record.seq).upper()
        n = len(seq)
        if n < trio_len:
            continue

        # Encode entire sequence at once via lookup table
        ords = torch.tensor([ord(c) for c in seq], dtype=torch.long)
        encoded = _LOOKUP_TENSOR[ords]

        # Extract all sliding windows of size trio_len, then reshape into 3 kmers
        windows = encoded.unfold(0, trio_len, 1)  # [num_windows, trio_len]
        chunks.append(windows.reshape(-1, 3, kmer_size))

    if not chunks:
        return torch.empty(0, 3, kmer_size, dtype=torch.long)

    return torch.cat(chunks, dim=0)


# For each standard amino index, list the appropriate ambiguous replacements.
# B(20) covers D(2),N(11); Z(21) covers E(3),Q(13); J(22) covers I(7),L(9);
# X(23) covers all standard aminos.
# Rows: amino index, Cols: up to 2 replacement options (padded with -1).
# For each of 26 amino indices, a single replacement value.
# Aminos with a specific ambiguous code get that; others get X(23).
# A second pass randomly picks between the specific code and X for those that
# have both, but we avoid per-element branching by building a flat table
# and using a random coin flip.
_AMBIG_SPECIFIC = [23] * 26  # default: X
_AMBIG_SPECIFIC[2] = 20  # D -> B
_AMBIG_SPECIFIC[11] = 20  # N -> B
_AMBIG_SPECIFIC[3] = 21  # E -> Z
_AMBIG_SPECIFIC[13] = 21  # Q -> Z
_AMBIG_SPECIFIC[7] = 22  # I -> J
_AMBIG_SPECIFIC[9] = 22  # L -> J
_AMBIG_SPECIFIC = torch.tensor(_AMBIG_SPECIFIC, dtype=torch.long)

# Which aminos have two options (specific + X)?
_HAS_SPECIFIC = torch.zeros(26, dtype=torch.bool)
for _idx in (2, 3, 7, 9, 11, 13):
    _HAS_SPECIFIC[_idx] = True

_ambig_cache = {}


def _get_ambig_tables(device):
    if device not in _ambig_cache:
        _ambig_cache[device] = (
            _AMBIG_SPECIFIC.to(device),
            _HAS_SPECIFIC.to(device),
        )
    return _ambig_cache[device]


def ambiguate_kmer_batch(batch, rate):
    """Randomly replace amino acids with appropriate ambiguous codes.

    Args:
        batch: Long tensor of integer-encoded amino acids (any shape).
        rate: Probability of replacing each position.

    Returns:
        New tensor with the same shape, with random positions replaced by
        an appropriate ambiguous amino acid code.
    """
    dev = batch.device
    ambig_specific, has_specific = _get_ambig_tables(dev)

    mask = torch.rand(batch.shape, device=dev) < rate
    mask = mask & (batch < 20)

    result = batch.clone()
    masked_vals = batch[mask]

    # Look up the specific ambiguous code only for masked positions
    replacement = ambig_specific[masked_vals]

    # For aminos that have both a specific code and X, randomly pick X instead
    # with 50% probability
    use_x = has_specific[masked_vals] & (
        torch.rand(masked_vals.shape, device=dev) < 0.5
    )
    replacement[use_x] = 23

    result[mask] = replacement
    return result


def generate_all_kmers(kmer_size):
    """Generate all possible kmers for a given k.

    Args:
        kmer_size: Length of each kmer.

    Returns:
        Long tensor of shape [26**kmer_size, kmer_size] containing every
        combination of the 26 amino acid indices.
    """
    # Each column cycles through 0-25 at a different frequency
    total = 26**kmer_size
    cols = []
    for pos in range(kmer_size):
        repeat_each = 26 ** (kmer_size - 1 - pos)
        repeat_block = 26**pos
        col = torch.arange(26, dtype=torch.long).repeat_interleave(repeat_each)
        cols.append(col.repeat(repeat_block))
    return torch.stack(cols, dim=1)


class KmerNet(nn.Module):
    def __init__(self, dims, k):
        super().__init__()
        self.emb = nn.Embedding(26, dims)
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


def negative_skipgram_loss(anchor, targets, negatives):
    """Compute negative skip-gram loss.

    Args:
        anchor: [batch, dims] - the center embedding.
        targets: [num_targets, batch, dims] - positive target embeddings.
        negatives: [num_targets, batch, num_neg, dims] - negative embeddings.

    Returns:
        Scalar loss.
    """
    # anchor: [B, D], targets: [T, B, D] -> pos_score: [T, B]
    pos_score = (targets * anchor.unsqueeze(0)).sum(dim=-1)
    # negatives: [T, B, N, D], anchor: [B, D] -> neg_score: [T, B, N]
    neg_score = torch.einsum("tbnd,bd->tbn", negatives, anchor)
    loss = -F.logsigmoid(pos_score).mean() - F.logsigmoid(-neg_score).mean()
    return loss


def train(
    fasta_path=None,
    kmer_size=3,
    dims=16,
    num_neg=5,
    batch_size=4096,
    epochs=2,
    lr=1e-3,
    ambig_rate=0.000,
    device=None,
):
    if fasta_path is None:
        fasta_path = os.path.join(os.path.dirname(__file__), "uniprot_sprot.fasta")

    if device is None:
        if torch.cuda.is_available():
            device = "cuda"
        elif torch.backends.mps.is_available():
            device = "mps"
        else:
            device = "cpu"
    device = torch.device(device)
    log.info("Using device: %s", device)

    # Extract kmer trios
    log.info("Extracting kmer trios from %s (k=%d)...", fasta_path, kmer_size)
    trios = extract_kmer_trios(fasta_path, kmer_size)
    log.info("Extracted %d kmer trios", trios.shape[0])

    # Keep entire dataset on device to avoid per-batch transfers (big win on MPS)
    trios = trios.to(device)
    num_trios = trios.shape[0]
    num_batches_per_epoch = num_trios // batch_size

    # 3 separate KmerNets: left context, center, right context
    left_net = KmerNet(dims, kmer_size).to(device)
    center_net = KmerNet(dims, kmer_size).to(device)
    right_net = KmerNet(dims, kmer_size).to(device)
    params = (
        list(left_net.parameters())
        + list(center_net.parameters())
        + list(right_net.parameters())
    )
    optimizer = torch.optim.Adam(params, lr=lr)

    for epoch in range(epochs):
        total_loss = 0.0
        interval_loss = 0.0

        # Shuffle indices on device
        perm = torch.randperm(num_trios, device=device)

        for step in range(num_batches_per_epoch):
            idx = perm[step * batch_size : (step + 1) * batch_size]
            batch = trios[idx]  # [B, 3, kmer_size], already on device

            center_kmers = batch[:, 1, :]
            if ambig_rate > 0:
                ctx = ambiguate_kmer_batch(batch[:, [0, 2], :], ambig_rate)
                left_kmers = ctx[:, 0, :]
                right_kmers = ctx[:, 1, :]
            else:
                left_kmers = batch[:, 0, :]
                right_kmers = batch[:, 2, :]

            left_emb = center_net(left_kmers)  # [B, dims]
            center_emb = center_net(center_kmers)  # [B, dims]
            right_emb = center_net(right_kmers)  # [B, dims]

            left_emb = center_net.T(left_emb)
            center_emb = center_net.T(center_emb)
            right_emb = center_net.T(right_emb)

            # Sample negatives for left and right targets
            neg_left = left_emb[
                torch.randint(0, batch_size, (batch_size, num_neg), device=device)
            ]
            neg_right = right_emb[
                torch.randint(0, batch_size, (batch_size, num_neg), device=device)
            ]

            # Center predicts left and right contexts
            targets = torch.stack([left_emb, right_emb])  # [2, B, dims]
            negatives = torch.stack([neg_left, neg_right])  # [2, B, num_neg, dims]
            loss = negative_skipgram_loss(center_emb, targets, negatives)

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
    return left_net, center_net, right_net


def save_models(path, left_net, center_net, right_net):
    """Save the three KmerNet models to a single file."""
    dims = left_net.L.out_features
    k = left_net.L.in_features // dims
    torch.save(
        {
            "dims": dims,
            "k": k,
            "left": left_net.state_dict(),
            "center": center_net.state_dict(),
            "right": right_net.state_dict(),
        },
        path,
    )
    log.info("Saved models to %s", path)


def load_models(path, device=None):
    """Load three KmerNet models from a file.

    Args:
        path: Path to the saved checkpoint.
        device: Device to load onto. If None, auto-detects.

    Returns:
        Tuple of (left_net, center_net, right_net).
    """
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

    left_net = KmerNet(dims, k).to(device)
    center_net = KmerNet(dims, k).to(device)
    right_net = KmerNet(dims, k).to(device)

    left_net.load_state_dict(checkpoint["left"])
    center_net.load_state_dict(checkpoint["center"])
    right_net.load_state_dict(checkpoint["right"])

    left_net.eval()
    center_net.eval()
    right_net.eval()

    log.info("Loaded models from %s (dims=%d, k=%d)", path, dims, k)
    return left_net, center_net, right_net


if __name__ == "__main__":
    left, center, right = train()
    save_models(
        os.path.join(os.path.dirname(__file__), "kmer_nets.pt"),
        left,
        center,
        right,
    )
