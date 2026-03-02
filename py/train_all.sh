#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
FASTA_AA="${SCRIPT_DIR}/uniprot_sprot.fasta"
FASTA_DNA="${1:-}"
OUT_DIR="${SCRIPT_DIR}/models"

if [ -z "$FASTA_DNA" ]; then
    echo "Usage: $0 <dna_fasta_file> [aa_fasta_file]"
    echo "  aa_fasta_file defaults to ${FASTA_AA}"
    exit 1
fi

if [ "${2:-}" != "" ]; then
    FASTA_AA="$2"
fi

mkdir -p "$OUT_DIR"

MAX_STEPS=75000

for k in 1 2 3 4; do
    echo "=== Training amino k=${k} ==="
    python3 "${SCRIPT_DIR}/generate_kmer_hashes.py" "$FASTA_AA" \
        --alphabet amino \
        --kmer-size "$k" \
        --max-steps "$MAX_STEPS" \
        --save-model "${OUT_DIR}/amino_k${k}.pt" \
        --output "${OUT_DIR}/amino_k${k}.bin"

    echo "=== Training dna k=${k} ==="
    python3 "${SCRIPT_DIR}/generate_kmer_hashes.py" "$FASTA_DNA" \
        --alphabet dna \
        --kmer-size "$k" \
        --max-steps "$MAX_STEPS" \
        --save-model "${OUT_DIR}/dna_k${k}.pt" \
        --output "${OUT_DIR}/dna_k${k}.bin"
done

echo "Done. Models and bins written to ${OUT_DIR}/"
