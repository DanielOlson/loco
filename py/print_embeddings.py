#!/usr/bin/env python3
"""Print kmer embeddings from a .bin file."""

import argparse
import struct
import sys


def load_bin(path):
    with open(path, "rb") as f:
        magic = f.read(4)
        if magic != b"KMER":
            print("Error: not a KMER bin file", file=sys.stderr)
            sys.exit(1)

        alphabet_size, kmer_size, dims, alpha_len = struct.unpack("<IIII", f.read(16))
        alphabet = f.read(alpha_len).decode("ascii")

        num_kmers = alphabet_size ** kmer_size
        data = f.read(num_kmers * dims * 4)
        embeddings = struct.unpack("<%df" % (num_kmers * dims), data)

    return alphabet_size, kmer_size, dims, alphabet, num_kmers, embeddings


def index_to_kmer_str(index, kmer_size, alphabet_size, alphabet):
    chars = []
    for _ in range(kmer_size):
        chars.append(alphabet[index % alphabet_size])
        index //= alphabet_size
    return "".join(reversed(chars))


def main():
    p = argparse.ArgumentParser(description="Print kmer embeddings from a .bin file.")
    p.add_argument("bin_file", help="Path to .bin kmer hash file")
    p.add_argument("--tsv", action="store_true", help="Tab-separated output")
    args = p.parse_args()

    alphabet_size, kmer_size, dims, alphabet, num_kmers, embeddings = load_bin(
        args.bin_file
    )

    sep = "\t" if args.tsv else ","

    # Header
    print("# alphabet: %s  k: %d  dims: %d  num_kmers: %d" % (
        alphabet, kmer_size, dims, num_kmers
    ))
    print(sep.join(["kmer"] + ["dim_%d" % d for d in range(dims)]))

    # Rows
    for i in range(num_kmers):
        kmer_str = index_to_kmer_str(i, kmer_size, alphabet_size, alphabet)
        row = embeddings[i * dims : (i + 1) * dims]
        vals = sep.join("%.8f" % v for v in row)
        print("%s%s%s" % (kmer_str, sep, vals))


if __name__ == "__main__":
    main()
