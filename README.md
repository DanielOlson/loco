# loco

The **LO**w **CO**mplexity annotator.
This is **NOT** ready for broad use. This README is AI-slop. This **IS** ready for careful use. Contact with questions.


loco identifies and masks low-complexity regions in DNA and protein sequences
using multi-channel entropy analysis. It combines traditional k-mer frequency
channels with learned neural embedding channels, and supports metachannel
operations (OR, AND, SUM, W_MIN, W_MAX, W_BLUR) for flexible masking
strategies.

Pre-trained embedding models are available at:
https://huggingface.co/datasets/DanielROlson/loco

## Building

```
git clone https://github.com/DanielOlson/loco.git
cd loco
make
```

## Quick start

Using a config file (recommended):

```
# DNA masking with default config
loco --config dna_default.yaml input.fa

# Protein masking with default config
loco --config amino_default.yaml input.fa
```

Using command-line channels directly:

```
# Single k-mer channel: 2-mers, window of 15 k-mer units, threshold 2.05
loco --ww 2 15 2.05 input.fa

# Add an embedding channel
loco --ww 2 15 2.05 --wh models/dna_k3.bin 14 1.99 input.fa
```

## Usage

```
loco [options] <input.fa>

Options:
  -h, --help         Show help message
  -v, --version      Show version information
  --aa               Amino acid mode (default: DNA)
  -b, --bed          Output BED format (default: FASTA)
  -x                 Hard masking: N for DNA, X for protein (default: lowercase)
  -r, --raw          Output raw entropy values per channel
  -o, --out <path>   Output file (default: stdout)

Channel options:
  --ww <kmer_size> <kmer_units> <threshold> [spacing=N]
                     K-mer entropy channel (repeatable)
  --wh <model.bin> <kmer_units> <threshold> [spacing=N]
                     Embedding entropy channel (repeatable)
  --config <file>    Load channels and settings from a YAML config file
  --ratios <r1,r2,...>
                     Letter frequencies for fit mode
                     DNA order: A,C,G,T (default: 0.3,0.2,0.2,0.3)
                     AA order: A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y
```

## Config files

loco is driven by YAML config files that define channels and settings.
Pre-calibrated configs are provided: `dna_default.yaml` and
`amino_default.yaml`.

### Config format

```yaml
fit: false        # true to enter calibration mode
aa: false         # true for amino acid mode
bed: false        # true for BED output
x_mask: false     # true for hard masking (N/X)
raw: false        # true for raw entropy output
out: output.fa    # output file path (optional)

channels:
  # K-mer entropy channel
  - ww <kmer_size> <kmer_units> <threshold> [spacing=N] [name=xxx]

  # Embedding entropy channel
  - wh <model.bin> <kmer_units> <threshold> [spacing=N] [name=xxx]

  # Metachannels (operate on named channels)
  - OR <member1> <member2> ...
  - AND <member1> <member2> ...
  - SUM <threshold> <member1> <member2> ...
  - W_MIN <member> <window_size> <threshold>
  - W_MAX <member> <window_size> <threshold>
  - W_BLUR <member> <window_size> <threshold>
```

### Channel types

**ww (word/k-mer channels):** Compute Shannon entropy over a sliding window of
k-mer frequencies. Parameters are the k-mer size (e.g. 1, 2, 3, 4), the
window size in k-mer units, and the entropy threshold below which a position is
masked.

**wh (hash/embedding channels):** Compute entropy over learned k-mer
embeddings loaded from a binary model file. Parameters are the model path, the
window size in k-mer units, and the entropy threshold.

**spacing:** Optional spaced k-mer sampling. `spacing=5` samples k-mers at
every 5+1 positions within the window instead of consecutively, useful for
detecting periodic low-complexity patterns.

**name:** Optional channel name used to reference it in metachannels. Unnamed
channels are OR'd together by default.

### Metachannels

Metachannels combine named channels:

- **OR:** Union of masks from member channels
- **AND:** Intersection of masks from member channels
- **SUM:** Masks positions where the sum of member entropy values falls below a
  threshold
- **W_MIN / W_MAX / W_BLUR:** Apply a min/max/moving-average convolution to a
  channel's entropy before thresholding

## Calibration (fit mode)

To calibrate thresholds for a specific masking rate on random sequence:

```yaml
# fit_dna.yaml
fit: true
aa: false

channels:
  - ww 1 21 0.03       # target: mask 3% of random sequence
  - ww 2 15 0.03
  - wh models/dna_k3.bin 14 0.01
  - wh models/dna_k4.bin 13 0.01
```

```
loco --config fit_dna.yaml > dna_default.yaml
```

In fit mode, the threshold field specifies the target masking rate (0.0-1.0)
instead of an absolute entropy value. loco generates a 10M-base random sequence
(using default or custom letter frequencies), computes entropy for each
channel, and outputs a calibrated config with absolute thresholds.

To calibrate against a real sequence instead of random:

```
loco --config fit_dna.yaml reference.fa > calibrated.yaml
```

Custom letter frequencies:

```
loco --config fit_dna.yaml --ratios 0.3,0.2,0.2,0.3 > calibrated.yaml
```

## Output formats

**FASTA (default):** Masked sequence. Lowercase masking by default,
or hard masking with `-x` (N for DNA, X for protein).

```
>seq1
ACGTacgtACGT       # lowercase = low complexity
ACGTNNNNNACGT      # -x flag: N-masked (DNA)
```

**BED (`-b`):** Tab-delimited regions (0-based, half-open).

```
seq1    4    8
```

**Raw (`-r`):** Per-position entropy values for each channel.

```
>seq1
1.386294 2.079442 1.945910
1.386294 2.079442 1.945910
...
```

## Models

Pre-trained embedding models for both DNA and amino acid sequences are
available at:

https://huggingface.co/datasets/DanielROlson/loco

Download the `.bin` files and place them in a `models/` directory (or adjust
paths in your config file). The default configs expect models at
`./py/models/`.

### Training custom models

The `py/` directory contains the training pipeline:

```
# Train all k=1..4 models for DNA and amino acid
./py/train_all.sh <dna_fasta> [aa_fasta]
```

This trains skip-gram neural networks on k-mer contexts and exports binary
embedding files. See `py/generate_kmer_hashes.py` for details.
