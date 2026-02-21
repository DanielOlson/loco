# loco
The LOw COmplexity annotator.

This is not quite ready for broad use.
------------
To build:
```
git clone https://github.com/DanielOlson/loco.git
cd loco
make
```

Usage:
```
Usage: loco [options] <input.fa>

Options:
  -h, --help    Show this help message
  --w1 <value>  1-mer window size (default: 20)
  --w2 <value>  2-mer window size (default: 30)
  --w3 <value>  3-mer window size (default: 40)
  --w4 <value>  4-mer window size (default: 50)
  --aa  Amino acid mode (default: DNA)
  -1 <value>    1-mer entropy threshold (default: 1.00)
  -2 <value>    2-mer entropy threshold (default: 2.05)
  -3 <value>    3-mer entropy threshold (default: 2.90)
  -4 <value>    4-mer entropy threshold (default: 3.47)
  -b, --bed     Output BED format (default: FASTA)
  -x    X-masking (default: lowercase)
  -r, --raw     Output raw entropy values
  -o, --out <value>     Output file path (default: stdout)
```
Default thresholds are meant to mask 1% of random (nonbiological) sequence with 60% AT-richness.

For protein (masks 1%): `--aa -1 2.010597 -2 3.110749 -3 3.564625 -4 3.820651`
For protein (masks 2.5%): `--aa -1 2.067322 -2 3.158371 -3 3.601106 -4 3.850146`
for protein (masks 5%): `--aa -1 2.123229 -2 3.176482 -3 3.601106 -4 3.850146`
Maybe suggested for protein? `--aa -1 2.067322 -2 3.158371 -3 3.564625 -4 3.820651`
