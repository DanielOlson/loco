# loco
The LOw COmplexity annotator.

This is not quite ready for broad use.

```
Usage: ./loco [options] <input.fa>

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
