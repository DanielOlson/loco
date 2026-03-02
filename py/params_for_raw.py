import sys
import numpy as np

a = []
with open(sys.argv[1]) as file:
    next(file)
    for line in file:
        a.append([float(f) for f in line.rstrip().split(" ")])
a = np.array(a)
for pct in [1, 2.5, 5]:
    vals = np.percentile(a, pct, axis=0)
    print(f"{pct}%: {' '.join(f'{v:.6f}' for v in vals)}")
