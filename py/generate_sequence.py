import numpy as np

# Standard amino acids and their background frequencies from UniProt/Swiss-Prot
aminos = list("ACDEFGHIKLMNPQRSTVWY")
p = [
    0.0826,  # A
    0.0137,  # C
    0.0546,  # D
    0.0675,  # E
    0.0386,  # F
    0.0708,  # G
    0.0227,  # H
    0.0597,  # I
    0.0584,  # K
    0.0966,  # L
    0.0242,  # M
    0.0406,  # N
    0.0470,  # P
    0.0393,  # Q
    0.0553,  # R
    0.0655,  # S
    0.0534,  # T
    0.0687,  # V
    0.0108,  # W
    0.0292,  # Y
]
# Normalize so probabilities sum to exactly 1.0
p = np.array(p)
p = p / p.sum()


with open('random_prot.fa', 'w') as file:
    file.write(">seq\n")
    for _ in range(10000):
        file.write(''.join(np.random.choice(aminos,p=p,size=(60,))))
        file.write('\n')