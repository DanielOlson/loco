#include "common.h"

BioLetter char_to_amino[256];
BioLetter char_to_dna[256];

void init_lookup_tables(void) {
  for (int i = 0; i < 256; i++) {
    char_to_amino[i] = UNDEFINED;
    char_to_dna[i] = UNDEFINED;
  }

  const char *aminos = "ACDEFGHIKLMNPQRSTVWY";
  const BioLetter amino_vals[] = {A, C, D, E, F, G, H, I, K, L,
                                  M, N, P, Q, R, S, T, V, W, Y};
  for (int i = 0; aminos[i]; i++) {
    char_to_amino[(unsigned char)aminos[i]] = amino_vals[i];
    char_to_amino[(unsigned char)aminos[i] + 32] = amino_vals[i];
  }

  const char *bases = "ACGT";
  const BioLetter base_vals[] = {A, C, G, T};
  for (int i = 0; bases[i]; i++) {
    char_to_dna[(unsigned char)bases[i]] = base_vals[i];
    char_to_dna[(unsigned char)bases[i] + 32] = base_vals[i];
  }
  char_to_dna['U'] = T;
  char_to_dna['u'] = T;
}
