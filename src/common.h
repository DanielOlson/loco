#ifndef COMMON_H
#define COMMON_H

#include <stdint.h>

typedef enum {
  A,
  C,
  G,
  T,
  D,
  E,
  F,
  H,
  I,
  K,
  L,
  M,
  N,
  P,
  Q,
  R,
  S,
  V,
  W,
  Y,
  UNDEFINED
} BioLetter;

typedef struct BioSequence {
  BioLetter *sequence;
  char *raw_sequence;
  char *raw_block;
  char *sequence_name;
  uint64_t sequence_length;
  uint64_t raw_block_length;

  struct BioSequence *next_sequence;
} BioSequence;

extern BioLetter char_to_amino[256];
extern BioLetter char_to_dna[256];

void init_lookup_tables(void);

#endif
