#ifndef IO_H
#define IO_H
#include "common.h"
#include "low_complexity_models.h"
#include <stdbool.h>
#include <stdio.h>

BioSequence *read_fasta_aa(char *file_path);
BioSequence *read_fasta_dna(char *file_path);

void output_BED_for_entropies(FILE *out_file, BioSequence *sequence,
                              Entropy *complexity_entropies, Entropy cutoff);

void output_entropies(FILE *out_file, BioSequence *sequence,
                      Entropy *complexity_entropies);

void output_bed_from_mask(FILE *out_file, BioSequence *sequence, _Bool *mask);

void output_fasta_masked(FILE *out_file, BioSequence *sequence, _Bool *mask,
                         _Bool use_x, _Bool use_aa);

void output_raw_entropies(FILE *out_file, BioSequence *sequence,
                          Entropy **ent_channels, int num_channels);

#endif
