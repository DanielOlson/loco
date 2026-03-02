#ifndef LOW_COMPLEXITY_MODELS_H
#define LOW_COMPLEXITY_MODELS_H
#include "common.h"
#include "kmer_hash.h"
#include <stdint.h>

typedef float Entropy;

void fill_count_to_entropy_lookup_table_kmer(Entropy *table,
                                             uint64_t num_kmer_samples);

void fill_complexity_buffer_for_seq_kmer(BioLetter *sequence,
                                         uint64_t sequence_length,
                                         Entropy *complexity_entropy_buffer,
                                         uint64_t entropy_window_size,
                                         uint32_t kmer_size,
                                         Entropy *count_to_entropy_table,
                                         uint32_t spacing);

void fill_complexity_buffer_for_seq_hash(
    BioLetter *sequence, uint64_t sequence_length,
    Entropy *complexity_entropy_buffer, uint64_t entropy_window_size,
    const KmerHash *hash, const int *bio_to_alpha, uint32_t spacing);

#endif
