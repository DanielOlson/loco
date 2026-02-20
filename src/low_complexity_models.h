#ifndef LOW_COMPLEXITY_MODELS_H
#define LOW_COMPLEXITY_MODELS_H
#include "common.h"
#include <stdint.h>

typedef float Entropy;

void fill_count_to_entropy_lookup_table(Entropy *table,
                                        uint64_t entropy_window_size);

void fill_complexity_buffer_for_seq(BioLetter *sequence,
                                    uint64_t sequence_length,
                                    Entropy *complexity_entropy_buffer,
                                    uint64_t entropy_window_size,
                                    Entropy *count_to_entropy_table);

void fill_count_to_entropy_lookup_table_2mer(Entropy *table,
                                             uint64_t entropy_window_size);

void fill_count_to_entropy_lookup_table_3mer(Entropy *table,
                                             uint64_t entropy_window_size);

void fill_complexity_buffer_for_seq_2mer(BioLetter *sequence,
                                         uint64_t sequence_length,
                                         Entropy *complexity_entropy_buffer,
                                         uint64_t entropy_window_size,
                                         Entropy *count_to_entropy_table);

void fill_complexity_buffer_for_seq_3mer(BioLetter *sequence,
                                         uint64_t sequence_length,
                                         Entropy *complexity_entropy_buffer,
                                         uint64_t entropy_window_size,
                                         Entropy *count_to_entropy_table);

void fill_count_to_entropy_lookup_table_4mer(Entropy *table,
                                             uint64_t entropy_window_size);

void fill_complexity_buffer_for_seq_4mer(BioLetter *sequence,
                                         uint64_t sequence_length,
                                         Entropy *complexity_entropy_buffer,
                                         uint64_t entropy_window_size,
                                         Entropy *count_to_entropy_table);

#endif
