#include "low_complexity_models.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

#define NUM_LETTERS 21
#define NUM_4MERS (NUM_LETTERS * NUM_LETTERS * NUM_LETTERS * NUM_LETTERS)

void fill_count_to_entropy_lookup_table(Entropy *table,
                                        uint64_t entropy_window_size) {
  double number_of_letters = entropy_window_size;
  table[0] = 0.0;
  for (uint64_t i = 1; i < entropy_window_size; ++i) {
    double count = i;
    double p = count / number_of_letters;
    table[i] = -p * log(p);
  }
}

void fill_count_to_entropy_lookup_table_2mer(Entropy *table,
                                             uint64_t entropy_window_size) {
  double num_kmers = entropy_window_size - 1;
  table[0] = 0.0;
  for (uint64_t i = 1; i < entropy_window_size; ++i) {
    double p = (double)i / num_kmers;
    table[i] = -p * log(p);
  }
}

void fill_count_to_entropy_lookup_table_3mer(Entropy *table,
                                             uint64_t entropy_window_size) {
  double num_kmers = entropy_window_size - 2;
  table[0] = 0.0;
  for (uint64_t i = 1; i < entropy_window_size; ++i) {
    double p = (double)i / num_kmers;
    table[i] = -p * log(p);
  }
}

void fill_count_to_entropy_lookup_table_4mer(Entropy *table,
                                             uint64_t entropy_window_size) {
  double num_kmers = entropy_window_size - 3;
  table[0] = 0.0;
  for (uint64_t i = 1; i < entropy_window_size; ++i) {
    double p = (double)i / num_kmers;
    table[i] = -p * log(p);
  }
}

void fill_complexity_buffer_for_seq(BioLetter *sequence,
                                    uint64_t sequence_length,
                                    Entropy *complexity_entropy_buffer,
                                    uint64_t entropy_window_size,
                                    Entropy *count_to_entropy) {
  int letter_to_count[21] = {0};
  uint64_t half_window = entropy_window_size / 2;
  Entropy entropy_sum = 0.0;

  // Count the first entropy_window_size letters
  for (uint64_t i = 0; i < entropy_window_size; ++i) {
    letter_to_count[sequence[i]]++;
  }

  // Calculate the current entropy
  for (int i = 0; i < 21; ++i) {
    entropy_sum += count_to_entropy[letter_to_count[i]];
  }

  // Now that the first full window has been calculated
  // we can fill entropy[0] to entropy[half_window] as entropy_sum
  for (uint64_t i = 0; i <= half_window; ++i) {
    complexity_entropy_buffer[i] = entropy_sum;
  }

  // Fill the entropy buffer normally (add sequence[i] and remove
  // sequence[i-entropy_window_size])
  for (uint64_t i = entropy_window_size; i < sequence_length; ++i) {
    if (sequence[i] != sequence[i - entropy_window_size]) {
      // Remove the old letter
      int old_letter = sequence[i - entropy_window_size];
      entropy_sum -= count_to_entropy[letter_to_count[old_letter]];
      letter_to_count[old_letter]--;
      entropy_sum += count_to_entropy[letter_to_count[old_letter]];

      // Add the new letter
      int new_letter = sequence[i];
      entropy_sum -= count_to_entropy[letter_to_count[new_letter]];
      letter_to_count[new_letter]++;
      entropy_sum += count_to_entropy[letter_to_count[new_letter]];
      if (entropy_sum < 0.0f)
        entropy_sum = 0.0f;
    }
    complexity_entropy_buffer[i - half_window] = entropy_sum;
  }

  // Fill the last half_window of the entropy buffer
  for (uint64_t i = sequence_length - half_window; i < sequence_length; ++i) {
    complexity_entropy_buffer[i] = entropy_sum;
  }
}

void fill_complexity_buffer_for_seq_2mer(BioLetter *sequence,
                                         uint64_t sequence_length,
                                         Entropy *complexity_entropy_buffer,
                                         uint64_t entropy_window_size,
                                         Entropy *count_to_entropy) {
  int kmer_to_count[NUM_LETTERS * NUM_LETTERS];
  memset(kmer_to_count, 0, sizeof(kmer_to_count));
  uint64_t half_window = entropy_window_size / 2;
  Entropy entropy_sum = 0.0;

  // Count 2-mers in the first window
  for (uint64_t i = 0; i + 1 < entropy_window_size; ++i) {
    int idx = sequence[i] * NUM_LETTERS + sequence[i + 1];
    kmer_to_count[idx]++;
  }

  // Calculate the current entropy
  for (int i = 0; i < NUM_LETTERS * NUM_LETTERS; ++i) {
    entropy_sum += count_to_entropy[kmer_to_count[i]];
  }

  // Fill entropy[0] to entropy[half_window]
  for (uint64_t i = 0; i <= half_window; ++i) {
    complexity_entropy_buffer[i] = entropy_sum;
  }

  // Slide window: remove 2-mer starting at old left edge, add new 2-mer
  // ending at new right edge
  for (uint64_t i = entropy_window_size; i < sequence_length; ++i) {
    int old_kmer = sequence[i - entropy_window_size] * NUM_LETTERS +
                   sequence[i - entropy_window_size + 1];
    int new_kmer = sequence[i - 1] * NUM_LETTERS + sequence[i];

    if (old_kmer != new_kmer) {
      entropy_sum -= count_to_entropy[kmer_to_count[old_kmer]];
      kmer_to_count[old_kmer]--;
      entropy_sum += count_to_entropy[kmer_to_count[old_kmer]];

      entropy_sum -= count_to_entropy[kmer_to_count[new_kmer]];
      kmer_to_count[new_kmer]++;
      entropy_sum += count_to_entropy[kmer_to_count[new_kmer]];
      if (entropy_sum < 0.0f)
        entropy_sum = 0.0f;
    }
    complexity_entropy_buffer[i - half_window] = entropy_sum;
  }

  // Fill the last half_window of the entropy buffer
  for (uint64_t i = sequence_length - half_window; i < sequence_length; ++i) {
    complexity_entropy_buffer[i] = entropy_sum;
  }
}

void fill_complexity_buffer_for_seq_3mer(BioLetter *sequence,
                                         uint64_t sequence_length,
                                         Entropy *complexity_entropy_buffer,
                                         uint64_t entropy_window_size,
                                         Entropy *count_to_entropy) {
  int kmer_to_count[NUM_LETTERS * NUM_LETTERS * NUM_LETTERS];
  memset(kmer_to_count, 0, sizeof(kmer_to_count));
  uint64_t half_window = entropy_window_size / 2;
  Entropy entropy_sum = 0.0;

  // Count 3-mers in the first window
  for (uint64_t i = 0; i + 2 < entropy_window_size; ++i) {
    int idx = sequence[i] * NUM_LETTERS * NUM_LETTERS +
              sequence[i + 1] * NUM_LETTERS + sequence[i + 2];
    kmer_to_count[idx]++;
  }

  // Calculate the current entropy
  for (int i = 0; i < NUM_LETTERS * NUM_LETTERS * NUM_LETTERS; ++i) {
    entropy_sum += count_to_entropy[kmer_to_count[i]];
  }

  // Fill entropy[0] to entropy[half_window]
  for (uint64_t i = 0; i <= half_window; ++i) {
    complexity_entropy_buffer[i] = entropy_sum;
  }

  // Slide window: remove 3-mer starting at old left edge, add new 3-mer
  // ending at new right edge
  for (uint64_t i = entropy_window_size; i < sequence_length; ++i) {
    int old_kmer = sequence[i - entropy_window_size] * NUM_LETTERS * NUM_LETTERS +
                   sequence[i - entropy_window_size + 1] * NUM_LETTERS +
                   sequence[i - entropy_window_size + 2];
    int new_kmer = sequence[i - 2] * NUM_LETTERS * NUM_LETTERS +
                   sequence[i - 1] * NUM_LETTERS + sequence[i];

    if (old_kmer != new_kmer) {
      entropy_sum -= count_to_entropy[kmer_to_count[old_kmer]];
      kmer_to_count[old_kmer]--;
      entropy_sum += count_to_entropy[kmer_to_count[old_kmer]];

      entropy_sum -= count_to_entropy[kmer_to_count[new_kmer]];
      kmer_to_count[new_kmer]++;
      entropy_sum += count_to_entropy[kmer_to_count[new_kmer]];
      if (entropy_sum < 0.0f)
        entropy_sum = 0.0f;
    }
    complexity_entropy_buffer[i - half_window] = entropy_sum;
  }

  // Fill the last half_window of the entropy buffer
  for (uint64_t i = sequence_length - half_window; i < sequence_length; ++i) {
    complexity_entropy_buffer[i] = entropy_sum;
  }
}

void fill_complexity_buffer_for_seq_4mer(BioLetter *sequence,
                                         uint64_t sequence_length,
                                         Entropy *complexity_entropy_buffer,
                                         uint64_t entropy_window_size,
                                         Entropy *count_to_entropy) {
  int *kmer_to_count = (int *)calloc(NUM_4MERS, sizeof(int));
  uint64_t half_window = entropy_window_size / 2;
  Entropy entropy_sum = 0.0;

  // Count 4-mers in the first window
  for (uint64_t i = 0; i + 3 < entropy_window_size; ++i) {
    int idx = sequence[i] * NUM_LETTERS * NUM_LETTERS * NUM_LETTERS +
              sequence[i + 1] * NUM_LETTERS * NUM_LETTERS +
              sequence[i + 2] * NUM_LETTERS + sequence[i + 3];
    kmer_to_count[idx]++;
  }

  // Calculate the current entropy
  for (int i = 0; i < NUM_4MERS; ++i) {
    entropy_sum += count_to_entropy[kmer_to_count[i]];
  }

  // Fill entropy[0] to entropy[half_window]
  for (uint64_t i = 0; i <= half_window; ++i) {
    complexity_entropy_buffer[i] = entropy_sum;
  }

  // Slide window
  for (uint64_t i = entropy_window_size; i < sequence_length; ++i) {
    int old_kmer =
        sequence[i - entropy_window_size] * NUM_LETTERS * NUM_LETTERS * NUM_LETTERS +
        sequence[i - entropy_window_size + 1] * NUM_LETTERS * NUM_LETTERS +
        sequence[i - entropy_window_size + 2] * NUM_LETTERS +
        sequence[i - entropy_window_size + 3];
    int new_kmer = sequence[i - 3] * NUM_LETTERS * NUM_LETTERS * NUM_LETTERS +
                   sequence[i - 2] * NUM_LETTERS * NUM_LETTERS +
                   sequence[i - 1] * NUM_LETTERS + sequence[i];

    if (old_kmer != new_kmer) {
      entropy_sum -= count_to_entropy[kmer_to_count[old_kmer]];
      kmer_to_count[old_kmer]--;
      entropy_sum += count_to_entropy[kmer_to_count[old_kmer]];

      entropy_sum -= count_to_entropy[kmer_to_count[new_kmer]];
      kmer_to_count[new_kmer]++;
      entropy_sum += count_to_entropy[kmer_to_count[new_kmer]];
      if (entropy_sum < 0.0f)
        entropy_sum = 0.0f;
    }
    complexity_entropy_buffer[i - half_window] = entropy_sum;
  }

  // Fill the last half_window of the entropy buffer
  for (uint64_t i = sequence_length - half_window; i < sequence_length; ++i) {
    complexity_entropy_buffer[i] = entropy_sum;
  }

  free(kmer_to_count);
}
