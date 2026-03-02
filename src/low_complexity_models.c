#include "low_complexity_models.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

#define NUM_LETTERS 21

void fill_count_to_entropy_lookup_table_kmer(Entropy *table,
                                             uint64_t num_kmer_samples) {
  table[0] = 0.0;
  for (uint64_t i = 1; i <= num_kmer_samples; ++i) {
    double p = (double)i / (double)num_kmer_samples;
    table[i] = (Entropy)(-p * log(p));
  }
}

static inline uint32_t compute_kmer_index(BioLetter *seq, uint64_t pos,
                                           uint32_t kmer_size) {
  uint32_t idx = 0;
  for (uint32_t j = 0; j < kmer_size; j++)
    idx = idx * NUM_LETTERS + (uint32_t)seq[pos + j];
  return idx;
}

static uint32_t num_possible_kmers(uint32_t kmer_size) {
  uint32_t n = 1;
  for (uint32_t i = 0; i < kmer_size; i++)
    n *= NUM_LETTERS;
  return n;
}

void fill_complexity_buffer_for_seq_kmer(BioLetter *sequence,
                                         uint64_t sequence_length,
                                         Entropy *complexity_entropy_buffer,
                                         uint64_t entropy_window_size,
                                         uint32_t kmer_size,
                                         Entropy *count_to_entropy,
                                         uint32_t spacing) {
  uint32_t total_kmers = num_possible_kmers(kmer_size);
  uint64_t half_window = entropy_window_size / 2;

  if (spacing == 0) {
    /* Original code path: consecutive kmers */
    int *kmer_to_count = (int *)calloc(total_kmers, sizeof(int));
    Entropy entropy_sum = 0.0;
    uint64_t kmers_in_window = entropy_window_size - (kmer_size - 1);

    for (uint64_t i = 0; i < kmers_in_window; ++i) {
      uint32_t idx = compute_kmer_index(sequence, i, kmer_size);
      kmer_to_count[idx]++;
    }

    for (uint32_t i = 0; i < total_kmers; ++i) {
      entropy_sum += count_to_entropy[kmer_to_count[i]];
    }

    for (uint64_t i = 0; i <= half_window; ++i) {
      complexity_entropy_buffer[i] = entropy_sum;
    }

    for (uint64_t i = entropy_window_size; i < sequence_length; ++i) {
      uint32_t old_kmer = compute_kmer_index(sequence,
                                              i - entropy_window_size,
                                              kmer_size);
      uint32_t new_kmer = compute_kmer_index(sequence,
                                              i - kmer_size + 1, kmer_size);

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

    for (uint64_t i = sequence_length - half_window; i < sequence_length; ++i) {
      complexity_entropy_buffer[i] = entropy_sum;
    }

    free(kmer_to_count);
  } else {
    /* Multi-state sliding window for spaced kmers */
    uint32_t num_states = spacing + 1;
    uint32_t kmer_units = (uint32_t)((entropy_window_size - kmer_size) /
                                      num_states) + 1;

    for (uint64_t i = 0; i < sequence_length; i++)
      complexity_entropy_buffer[i] = 0.0f;

    for (uint32_t s = 0; s < num_states; s++) {
      int *kmer_to_count = (int *)calloc(total_kmers, sizeof(int));
      Entropy entropy_sum = 0.0;

      /* Check if this state has enough sequence for a full window */
      uint64_t last_kmer_start = s + (uint64_t)(kmer_units - 1) * num_states;
      if (last_kmer_start + kmer_size > sequence_length)
        goto next_state;

      /* Init first window: kmer starts at s, s+num_states, ... */
      for (uint32_t u = 0; u < kmer_units; u++) {
        uint64_t pos = s + (uint64_t)u * num_states;
        uint32_t idx = compute_kmer_index(sequence, pos, kmer_size);
        kmer_to_count[idx]++;
      }

      for (uint32_t i = 0; i < total_kmers; ++i)
        entropy_sum += count_to_entropy[kmer_to_count[i]];

      /* Fill early positions for this state up to center */
      {
        uint64_t center = s + half_window;
        for (uint64_t pos = s; pos <= center && pos < sequence_length;
             pos += num_states) {
          complexity_entropy_buffer[pos] = entropy_sum;
        }
      }

      /* Slide: for each step t, remove oldest kmer and add new one */
      uint32_t t = 1;
      while (1) {
        uint64_t new_kmer_start =
            s + (uint64_t)(t + kmer_units - 1) * num_states;
        if (new_kmer_start + kmer_size > sequence_length)
          break;

        uint64_t old_kmer_start = s + (uint64_t)(t - 1) * num_states;
        uint32_t old_kmer = compute_kmer_index(sequence, old_kmer_start,
                                                kmer_size);
        uint32_t new_kmer = compute_kmer_index(sequence, new_kmer_start,
                                                kmer_size);

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

        /* Write centered output */
        uint64_t out_pos = s + half_window + (uint64_t)t * num_states;
        if (out_pos < sequence_length)
          complexity_entropy_buffer[out_pos] = entropy_sum;

        t++;
      }

      /* Fill trailing positions for this state */
      {
        uint64_t last_center = s + half_window + (uint64_t)(t - 1) * num_states;
        for (uint64_t pos = last_center + num_states; pos < sequence_length;
             pos += num_states) {
          complexity_entropy_buffer[pos] = entropy_sum;
        }
      }

    next_state:
      free(kmer_to_count);
    }

  }
}

void fill_complexity_buffer_for_seq_hash(
    BioLetter *sequence, uint64_t sequence_length,
    Entropy *complexity_entropy_buffer, uint64_t entropy_window_size,
    const KmerHash *hash, const int *bio_to_alpha, uint32_t spacing) {

  uint32_t k = hash->kmer_size;
  uint32_t dims = hash->dims;
  uint64_t half_window = entropy_window_size / 2;

  if (spacing == 0) {
    /* Original code path */
    uint64_t kmers_in_window = entropy_window_size - (k - 1);

    if (sequence_length < entropy_window_size || kmers_in_window == 0) {
      for (uint64_t i = 0; i < sequence_length; i++)
        complexity_entropy_buffer[i] = 0.0f;
      return;
    }

    float *state = (float *)calloc(dims, sizeof(float));
    float inv_kmers = 1.0f / (float)kmers_in_window;

    uint64_t num_kmer_positions = sequence_length - (k - 1);
    uint32_t *kmer_indices =
        (uint32_t *)malloc(sizeof(uint32_t) * num_kmer_positions);
    for (uint64_t p = 0; p < num_kmer_positions; p++) {
      uint32_t idx = 0;
      for (uint32_t j = 0; j < k; j++) {
        idx = idx * hash->alphabet_size +
              (uint32_t)bio_to_alpha[sequence[p + j]];
      }
      kmer_indices[p] = idx;
    }

    for (uint64_t i = 0; i < kmers_in_window; i++) {
      const float *emb = kmer_hash_get(hash, kmer_indices[i]);
      for (uint32_t d = 0; d < dims; d++)
        state[d] += emb[d];
    }
    for (uint32_t d = 0; d < dims; d++)
      state[d] *= inv_kmers;

    Entropy ent = 0.0f;
    for (uint32_t d = 0; d < dims; d++) {
      if (state[d] > 0.0f)
        ent -= state[d] * logf(state[d]);
    }

    for (uint64_t i = 0; i <= half_window; i++)
      complexity_entropy_buffer[i] = ent;

    for (uint64_t i = entropy_window_size; i < sequence_length; i++) {
      uint64_t new_kmer_pos = i - (k - 1);
      uint64_t old_kmer_pos = i - entropy_window_size;

      const float *new_emb = kmer_hash_get(hash, kmer_indices[new_kmer_pos]);
      const float *old_emb = kmer_hash_get(hash, kmer_indices[old_kmer_pos]);

      for (uint32_t d = 0; d < dims; d++)
        state[d] += (new_emb[d] - old_emb[d]) * inv_kmers;

      ent = 0.0f;
      for (uint32_t d = 0; d < dims; d++) {
        if (state[d] > 0.0f)
          ent -= state[d] * logf(state[d]);
      }

      complexity_entropy_buffer[i - half_window] = ent;
    }

    for (uint64_t i = sequence_length - half_window; i < sequence_length; i++)
      complexity_entropy_buffer[i] = ent;

    free(state);
    free(kmer_indices);
  } else {
    /* Multi-state sliding window for spaced embeddings */
    uint32_t num_states = spacing + 1;
    uint32_t kmer_units = (uint32_t)((entropy_window_size - k) /
                                      num_states) + 1;
    float inv_kmers = 1.0f / (float)kmer_units;

    /* Precompute kmer indices for all positions */
    uint64_t num_kmer_positions = (sequence_length >= k) ?
                                   sequence_length - (k - 1) : 0;
    if (num_kmer_positions == 0) {
      for (uint64_t i = 0; i < sequence_length; i++)
        complexity_entropy_buffer[i] = 0.0f;
      return;
    }

    uint32_t *kmer_indices =
        (uint32_t *)malloc(sizeof(uint32_t) * num_kmer_positions);
    for (uint64_t p = 0; p < num_kmer_positions; p++) {
      uint32_t idx = 0;
      for (uint32_t j = 0; j < k; j++) {
        idx = idx * hash->alphabet_size +
              (uint32_t)bio_to_alpha[sequence[p + j]];
      }
      kmer_indices[p] = idx;
    }

    for (uint64_t i = 0; i < sequence_length; i++)
      complexity_entropy_buffer[i] = 0.0f;

    for (uint32_t s = 0; s < num_states; s++) {
      float *embed_state = (float *)calloc(dims, sizeof(float));

      /* Check if this state has enough sequence for a full window */
      uint64_t last_kmer_start = s + (uint64_t)(kmer_units - 1) * num_states;
      if (last_kmer_start + k > sequence_length)
        goto next_hash_state;

      /* Init first window */
      for (uint32_t u = 0; u < kmer_units; u++) {
        uint64_t pos = s + (uint64_t)u * num_states;
        const float *emb = kmer_hash_get(hash, kmer_indices[pos]);
        for (uint32_t d = 0; d < dims; d++)
          embed_state[d] += emb[d];
      }
      for (uint32_t d = 0; d < dims; d++)
        embed_state[d] *= inv_kmers;

      /* Compute initial entropy */
      Entropy ent = 0.0f;
      for (uint32_t d = 0; d < dims; d++) {
        if (embed_state[d] > 0.0f)
          ent -= embed_state[d] * logf(embed_state[d]);
      }

      /* Fill early positions for this state */
      {
        uint64_t center = s + half_window;
        for (uint64_t pos = s; pos <= center && pos < sequence_length;
             pos += num_states) {
          complexity_entropy_buffer[pos] = ent;
        }
      }

      /* Slide */
      {
        uint32_t t = 1;
        while (1) {
          uint64_t new_kmer_start =
              s + (uint64_t)(t + kmer_units - 1) * num_states;
          if (new_kmer_start + k > sequence_length)
            break;

          uint64_t old_kmer_start = s + (uint64_t)(t - 1) * num_states;
          const float *new_emb =
              kmer_hash_get(hash, kmer_indices[new_kmer_start]);
          const float *old_emb =
              kmer_hash_get(hash, kmer_indices[old_kmer_start]);

          for (uint32_t d = 0; d < dims; d++)
            embed_state[d] += (new_emb[d] - old_emb[d]) * inv_kmers;

          ent = 0.0f;
          for (uint32_t d = 0; d < dims; d++) {
            if (embed_state[d] > 0.0f)
              ent -= embed_state[d] * logf(embed_state[d]);
          }

          uint64_t out_pos = s + half_window + (uint64_t)t * num_states;
          if (out_pos < sequence_length)
            complexity_entropy_buffer[out_pos] = ent;

          t++;
        }

        /* Fill trailing positions for this state */
        uint64_t last_center =
            s + half_window + (uint64_t)(t - 1) * num_states;
        for (uint64_t pos = last_center + num_states; pos < sequence_length;
             pos += num_states) {
          complexity_entropy_buffer[pos] = ent;
        }
      }

    next_hash_state:
      free(embed_state);
    }

    free(kmer_indices);
  }
}
