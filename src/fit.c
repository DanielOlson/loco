#include "fit.h"
#include "common.h"
#include "io.h"
#include "low_complexity_models.h"
#include "metachannel.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

#define FIT_SEQ_LEN 10000000

/* Default DNA ratios: A=0.3, C=0.2, G=0.2, T=0.3 (60% AT) */
static const float default_dna_ratios[] = {0.3f, 0.2f, 0.2f, 0.3f};

/* Default AA ratios (Swiss-Prot background, ACDEFGHIKLMNPQRSTVWY order) */
static const float default_aa_ratios[] = {
    0.0825f, 0.0137f, 0.0545f, 0.0675f, 0.0386f, 0.0707f, 0.0227f,
    0.0596f, 0.0584f, 0.0966f, 0.0242f, 0.0406f, 0.0470f, 0.0393f,
    0.0553f, 0.0656f, 0.0534f, 0.0687f, 0.0108f, 0.0292f};

/* BioLetter values for DNA letters in ACGT order */
static const BioLetter dna_letters[] = {A, C, G, T};

/* BioLetter values for AA letters in ACDEFGHIKLMNPQRSTVWY order */
static const BioLetter aa_letters[] = {A, C, D, E, F, G, H, I, K, L,
                                       M, N, P, Q, R, S, T, V, W, Y};

/*
 * Parse comma-separated floats from str into out[].
 * Returns the number of values parsed, or -1 on error.
 */
static int parse_ratios(const char *str, float *out, int max_count) {
  int count = 0;
  const char *p = str;
  while (*p && count < max_count) {
    char *end;
    out[count] = strtof(p, &end);
    if (end == p)
      return -1;
    count++;
    if (*end == ',')
      p = end + 1;
    else if (*end == '\0')
      break;
    else
      return -1;
  }
  return count;
}

/*
 * Generate random BioLetter sequence using cumulative distribution.
 */
static void generate_random_sequence(BioLetter *seq, uint64_t len,
                                     const float *ratios,
                                     const BioLetter *letters,
                                     int alphabet_size) {
  /* Build cumulative distribution */
  float *cum = (float *)malloc(sizeof(float) * alphabet_size);
  cum[0] = ratios[0];
  for (int i = 1; i < alphabet_size; i++)
    cum[i] = cum[i - 1] + ratios[i];

  srand(42);
  for (uint64_t i = 0; i < len; i++) {
    float r = (float)rand() / (float)RAND_MAX;
    int idx = alphabet_size - 1;
    for (int j = 0; j < alphabet_size - 1; j++) {
      if (r < cum[j]) {
        idx = j;
        break;
      }
    }
    seq[i] = letters[idx];
  }

  free(cum);
}

/*
 * Comparison function for qsort on Entropy (float) values.
 */
static int compare_entropy(const void *a, const void *b) {
  float fa = *(const float *)a;
  float fb = *(const float *)b;
  if (fa < fb)
    return -1;
  if (fa > fb)
    return 1;
  return 0;
}

int fit_channels(FILE *out, _Bool use_aa, WwEntry *ww, int ww_count,
                 WhEntry *wh, int wh_count, const char *ratios_str,
                 const char *fasta_path,
                 MetaChannel *meta, int meta_count) {
  BioLetter *seq = NULL;
  uint64_t seq_len = 0;

  init_lookup_tables();

  if (fasta_path) {
    /* Load FASTA and concatenate all sequences into one flat array */
    BioSequence *sequences =
        use_aa ? read_fasta_aa((char *)fasta_path)
               : read_fasta_dna((char *)fasta_path);
    if (!sequences) {
      fprintf(stderr, "Error: could not read FASTA file: %s\n", fasta_path);
      return -1;
    }

    /* First pass: total length */
    for (BioSequence *s = sequences; s; s = (BioSequence *)s->next_sequence)
      seq_len += s->sequence_length;

    if (seq_len == 0) {
      fprintf(stderr, "Error: FASTA file contains no sequence data\n");
      return -1;
    }

    seq = (BioLetter *)malloc(sizeof(BioLetter) * seq_len);
    uint64_t offset = 0;
    for (BioSequence *s = sequences; s; s = (BioSequence *)s->next_sequence) {
      memcpy(seq + offset, s->sequence,
             sizeof(BioLetter) * s->sequence_length);
      offset += s->sequence_length;
    }
  } else {
    /* Generate random sequence from ratios */
    int alphabet_size = use_aa ? 20 : 4;
    const float *ratios;
    float parsed_ratios[20];

    if (ratios_str) {
      int n = parse_ratios(ratios_str, parsed_ratios, 20);
      if (n != alphabet_size) {
        fprintf(stderr,
                "Error: --ratios expects %d values for %s mode, got %d\n",
                alphabet_size, use_aa ? "AA" : "DNA", n);
        return -1;
      }
      float sum = 0.0f;
      for (int i = 0; i < n; i++) {
        if (parsed_ratios[i] < 0.0f) {
          fprintf(stderr, "Error: ratios must be non-negative\n");
          return -1;
        }
        sum += parsed_ratios[i];
      }
      if (fabsf(sum - 1.0f) > 0.01f) {
        fprintf(stderr, "Error: ratios must sum to ~1.0 (got %.4f)\n", sum);
        return -1;
      }
      ratios = parsed_ratios;
    } else {
      ratios = use_aa ? default_aa_ratios : default_dna_ratios;
    }

    const BioLetter *letters = use_aa ? aa_letters : dna_letters;
    seq_len = FIT_SEQ_LEN;
    seq = (BioLetter *)malloc(sizeof(BioLetter) * seq_len);
    if (!seq) {
      fprintf(stderr, "Error: could not allocate random sequence\n");
      return -1;
    }
    generate_random_sequence(seq, seq_len, ratios, letters, alphabet_size);
  }

  int total_channels = ww_count + wh_count;
  float *fitted_thresholds = (float *)malloc(sizeof(float) * total_channels);

  /* Compute entropy arrays for all channels (kept for global mask) */
  Entropy **all_entropies =
      (Entropy **)malloc(sizeof(Entropy *) * total_channels);
  int ch = 0;

  /* Process --ww channels */
  for (int w = 0; w < ww_count; w++) {
    all_entropies[ch] = (Entropy *)malloc(sizeof(Entropy) * seq_len);

    int num_kmer_samples;
    int ws;
    if (ww[w].spacing > 0) {
      int num_states = ww[w].spacing + 1;
      num_kmer_samples = ww[w].kmer_units;
      ws = (ww[w].kmer_units - 1) * num_states + ww[w].kmer_size;
    } else {
      ws = ww[w].kmer_units * ww[w].kmer_size;
      num_kmer_samples = ws - (ww[w].kmer_size - 1);
    }
    Entropy *table = (Entropy *)malloc(sizeof(Entropy) * (num_kmer_samples + 1));
    fill_count_to_entropy_lookup_table_kmer(table,
                                            (uint64_t)num_kmer_samples);
    fill_complexity_buffer_for_seq_kmer(seq, seq_len, all_entropies[ch],
                                        (uint64_t)ws,
                                        (uint32_t)ww[w].kmer_size, table,
                                        (uint32_t)ww[w].spacing);
    free(table);

    /* Sort a copy and pick threshold at target quantile */
    Entropy *sorted = (Entropy *)malloc(sizeof(Entropy) * seq_len);
    memcpy(sorted, all_entropies[ch], sizeof(Entropy) * seq_len);
    qsort(sorted, seq_len, sizeof(Entropy), compare_entropy);

    float target = ww[w].threshold;
    uint64_t idx = (uint64_t)floorf(target * (float)seq_len);
    if (idx >= seq_len)
      idx = seq_len - 1;
    fitted_thresholds[ch] = sorted[idx];

    /* Cap threshold to 99.9% of theoretical max entropy ln(num_kmer_samples)
     * so that channels with near-degenerate distributions on the calibration
     * sequence do not produce thresholds at the entropy ceiling. */
    float max_ent = logf((float)num_kmer_samples);
    float cap = max_ent * 0.999f;
    if (fitted_thresholds[ch] > cap)
      fitted_thresholds[ch] = cap;

    free(sorted);
    ch++;
  }

  /* Process --wh channels */
  BioLetter *lookup = use_aa ? char_to_amino : char_to_dna;
  for (int w = 0; w < wh_count; w++) {
    all_entropies[ch] = (Entropy *)malloc(sizeof(Entropy) * seq_len);

    if (kmer_hash_load(&wh[w].hash, wh[w].path) != 0) {
      fprintf(stderr, "Error: could not load bin file: %s\n", wh[w].path);
      for (int c = 0; c <= ch; c++)
        free(all_entropies[c]);
      free(all_entropies);
      free(fitted_thresholds);
      free(seq);
      return -1;
    }

    /* Build bio_to_alpha mapping */
    for (int i = 0; i < 21; i++)
      wh[w].bio_to_alpha[i] = 0;
    const char *alpha = wh[w].hash.alphabet;
    for (uint32_t i = 0; i < wh[w].hash.alphabet_size; i++) {
      BioLetter bl = lookup[(unsigned char)alpha[i]];
      wh[w].bio_to_alpha[bl] = (int)i;
    }

    int ws;
    if (wh[w].spacing > 0) {
      int num_states = wh[w].spacing + 1;
      ws = (wh[w].kmer_units - 1) * num_states + (int)wh[w].hash.kmer_size;
    } else {
      ws = wh[w].kmer_units * (int)wh[w].hash.kmer_size;
    }
    fill_complexity_buffer_for_seq_hash(seq, seq_len, all_entropies[ch],
                                        (uint64_t)ws,
                                        &wh[w].hash, wh[w].bio_to_alpha,
                                        (uint32_t)wh[w].spacing);

    /* Sort a copy and pick threshold at target quantile */
    Entropy *sorted = (Entropy *)malloc(sizeof(Entropy) * seq_len);
    memcpy(sorted, all_entropies[ch], sizeof(Entropy) * seq_len);
    qsort(sorted, seq_len, sizeof(Entropy), compare_entropy);

    float target = wh[w].threshold;
    uint64_t idx = (uint64_t)floorf(target * (float)seq_len);
    if (idx >= seq_len)
      idx = seq_len - 1;
    fitted_thresholds[ch] = sorted[idx];

    free(sorted);
    ch++;
  }

  /* Fit metachannel thresholds (single pass, topo order).
   * W_* convolved buffers are kept so SUM can reference them. */
  Entropy **meta_ent = (Entropy **)calloc(
      meta_count > 0 ? meta_count : 1, sizeof(Entropy *));
  for (int m = 0; m < meta_count; m++) {
    if (meta[m].op == META_W_MIN || meta[m].op == META_W_MAX ||
        meta[m].op == META_W_BLUR) {
      int ci = find_channel_by_name(meta[m].members[0],
                                    ww, ww_count, wh, wh_count);
      if (ci < 0)
        continue;

      Entropy *convolved = (Entropy *)malloc(sizeof(Entropy) * seq_len);
      if (meta[m].op == META_W_MIN)
        convolve_min(all_entropies[ci], convolved, seq_len,
                     meta[m].window_size);
      else if (meta[m].op == META_W_MAX)
        convolve_max(all_entropies[ci], convolved, seq_len,
                     meta[m].window_size);
      else
        convolve_blur(all_entropies[ci], convolved, seq_len,
                      meta[m].window_size);

      Entropy *sorted = (Entropy *)malloc(sizeof(Entropy) * seq_len);
      memcpy(sorted, convolved, sizeof(Entropy) * seq_len);
      qsort(sorted, seq_len, sizeof(Entropy), compare_entropy);

      float target = meta[m].threshold;
      uint64_t idx = (uint64_t)floorf(target * (float)seq_len);
      if (idx >= seq_len)
        idx = seq_len - 1;
      meta[m].threshold = sorted[idx];

      free(sorted);
      meta_ent[m] = convolved;
    } else if (meta[m].op == META_SUM) {
      /* Sum entropies from base channels or W_* metachannels */
      Entropy *sum_entropies = (Entropy *)calloc(seq_len, sizeof(Entropy));
      for (int mb = 0; mb < meta[m].member_count; mb++) {
        int ci = find_channel_by_name(meta[m].members[mb],
                                      ww, ww_count, wh, wh_count);
        if (ci >= 0) {
          for (uint64_t i = 0; i < seq_len; i++)
            sum_entropies[i] += all_entropies[ci][i];
        } else {
          int mi = find_meta_by_name(meta[m].members[mb],
                                     meta, meta_count);
          if (mi >= 0 && meta_ent[mi] != NULL) {
            for (uint64_t i = 0; i < seq_len; i++)
              sum_entropies[i] += meta_ent[mi][i];
          }
        }
      }

      /* Sort and pick threshold at target quantile */
      Entropy *sorted = (Entropy *)malloc(sizeof(Entropy) * seq_len);
      memcpy(sorted, sum_entropies, sizeof(Entropy) * seq_len);
      qsort(sorted, seq_len, sizeof(Entropy), compare_entropy);

      float target = meta[m].threshold;
      uint64_t idx = (uint64_t)floorf(target * (float)seq_len);
      if (idx >= seq_len)
        idx = seq_len - 1;
      meta[m].threshold = sorted[idx];

      free(sorted);
      free(sum_entropies);
    }
  }
  for (int m = 0; m < meta_count; m++)
    free(meta_ent[m]);
  free(meta_ent);

  /* Bulk-compute metachannel masks (each computed once, topo order) */
  _Bool **meta_masks = NULL;
  if (meta_count > 0) {
    meta_masks = compute_meta_masks(meta, meta_count,
                                    ww, ww_count, wh, wh_count,
                                    all_entropies, fitted_thresholds,
                                    seq_len);
  }

  /* Compute global OR mask and per-source unique masking rates.
   * Sources: base channels 0..total_channels-1, then metachannels. */
  int num_sources = total_channels + meta_count;
  uint64_t masked_count = 0;
  uint64_t *unique_counts = (uint64_t *)calloc(num_sources, sizeof(uint64_t));

  for (uint64_t i = 0; i < seq_len; i++) {
    int src_mask_count = 0;
    int src_mask_idx = -1;

    for (int c = 0; c < total_channels; c++) {
      if (all_entropies[c][i] < fitted_thresholds[c]) {
        src_mask_count++;
        src_mask_idx = c;
      }
    }

    for (int m = 0; m < meta_count; m++) {
      if (meta_masks[m][i]) {
        src_mask_count++;
        src_mask_idx = total_channels + m;
      }
    }

    if (src_mask_count > 0)
      masked_count++;
    if (src_mask_count == 1)
      unique_counts[src_mask_idx]++;
  }

  /* Per-metachannel masking rates */
  float *meta_rates = (float *)calloc(meta_count > 0 ? meta_count : 1,
                                       sizeof(float));
  if (meta_masks) {
    for (int m = 0; m < meta_count; m++) {
      uint64_t cnt = 0;
      for (uint64_t i = 0; i < seq_len; i++)
        if (meta_masks[m][i])
          cnt++;
      meta_rates[m] = (float)cnt / (float)seq_len;
    }
    for (int m = 0; m < meta_count; m++)
      free(meta_masks[m]);
    free(meta_masks);
  }

  /* Per-channel individual masking rates (for all base channels) */
  float *channel_rates = (float *)calloc(total_channels, sizeof(float));
  for (int c = 0; c < total_channels; c++) {
    uint64_t cnt = 0;
    for (uint64_t i = 0; i < seq_len; i++)
      if (all_entropies[c][i] < fitted_thresholds[c])
        cnt++;
    channel_rates[c] = (float)cnt / (float)seq_len;
  }

  float global_rate = (float)masked_count / (float)seq_len;

  /* Write YAML config */
  fprintf(out, "# loco --fit calibration\n");
  fprintf(out, "# global masking rate: %.4f\n", global_rate);
  fprintf(out, "fit: false\n");
  fprintf(out, "aa: %s\n", use_aa ? "true" : "false");
  fprintf(out, "\nchannels:\n");

  ch = 0;
  for (int w = 0; w < ww_count; w++) {
    fprintf(out, "  - ww %d %d %.6g",
            ww[w].kmer_size, ww[w].kmer_units, fitted_thresholds[ch]);
    if (ww[w].spacing != 0)
      fprintf(out, " spacing=%d", ww[w].spacing);
    if (ww[w].name[0] != '\0')
      fprintf(out, " name=%s", ww[w].name);
    fprintf(out, "  # rate: %.4f, unique: %.4f\n",
            channel_rates[ch],
            (float)unique_counts[ch] / (float)seq_len);
    ch++;
  }
  for (int w = 0; w < wh_count; w++) {
    fprintf(out, "  - wh %s %d %.6g",
            wh[w].path, wh[w].kmer_units, fitted_thresholds[ch]);
    if (wh[w].spacing != 0)
      fprintf(out, " spacing=%d", wh[w].spacing);
    if (wh[w].name[0] != '\0')
      fprintf(out, " name=%s", wh[w].name);
    fprintf(out, "  # rate: %.4f, unique: %.4f\n",
            channel_rates[ch],
            (float)unique_counts[ch] / (float)seq_len);
    ch++;
  }

  /* Write metachannels inside channels block */
  for (int m = 0; m < meta_count; m++) {
    const char *op_str = "AND";
    if (meta[m].op == META_OR)
      op_str = "OR";
    else if (meta[m].op == META_SUM)
      op_str = "SUM";
    else if (meta[m].op == META_W_MIN)
      op_str = "W_MIN";
    else if (meta[m].op == META_W_MAX)
      op_str = "W_MAX";
    else if (meta[m].op == META_W_BLUR)
      op_str = "W_BLUR";

    fprintf(out, "  - %s", op_str);

    if (meta[m].op == META_W_MIN || meta[m].op == META_W_MAX ||
        meta[m].op == META_W_BLUR) {
      /* W_* format: <op> <member> <window_size> <threshold> [name=xxx] */
      fprintf(out, " %s %d %.6g", meta[m].members[0],
              meta[m].window_size, meta[m].threshold);
      if (meta[m].name[0] != '\0')
        fprintf(out, " name=%s", meta[m].name);
    } else {
      if (meta[m].op == META_SUM)
        fprintf(out, " %.6g", meta[m].threshold);

      if (meta[m].name[0] != '\0')
        fprintf(out, " name=%s", meta[m].name);

      for (int mb = 0; mb < meta[m].member_count; mb++)
        fprintf(out, " %s", meta[m].members[mb]);
    }

    fprintf(out, "  # rate: %.4f, unique: %.4f\n",
            meta_rates[m],
            (float)unique_counts[total_channels + m] / (float)seq_len);
  }
  free(unique_counts);
  free(channel_rates);
  free(meta_rates);

  /* Cleanup */
  for (int c = 0; c < total_channels; c++)
    free(all_entropies[c]);
  free(all_entropies);
  free(fitted_thresholds);
  free(seq);

  for (int w = 0; w < wh_count; w++)
    kmer_hash_free(&wh[w].hash);

  return 0;
}
