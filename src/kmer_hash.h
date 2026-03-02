#ifndef KMER_HASH_H
#define KMER_HASH_H

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct {
    uint32_t alphabet_size;
    uint32_t kmer_size;
    uint32_t dims;
    char *alphabet;
    float *embeddings; /* [alphabet_size^kmer_size][dims], row-major */
    uint32_t num_kmers;
} KmerHash;

/*
 * Compute the flat index for a kmer given as an array of alphabet indices.
 * kmer[i] should be in [0, h->alphabet_size).
 */
static inline uint32_t kmer_hash_index(const KmerHash *h, const uint8_t *kmer) {
    uint32_t idx = 0;
    for (uint32_t i = 0; i < h->kmer_size; i++)
        idx = idx * h->alphabet_size + kmer[i];
    return idx;
}

/*
 * Return a pointer to the dims-length embedding for the given kmer index.
 */
static inline const float *kmer_hash_get(const KmerHash *h, uint32_t kmer_idx) {
    return h->embeddings + (size_t)kmer_idx * h->dims;
}

/*
 * Look up the embedding for a kmer given as alphabet indices.
 */
static inline const float *kmer_hash_lookup(const KmerHash *h,
                                            const uint8_t *kmer) {
    return kmer_hash_get(h, kmer_hash_index(h, kmer));
}

/*
 * Load a .bin kmer hash file. Returns 0 on success, -1 on error.
 * The caller must eventually call kmer_hash_free().
 */
int kmer_hash_load(KmerHash *h, const char *path);

void kmer_hash_free(KmerHash *h);

#endif
