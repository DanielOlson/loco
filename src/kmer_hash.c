#include "kmer_hash.h"

int kmer_hash_load(KmerHash *h, const char *path) {
    FILE *f = fopen(path, "rb");
    if (!f) return -1;

    char magic[4];
    if (fread(magic, 1, 4, f) != 4 || memcmp(magic, "KMER", 4) != 0) {
        fclose(f);
        return -1;
    }

    uint32_t header[4];
    if (fread(header, sizeof(uint32_t), 4, f) != 4) {
        fclose(f);
        return -1;
    }
    h->alphabet_size = header[0];
    h->kmer_size     = header[1];
    h->dims          = header[2];
    uint32_t alpha_len = header[3];

    h->alphabet = (char *)malloc(alpha_len + 1);
    if (!h->alphabet) { fclose(f); return -1; }
    if (fread(h->alphabet, 1, alpha_len, f) != alpha_len) {
        free(h->alphabet);
        fclose(f);
        return -1;
    }
    h->alphabet[alpha_len] = '\0';

    /* Compute num_kmers = alphabet_size ^ kmer_size */
    uint32_t num = 1;
    for (uint32_t i = 0; i < h->kmer_size; i++)
        num *= h->alphabet_size;
    h->num_kmers = num;

    size_t data_bytes = (size_t)num * h->dims * sizeof(float);
    h->embeddings = (float *)malloc(data_bytes);
    if (!h->embeddings) {
        free(h->alphabet);
        fclose(f);
        return -1;
    }
    if (fread(h->embeddings, 1, data_bytes, f) != data_bytes) {
        free(h->embeddings);
        free(h->alphabet);
        fclose(f);
        return -1;
    }

    fclose(f);
    return 0;
}

void kmer_hash_free(KmerHash *h) {
    free(h->alphabet);
    free(h->embeddings);
    h->alphabet = NULL;
    h->embeddings = NULL;
}
