#ifndef CONFIG_H
#define CONFIG_H

#include "kmer_hash.h"

#define MAX_CHANNELS 16
#define MAX_META_MEMBERS 16

typedef struct {
  int kmer_size;
  int kmer_units;
  float threshold;
  int spacing;
  char name[64];
} WwEntry;

typedef struct {
  char path[1024];
  int kmer_units;
  float threshold;
  int spacing;
  KmerHash hash;
  int bio_to_alpha[21];
  char name[64];
} WhEntry;

typedef enum {
  META_AND,
  META_OR,
  META_SUM,
  META_W_MIN,
  META_W_MAX,
  META_W_BLUR
} MetaOp;

typedef struct {
  MetaOp op;
  char name[64];
  float threshold;
  char members[MAX_META_MEMBERS][64];
  int member_count;
  int window_size;
} MetaChannel;

typedef struct {
  /* Settings (-1 = unset, let CLI decide) */
  int aa;
  int bed;
  int raw;
  int x_mask;
  int fit;
  char out[1024];

  /* Channels */
  WwEntry ww_entries[MAX_CHANNELS];
  int ww_count;
  WhEntry wh_entries[MAX_CHANNELS];
  int wh_count;

  /* Metachannels */
  MetaChannel meta_entries[MAX_CHANNELS];
  int meta_count;
} Config;

/*
 * Load a YAML-like config file into cfg.
 * Returns 0 on success, -1 on parse error (with stderr message).
 */
int config_load(Config *cfg, const char *path);

#endif
