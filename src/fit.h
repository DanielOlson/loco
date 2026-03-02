#ifndef FIT_H
#define FIT_H

#include "config.h"
#include <stdio.h>

/*
 * Generate random sequence (or use supplied FASTA), compute entropy for each
 * channel, find thresholds that achieve the target masking rates (stored in
 * each entry's .threshold), and write a YAML config to out.
 *
 * If fasta_path is non-NULL, reads that file instead of generating random
 * sequence. ratios_str and fasta_path are mutually exclusive.
 *
 * meta/meta_count provide metachannels for SUM fitting and YAML output.
 *
 * Returns 0 on success, -1 on error.
 */
int fit_channels(FILE *out, _Bool use_aa,
                 WwEntry *ww, int ww_count,
                 WhEntry *wh, int wh_count,
                 const char *ratios_str,
                 const char *fasta_path,
                 MetaChannel *meta, int meta_count);

#endif
