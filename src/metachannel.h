#ifndef METACHANNEL_H
#define METACHANNEL_H

#include "config.h"
#include "low_complexity_models.h"

/*
 * Find a base channel (ww or wh) by name.
 * Returns flat channel index (0..ww_count+wh_count-1) or -1 if not found.
 * ww channels occupy indices 0..ww_count-1, wh channels ww_count..total-1.
 */
int find_channel_by_name(const char *name,
                         const WwEntry *ww, int ww_count,
                         const WhEntry *wh, int wh_count);

/*
 * Find a metachannel by name.
 * Returns index into meta array, or -1 if not found.
 */
int find_meta_by_name(const char *name,
                      const MetaChannel *meta, int meta_count);

/*
 * Validate all metachannels: check member names resolve, SUM members are
 * base channels only, no duplicate names, detect cycles, warn on orphaned
 * named channels.
 * Returns 0 on success, -1 on error (with stderr messages).
 */
int validate_metachannels(const MetaChannel *meta, int meta_count,
                          const WwEntry *ww, int ww_count,
                          const WhEntry *wh, int wh_count);

/*
 * Topologically sort metachannels in-place so dependencies come before
 * dependents. Must be called after validate_metachannels (assumes no cycles).
 * Returns 0 on success, -1 on error.
 */
int sort_metachannels(MetaChannel *meta, int meta_count,
                      const WwEntry *ww, int ww_count,
                      const WhEntry *wh, int wh_count);

/*
 * Bulk-compute per-position mask arrays for all metachannels.
 * meta_masks[m] is a _Bool array of length len for metachannel m.
 * Metachannels must be topologically sorted so that dependencies are
 * computed before dependents. Each metachannel is evaluated once over
 * the full sequence, reusing already-computed masks for member
 * metachannels.
 *
 * Caller must free each meta_masks[m] and the meta_masks array itself.
 * Returns the allocated array, or NULL on allocation failure.
 */
_Bool **compute_meta_masks(const MetaChannel *meta, int meta_count,
                           const WwEntry *ww, int ww_count,
                           const WhEntry *wh, int wh_count,
                           Entropy **channels, const float *thresholds,
                           uint64_t len);

/*
 * Sliding-window convolution helpers for W_MIN/W_MAX/W_BLUR metachannels.
 * Each applies a centered window of size ws over src, writing results to dst.
 * Edges are clamped. len is the array length.
 */
void convolve_min(const Entropy *src, Entropy *dst, uint64_t len, int ws);
void convolve_max(const Entropy *src, Entropy *dst, uint64_t len, int ws);
void convolve_blur(const Entropy *src, Entropy *dst, uint64_t len, int ws);

#endif
