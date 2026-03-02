#include "metachannel.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void convolve_min(const Entropy *src, Entropy *dst, uint64_t len, int ws) {
  int half = ws / 2;
  for (uint64_t i = 0; i < len; i++) {
    int64_t lo = (int64_t)i - half;
    int64_t hi = lo + ws;
    if (lo < 0) lo = 0;
    if (hi > (int64_t)len) hi = (int64_t)len;
    Entropy val = src[lo];
    for (int64_t j = lo + 1; j < hi; j++)
      if (src[j] < val) val = src[j];
    dst[i] = val;
  }
}

void convolve_max(const Entropy *src, Entropy *dst, uint64_t len, int ws) {
  int half = ws / 2;
  for (uint64_t i = 0; i < len; i++) {
    int64_t lo = (int64_t)i - half;
    int64_t hi = lo + ws;
    if (lo < 0) lo = 0;
    if (hi > (int64_t)len) hi = (int64_t)len;
    Entropy val = src[lo];
    for (int64_t j = lo + 1; j < hi; j++)
      if (src[j] > val) val = src[j];
    dst[i] = val;
  }
}

void convolve_blur(const Entropy *src, Entropy *dst, uint64_t len, int ws) {
  int half = ws / 2;
  for (uint64_t i = 0; i < len; i++) {
    int64_t lo = (int64_t)i - half;
    int64_t hi = lo + ws;
    if (lo < 0) lo = 0;
    if (hi > (int64_t)len) hi = (int64_t)len;
    Entropy sum = 0.0f;
    for (int64_t j = lo; j < hi; j++)
      sum += src[j];
    dst[i] = sum / (Entropy)(hi - lo);
  }
}

int find_channel_by_name(const char *name,
                         const WwEntry *ww, int ww_count,
                         const WhEntry *wh, int wh_count) {
  for (int i = 0; i < ww_count; i++) {
    if (ww[i].name[0] != '\0' && strcmp(ww[i].name, name) == 0)
      return i;
  }
  for (int i = 0; i < wh_count; i++) {
    if (wh[i].name[0] != '\0' && strcmp(wh[i].name, name) == 0)
      return ww_count + i;
  }
  return -1;
}

int find_meta_by_name(const char *name,
                      const MetaChannel *meta, int meta_count) {
  for (int i = 0; i < meta_count; i++) {
    if (meta[i].name[0] != '\0' && strcmp(meta[i].name, name) == 0)
      return i;
  }
  return -1;
}

/*
 * Get the per-position mask for a member by name.
 * If the member is a base channel, build a mask from entropy < threshold.
 * If the member is a metachannel, return its already-computed mask.
 * Returns a pointer to a _Bool array of length len, or NULL.
 * Sets *allocated = 1 if the caller must free the returned array.
 */
static _Bool *get_member_mask(const char *name,
                              const MetaChannel *meta, int meta_count,
                              const WwEntry *ww, int ww_count,
                              const WhEntry *wh, int wh_count,
                              Entropy **channels, const float *thresholds,
                              _Bool **meta_masks, uint64_t len,
                              int *allocated) {
  *allocated = 0;

  /* Check if it's a metachannel first (meta names checked before base) */
  int mi = find_meta_by_name(name, meta, meta_count);
  if (mi >= 0 && meta_masks[mi] != NULL) {
    return meta_masks[mi];
  }

  /* Must be a base channel */
  int ch = find_channel_by_name(name, ww, ww_count, wh, wh_count);
  if (ch < 0)
    return NULL;

  _Bool *mask = (_Bool *)malloc(sizeof(_Bool) * len);
  if (!mask)
    return NULL;
  *allocated = 1;

  float thr = thresholds[ch];
  Entropy *ent = channels[ch];
  for (uint64_t i = 0; i < len; i++)
    mask[i] = ent[i] < thr;

  return mask;
}

_Bool **compute_meta_masks(const MetaChannel *meta, int meta_count,
                           const WwEntry *ww, int ww_count,
                           const WhEntry *wh, int wh_count,
                           Entropy **channels, const float *thresholds,
                           uint64_t len) {
  _Bool **meta_masks = (_Bool **)calloc(meta_count, sizeof(_Bool *));
  if (!meta_masks)
    return NULL;

  /* Entropy buffers for W_* metachannels so SUM can reference them */
  Entropy **meta_ent = (Entropy **)calloc(meta_count, sizeof(Entropy *));
  if (!meta_ent) {
    free(meta_masks);
    return NULL;
  }

  /* Process in order (topologically sorted: dependencies first) */
  for (int m = 0; m < meta_count; m++) {
    meta_masks[m] = (_Bool *)calloc(len, sizeof(_Bool));
    if (!meta_masks[m])
      goto fail;

    const MetaChannel *mc = &meta[m];

    if (mc->op == META_AND) {
      /* Start with all 1s, AND in each member */
      memset(meta_masks[m], 1, sizeof(_Bool) * len);
      for (int mb = 0; mb < mc->member_count; mb++) {
        int alloc = 0;
        _Bool *member = get_member_mask(mc->members[mb], meta, meta_count,
                                        ww, ww_count, wh, wh_count,
                                        channels, thresholds,
                                        meta_masks, len, &alloc);
        if (!member)
          goto fail;
        for (uint64_t i = 0; i < len; i++)
          meta_masks[m][i] &= member[i];
        if (alloc)
          free(member);
      }
    } else if (mc->op == META_OR) {
      /* Start with all 0s, OR in each member */
      for (int mb = 0; mb < mc->member_count; mb++) {
        int alloc = 0;
        _Bool *member = get_member_mask(mc->members[mb], meta, meta_count,
                                        ww, ww_count, wh, wh_count,
                                        channels, thresholds,
                                        meta_masks, len, &alloc);
        if (!member)
          goto fail;
        for (uint64_t i = 0; i < len; i++)
          meta_masks[m][i] |= member[i];
        if (alloc)
          free(member);
      }
    } else if (mc->op == META_SUM) {
      /* Sum entropies from base channels or W_* metachannels */
      float *sums = (float *)calloc(len, sizeof(float));
      if (!sums)
        goto fail;
      for (int mb = 0; mb < mc->member_count; mb++) {
        int ch = find_channel_by_name(mc->members[mb],
                                      ww, ww_count, wh, wh_count);
        if (ch >= 0) {
          Entropy *ent = channels[ch];
          for (uint64_t i = 0; i < len; i++)
            sums[i] += ent[i];
        } else {
          int mi = find_meta_by_name(mc->members[mb], meta, meta_count);
          if (mi >= 0 && meta_ent[mi] != NULL) {
            for (uint64_t i = 0; i < len; i++)
              sums[i] += meta_ent[mi][i];
          }
        }
      }
      float thr = mc->threshold;
      for (uint64_t i = 0; i < len; i++)
        meta_masks[m][i] = sums[i] < thr;
      free(sums);
    } else if (mc->op == META_W_MIN || mc->op == META_W_MAX ||
               mc->op == META_W_BLUR) {
      int ch = find_channel_by_name(mc->members[0],
                                    ww, ww_count, wh, wh_count);
      if (ch < 0)
        goto fail;
      Entropy *convolved = (Entropy *)malloc(sizeof(Entropy) * len);
      if (!convolved)
        goto fail;
      if (mc->op == META_W_MIN)
        convolve_min(channels[ch], convolved, len, mc->window_size);
      else if (mc->op == META_W_MAX)
        convolve_max(channels[ch], convolved, len, mc->window_size);
      else
        convolve_blur(channels[ch], convolved, len, mc->window_size);
      float thr = mc->threshold;
      for (uint64_t i = 0; i < len; i++)
        meta_masks[m][i] = convolved[i] < thr;
      /* Keep convolved buffer so SUM can reference it */
      meta_ent[m] = convolved;
    }
  }

  for (int i = 0; i < meta_count; i++)
    free(meta_ent[i]);
  free(meta_ent);
  return meta_masks;

fail:
  for (int i = 0; i < meta_count; i++) {
    free(meta_masks[i]);
    free(meta_ent[i]);
  }
  free(meta_masks);
  free(meta_ent);
  return NULL;
}

/*
 * DFS cycle detection + topological sort.
 * color: 0=white (unvisited), 1=gray (in progress), 2=black (done).
 * order[] is filled in topological order (dependencies first).
 */
static int topo_visit(int node, const MetaChannel *meta, int meta_count,
                      const WwEntry *ww, int ww_count,
                      const WhEntry *wh, int wh_count,
                      int *color, int *order, int *order_pos) {
  if (color[node] == 2)
    return 0;
  if (color[node] == 1) {
    const char *name = meta[node].name[0] ? meta[node].name : "(unnamed)";
    fprintf(stderr, "Error: cycle detected involving metachannel '%s'\n", name);
    return -1;
  }

  color[node] = 1;

  for (int m = 0; m < meta[node].member_count; m++) {
    /* Skip base channels */
    if (find_channel_by_name(meta[node].members[m],
                             ww, ww_count, wh, wh_count) >= 0)
      continue;
    int mi = find_meta_by_name(meta[node].members[m], meta, meta_count);
    if (mi < 0)
      continue;
    if (topo_visit(mi, meta, meta_count, ww, ww_count, wh, wh_count,
                   color, order, order_pos) != 0)
      return -1;
  }

  color[node] = 2;
  order[(*order_pos)++] = node;
  return 0;
}

int validate_metachannels(const MetaChannel *meta, int meta_count,
                          const WwEntry *ww, int ww_count,
                          const WhEntry *wh, int wh_count) {
  int total_base = ww_count + wh_count;

  /* Check for duplicate names among base channels */
  for (int i = 0; i < total_base; i++) {
    const char *name_i;
    if (i < ww_count)
      name_i = ww[i].name;
    else
      name_i = wh[i - ww_count].name;
    if (name_i[0] == '\0')
      continue;

    for (int j = i + 1; j < total_base; j++) {
      const char *name_j;
      if (j < ww_count)
        name_j = ww[j].name;
      else
        name_j = wh[j - ww_count].name;
      if (name_j[0] == '\0')
        continue;
      if (strcmp(name_i, name_j) == 0) {
        fprintf(stderr, "Error: duplicate channel name: %s\n", name_i);
        return -1;
      }
    }

    /* Also check against metachannel names */
    for (int j = 0; j < meta_count; j++) {
      if (meta[j].name[0] != '\0' && strcmp(name_i, meta[j].name) == 0) {
        fprintf(stderr,
                "Error: channel and metachannel share name: %s\n", name_i);
        return -1;
      }
    }
  }

  /* Check for duplicate metachannel names */
  for (int i = 0; i < meta_count; i++) {
    if (meta[i].name[0] == '\0')
      continue;
    for (int j = i + 1; j < meta_count; j++) {
      if (meta[j].name[0] != '\0' &&
          strcmp(meta[i].name, meta[j].name) == 0) {
        fprintf(stderr, "Error: duplicate metachannel name: %s\n",
                meta[i].name);
        return -1;
      }
    }
  }

  /* Validate each metachannel's members */
  for (int i = 0; i < meta_count; i++) {
    for (int m = 0; m < meta[i].member_count; m++) {
      const char *mname = meta[i].members[m];
      int ch = find_channel_by_name(mname, ww, ww_count, wh, wh_count);
      int mi = find_meta_by_name(mname, meta, meta_count);

      if (ch < 0 && mi < 0) {
        fprintf(stderr, "Error: metachannel member '%s' not found\n", mname);
        return -1;
      }

      /* SUM members must be base channels or W_* metachannels */
      if (meta[i].op == META_SUM && ch < 0) {
        if (mi < 0 || (meta[mi].op != META_W_MIN &&
                        meta[mi].op != META_W_MAX &&
                        meta[mi].op != META_W_BLUR)) {
          fprintf(stderr,
                  "Error: SUM metachannel member '%s' must be a base channel"
                  " or a W_* metachannel\n", mname);
          return -1;
        }
      }

      /* W_MIN/W_MAX/W_BLUR can only reference base channels */
      if ((meta[i].op == META_W_MIN || meta[i].op == META_W_MAX ||
           meta[i].op == META_W_BLUR) && ch < 0) {
        fprintf(stderr,
                "Error: W_* metachannel member '%s' must be a base channel\n",
                mname);
        return -1;
      }
    }
  }

  /* Validate W_* metachannels: exactly 1 member, window_size >= 1 */
  for (int i = 0; i < meta_count; i++) {
    if (meta[i].op == META_W_MIN || meta[i].op == META_W_MAX ||
        meta[i].op == META_W_BLUR) {
      if (meta[i].member_count != 1) {
        fprintf(stderr,
                "Error: W_* metachannel requires exactly 1 member\n");
        return -1;
      }
      if (meta[i].window_size < 1) {
        fprintf(stderr,
                "Error: W_* metachannel window_size must be >= 1\n");
        return -1;
      }
    }
  }

  /* Cycle detection via DFS */
  if (meta_count > 0) {
    int color[MAX_CHANNELS] = {0};
    int order[MAX_CHANNELS];
    int order_pos = 0;

    for (int i = 0; i < meta_count; i++) {
      if (color[i] == 0) {
        if (topo_visit(i, meta, meta_count, ww, ww_count, wh, wh_count,
                       color, order, &order_pos) != 0)
          return -1;
      }
    }
  }

  /* Warn about orphaned named channels (named but not in any metachannel) */
  for (int i = 0; i < total_base; i++) {
    const char *name;
    if (i < ww_count)
      name = ww[i].name;
    else
      name = wh[i - ww_count].name;
    if (name[0] == '\0')
      continue;

    int found = 0;
    for (int j = 0; j < meta_count && !found; j++) {
      for (int m = 0; m < meta[j].member_count; m++) {
        if (strcmp(meta[j].members[m], name) == 0) {
          found = 1;
          break;
        }
      }
    }
    if (!found) {
      fprintf(stderr, "Warning: named channel '%s' is not used by any "
                      "metachannel (will not contribute to masking)\n", name);
    }
  }

  return 0;
}

int sort_metachannels(MetaChannel *meta, int meta_count,
                      const WwEntry *ww, int ww_count,
                      const WhEntry *wh, int wh_count) {
  if (meta_count <= 1)
    return 0;

  int color[MAX_CHANNELS] = {0};
  int order[MAX_CHANNELS];
  int order_pos = 0;

  for (int i = 0; i < meta_count; i++) {
    if (color[i] == 0) {
      if (topo_visit(i, meta, meta_count, ww, ww_count, wh, wh_count,
                     color, order, &order_pos) != 0)
        return -1;
    }
  }

  MetaChannel sorted[MAX_CHANNELS];
  for (int i = 0; i < meta_count; i++)
    sorted[i] = meta[order[i]];
  memcpy(meta, sorted, sizeof(MetaChannel) * meta_count);

  return 0;
}
