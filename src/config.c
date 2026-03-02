#include "config.h"
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static char *trim(char *s) {
  while (*s && isspace((unsigned char)*s))
    s++;
  char *end = s + strlen(s);
  while (end > s && isspace((unsigned char)end[-1]))
    end--;
  *end = '\0';
  return s;
}

static int parse_bool(const char *val, int *out) {
  if (strcmp(val, "true") == 0 || strcmp(val, "1") == 0 ||
      strcmp(val, "yes") == 0) {
    *out = 1;
    return 0;
  }
  if (strcmp(val, "false") == 0 || strcmp(val, "0") == 0 ||
      strcmp(val, "no") == 0) {
    *out = 0;
    return 0;
  }
  return -1;
}

/*
 * Scan the remainder of a line for a name=xxx token.
 * If found, copy the name into dst (up to dst_size-1 chars).
 * Otherwise set dst[0] = '\0'.
 */
static void parse_optional_name(const char *line, char *dst, int dst_size) {
  dst[0] = '\0';
  const char *p = line;
  while (*p) {
    while (*p && isspace((unsigned char)*p))
      p++;
    if (strncmp(p, "name=", 5) == 0) {
      p += 5;
      int i = 0;
      while (*p && !isspace((unsigned char)*p) && i < dst_size - 1) {
        dst[i++] = *p++;
      }
      dst[i] = '\0';
      return;
    }
    /* skip this token */
    while (*p && !isspace((unsigned char)*p))
      p++;
  }
}

/*
 * Scan the remainder of a line for a spacing=N token.
 * If found, return the integer value. Otherwise return 0.
 */
static int parse_optional_int(const char *line, const char *key) {
  size_t klen = strlen(key);
  const char *p = line;
  while (*p) {
    while (*p && isspace((unsigned char)*p))
      p++;
    if (strncmp(p, key, klen) == 0 && p[klen] == '=') {
      return atoi(p + klen + 1);
    }
    while (*p && !isspace((unsigned char)*p))
      p++;
  }
  return 0;
}

int config_load(Config *cfg, const char *path) {
  FILE *f = fopen(path, "r");
  if (!f) {
    fprintf(stderr, "Error: could not open config file: %s\n", path);
    return -1;
  }

  cfg->aa = -1;
  cfg->bed = -1;
  cfg->raw = -1;
  cfg->x_mask = -1;
  cfg->fit = -1;
  cfg->out[0] = '\0';
  cfg->ww_count = 0;
  cfg->wh_count = 0;
  cfg->meta_count = 0;

  char line[2048];
  int in_channels = 0;
  int lineno = 0;

  while (fgets(line, sizeof(line), f)) {
    lineno++;
    char *s = trim(line);

    if (*s == '\0' || *s == '#')
      continue;

    /* Check for list item */
    if (s[0] == '-' && s[1] == ' ') {
      if (!in_channels) {
        fprintf(stderr, "Error: %s:%d: list item outside channels block\n",
                path, lineno);
        fclose(f);
        return -1;
      }

      char *item = trim(s + 2);

      /* Strip trailing comment (# ...) */
      char *hash = strchr(item, '#');
      if (hash) {
        *hash = '\0';
        item = trim(item);
      }

      /* Metachannel: AND, OR, SUM, W_MIN, W_MAX, W_BLUR */
      if (strncmp(item, "AND ", 4) == 0 ||
          strncmp(item, "OR ", 3) == 0 ||
          strncmp(item, "SUM ", 4) == 0 ||
          strncmp(item, "W_MIN ", 6) == 0 ||
          strncmp(item, "W_MAX ", 6) == 0 ||
          strncmp(item, "W_BLUR ", 7) == 0) {
        if (cfg->meta_count >= MAX_CHANNELS) {
          fprintf(stderr, "Error: %s:%d: too many metachannels (max %d)\n",
                  path, lineno, MAX_CHANNELS);
          fclose(f);
          return -1;
        }

        MetaChannel *mc = &cfg->meta_entries[cfg->meta_count];
        mc->name[0] = '\0';
        mc->threshold = 0.0f;
        mc->member_count = 0;
        mc->window_size = 0;

        /* Tokenize item */
        char buf[2048];
        strncpy(buf, item, sizeof(buf) - 1);
        buf[sizeof(buf) - 1] = '\0';

        char *tokens[64];
        int ntokens = 0;
        char *tok = strtok(buf, " \t");
        while (tok && ntokens < 64) {
          tokens[ntokens++] = tok;
          tok = strtok(NULL, " \t");
        }

        /* First token: operation */
        if (strcmp(tokens[0], "AND") == 0)
          mc->op = META_AND;
        else if (strcmp(tokens[0], "OR") == 0)
          mc->op = META_OR;
        else if (strcmp(tokens[0], "SUM") == 0)
          mc->op = META_SUM;
        else if (strcmp(tokens[0], "W_MIN") == 0)
          mc->op = META_W_MIN;
        else if (strcmp(tokens[0], "W_MAX") == 0)
          mc->op = META_W_MAX;
        else
          mc->op = META_W_BLUR;

        /* Remaining tokens: name=xxx, threshold (for SUM), member names */
        int ti = 1;

        /* Check for name=xxx among remaining tokens */
        for (int t = 1; t < ntokens; t++) {
          if (strncmp(tokens[t], "name=", 5) == 0) {
            strncpy(mc->name, tokens[t] + 5, sizeof(mc->name) - 1);
            mc->name[sizeof(mc->name) - 1] = '\0';
            for (int u = t; u < ntokens - 1; u++)
              tokens[u] = tokens[u + 1];
            ntokens--;
            break;
          }
        }

        /* For W_MIN/W_MAX/W_BLUR: <op> <member> <window_size> <threshold> */
        if (mc->op == META_W_MIN || mc->op == META_W_MAX ||
            mc->op == META_W_BLUR) {
          if (ntokens - ti < 3) {
            fprintf(stderr,
                    "Error: %s:%d: %s requires <member> <window_size>"
                    " <threshold>\n",
                    path, lineno, tokens[0]);
            fclose(f);
            return -1;
          }
          /* token[ti] = member name */
          strncpy(mc->members[0], tokens[ti], 63);
          mc->members[0][63] = '\0';
          mc->member_count = 1;
          ti++;
          /* token[ti] = window_size */
          mc->window_size = atoi(tokens[ti]);
          ti++;
          /* token[ti] = threshold */
          char *endp;
          mc->threshold = strtof(tokens[ti], &endp);
          if (endp == tokens[ti]) {
            fprintf(stderr,
                    "Error: %s:%d: %s: invalid threshold '%s'\n",
                    path, lineno, tokens[0], tokens[ti]);
            fclose(f);
            return -1;
          }
          ti++;
          cfg->meta_count++;
          continue;
        }

        /* For SUM, next token after op is threshold */
        if (mc->op == META_SUM) {
          if (ti >= ntokens) {
            fprintf(stderr,
                    "Error: %s:%d: SUM requires a threshold\n",
                    path, lineno);
            fclose(f);
            return -1;
          }
          char *endp;
          mc->threshold = strtof(tokens[ti], &endp);
          if (endp == tokens[ti]) {
            fprintf(stderr,
                    "Error: %s:%d: SUM: invalid threshold '%s'\n",
                    path, lineno, tokens[ti]);
            fclose(f);
            return -1;
          }
          ti++;
        }

        /* Remaining tokens are member names */
        for (; ti < ntokens; ti++) {
          if (mc->member_count >= MAX_META_MEMBERS) {
            fprintf(stderr,
                    "Error: %s:%d: too many metachannel members (max %d)\n",
                    path, lineno, MAX_META_MEMBERS);
            fclose(f);
            return -1;
          }
          strncpy(mc->members[mc->member_count], tokens[ti], 63);
          mc->members[mc->member_count][63] = '\0';
          mc->member_count++;
        }

        if (mc->member_count < 2) {
          fprintf(stderr,
                  "Error: %s:%d: metachannel requires at least 2 members\n",
                  path, lineno);
          fclose(f);
          return -1;
        }

        cfg->meta_count++;
      } else if (strncmp(item, "ww ", 3) == 0) {
        if (cfg->ww_count + cfg->wh_count >= MAX_CHANNELS) {
          fprintf(stderr, "Error: %s:%d: too many channels (max %d)\n", path,
                  lineno, MAX_CHANNELS);
          fclose(f);
          return -1;
        }
        WwEntry *e = &cfg->ww_entries[cfg->ww_count];
        if (sscanf(item + 3, "%d %d %f", &e->kmer_size, &e->kmer_units,
                   &e->threshold) != 3) {
          fprintf(stderr,
                  "Error: %s:%d: ww requires 3 values:"
                  " <kmer_size> <kmer_units> <threshold>\n",
                  path, lineno);
          fclose(f);
          return -1;
        }
        if (e->kmer_size < 1) {
          fprintf(stderr, "Error: %s:%d: ww kmer_size must be >= 1\n", path,
                  lineno);
          fclose(f);
          return -1;
        }
        parse_optional_name(item + 3, e->name, sizeof(e->name));
        e->spacing = parse_optional_int(item + 3, "spacing");
        cfg->ww_count++;
      } else if (strncmp(item, "wh ", 3) == 0) {
        if (cfg->ww_count + cfg->wh_count >= MAX_CHANNELS) {
          fprintf(stderr, "Error: %s:%d: too many channels (max %d)\n", path,
                  lineno, MAX_CHANNELS);
          fclose(f);
          return -1;
        }
        WhEntry *e = &cfg->wh_entries[cfg->wh_count];
        if (sscanf(item + 3, "%1023s %d %f", e->path, &e->kmer_units,
                   &e->threshold) != 3) {
          fprintf(stderr,
                  "Error: %s:%d: wh requires 3 values:"
                  " <path> <kmer_units> <threshold>\n",
                  path, lineno);
          fclose(f);
          return -1;
        }
        parse_optional_name(item + 3, e->name, sizeof(e->name));
        e->spacing = parse_optional_int(item + 3, "spacing");
        cfg->wh_count++;
      } else {
        fprintf(stderr, "Error: %s:%d: unknown channel type: %s\n", path,
                lineno, item);
        fclose(f);
        return -1;
      }
      continue;
    }

    /* Key: value pair */
    char *colon = strchr(s, ':');
    if (!colon) {
      fprintf(stderr, "Error: %s:%d: expected 'key: value'\n", path, lineno);
      fclose(f);
      return -1;
    }

    *colon = '\0';
    char *key = trim(s);
    char *val = trim(colon + 1);

    if (strcmp(key, "channels") == 0) {
      in_channels = 1;
      continue;
    }

    in_channels = 0;

    if (strcmp(key, "aa") == 0) {
      if (parse_bool(val, &cfg->aa) != 0) {
        fprintf(stderr, "Error: %s:%d: invalid bool for aa: %s\n", path,
                lineno, val);
        fclose(f);
        return -1;
      }
    } else if (strcmp(key, "bed") == 0) {
      if (parse_bool(val, &cfg->bed) != 0) {
        fprintf(stderr, "Error: %s:%d: invalid bool for bed: %s\n", path,
                lineno, val);
        fclose(f);
        return -1;
      }
    } else if (strcmp(key, "raw") == 0) {
      if (parse_bool(val, &cfg->raw) != 0) {
        fprintf(stderr, "Error: %s:%d: invalid bool for raw: %s\n", path,
                lineno, val);
        fclose(f);
        return -1;
      }
    } else if (strcmp(key, "fit") == 0) {
      if (parse_bool(val, &cfg->fit) != 0) {
        fprintf(stderr, "Error: %s:%d: invalid bool for fit: %s\n", path,
                lineno, val);
        fclose(f);
        return -1;
      }
    } else if (strcmp(key, "x_mask") == 0) {
      if (parse_bool(val, &cfg->x_mask) != 0) {
        fprintf(stderr, "Error: %s:%d: invalid bool for x_mask: %s\n", path,
                lineno, val);
        fclose(f);
        return -1;
      }
    } else if (strcmp(key, "out") == 0) {
      strncpy(cfg->out, val, sizeof(cfg->out) - 1);
      cfg->out[sizeof(cfg->out) - 1] = '\0';
    } else {
      fprintf(stderr, "Error: %s:%d: unknown key: %s\n", path, lineno, key);
      fclose(f);
      return -1;
    }
  }

  fclose(f);
  return 0;
}
