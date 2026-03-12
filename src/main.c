#include "arg.h"
#include "config.h"
#include "fit.h"
#include "io.h"
#include "kmer_hash.h"
#include "low_complexity_models.h"
#include "metachannel.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*
 * Pre-scan argv for --ww, --wh, and --config entries, stripping them out.
 * Each --ww consumes 3 subsequent args: kmer_size kmer_units threshold
 * Each --wh consumes 3 subsequent args: path kmer_units threshold
 * --config consumes 1 subsequent arg: path
 * Returns 0 on success, -1 on error.
 */
static int prescan_channel_args(int argc, const char **argv,
                                WwEntry *ww_entries, int *ww_count,
                                WhEntry *wh_entries, int *wh_count,
                                char *config_path,
                                char *ratios_str,
                                int *new_argc, const char ***new_argv) {
  *ww_count = 0;
  *wh_count = 0;
  config_path[0] = '\0';
  ratios_str[0] = '\0';
  const char **out = (const char **)malloc(sizeof(const char *) * argc);
  int out_idx = 0;

  for (int i = 0; i < argc; i++) {
    if (strcmp(argv[i], "--ww") == 0) {
      if (i + 3 >= argc) {
        fprintf(stderr,
                "Error: --ww requires 3 arguments:"
                " <kmer_size> <kmer_units> <threshold>\n");
        free(out);
        return -1;
      }
      if (*ww_count + *wh_count >= MAX_CHANNELS) {
        fprintf(stderr, "Error: too many channels (max %d)\n", MAX_CHANNELS);
        free(out);
        return -1;
      }
      WwEntry *e = &ww_entries[*ww_count];
      e->kmer_size = atoi(argv[i + 1]);
      e->kmer_units = atoi(argv[i + 2]);
      e->threshold = (float)atof(argv[i + 3]);
      e->name[0] = '\0';
      e->spacing = 0;
      if (e->kmer_size < 1) {
        fprintf(stderr, "Error: --ww kmer_size must be >= 1\n");
        free(out);
        return -1;
      }
      (*ww_count)++;
      i += 3;
      /* Peek at next arg for spacing=N */
      if (i + 1 < argc && strncmp(argv[i + 1], "spacing=", 8) == 0) {
        e->spacing = atoi(argv[i + 1] + 8);
        i += 1;
      }
    } else if (strcmp(argv[i], "--wh") == 0) {
      if (i + 3 >= argc) {
        fprintf(stderr,
                "Error: --wh requires 3 arguments:"
                " <path> <kmer_units> <threshold>\n");
        free(out);
        return -1;
      }
      if (*ww_count + *wh_count >= MAX_CHANNELS) {
        fprintf(stderr, "Error: too many channels (max %d)\n", MAX_CHANNELS);
        free(out);
        return -1;
      }
      WhEntry *e = &wh_entries[*wh_count];
      strncpy(e->path, argv[i + 1], sizeof(e->path) - 1);
      e->path[sizeof(e->path) - 1] = '\0';
      e->kmer_units = atoi(argv[i + 2]);
      e->threshold = (float)atof(argv[i + 3]);
      e->name[0] = '\0';
      e->spacing = 0;
      (*wh_count)++;
      i += 3;
      /* Peek at next arg for spacing=N */
      if (i + 1 < argc && strncmp(argv[i + 1], "spacing=", 8) == 0) {
        e->spacing = atoi(argv[i + 1] + 8);
        i += 1;
      }
    } else if (strcmp(argv[i], "--config") == 0) {
      if (i + 1 >= argc) {
        fprintf(stderr, "Error: --config requires 1 argument: <path>\n");
        free(out);
        return -1;
      }
      strncpy(config_path, argv[i + 1], 1023);
      config_path[1023] = '\0';
      i += 1;
    } else if (strcmp(argv[i], "--ratios") == 0) {
      if (i + 1 >= argc) {
        fprintf(stderr, "Error: --ratios requires 1 argument\n");
        free(out);
        return -1;
      }
      strncpy(ratios_str, argv[i + 1], 1023);
      ratios_str[1023] = '\0';
      i += 1;
    } else {
      out[out_idx++] = argv[i];
    }
  }

  *new_argc = out_idx;
  *new_argv = out;
  return 0;
}

int main(int argc, const char **argv) {
  WwEntry ww_entries[MAX_CHANNELS];
  WhEntry wh_entries[MAX_CHANNELS];
  int ww_count = 0;
  int wh_count = 0;
  char config_path[1024] = {0};
  char ratios_str[1024] = {0};
  int filtered_argc;
  const char **filtered_argv;

  MetaChannel *meta_entries = NULL;
  int meta_count = 0;

  if (prescan_channel_args(argc, argv, ww_entries, &ww_count, wh_entries,
                           &wh_count, config_path, ratios_str,
                           &filtered_argc, &filtered_argv) != 0) {
    return 1;
  }

  /* Load config file if specified, then merge CLI channels on top */
  Config cfg;
  int have_config = 0;
  if (config_path[0] != '\0') {
    if (config_load(&cfg, config_path) != 0) {
      free((void *)filtered_argv);
      return 1;
    }
    have_config = 1;
  }

  _Bool use_aa = 0;
  _Bool use_bed = 0;
  _Bool use_x = 0;
  _Bool use_raw = 0;
  _Bool show_help = 0;
  _Bool show_version = 0;
  char file_path[1024] = {0};
  char out_path[1024] = {0};

  ArgumentList arguments = create_argument_list();

  add_argument(&arguments, (Arg){.short_name = 'h',
                                 .long_name = "help",
                                 .bool_buffer = &show_help,
                                 .help_string = "Show this help message"});

  add_argument(&arguments, (Arg){.short_name = 'v',
                                 .long_name = "version",
                                 .bool_buffer = &show_version,
                                 .help_string = "Show version information"});

  add_argument(&arguments,
               (Arg){.long_name = "aa",
                     .bool_buffer = &use_aa,
                     .help_string = "Amino acid mode (default: DNA)"});

  add_argument(&arguments,
               (Arg){.short_name = 'b',
                     .long_name = "bed",
                     .bool_buffer = &use_bed,
                     .help_string = "Output BED format (default: FASTA)"});

  add_argument(&arguments,
               (Arg){.short_name = 'x',
                     .bool_buffer = &use_x,
                     .help_string = "Hard-masking: N for DNA, X for amino acids (default: lowercase soft-masking)"});

  add_argument(&arguments, (Arg){.short_name = 'r',
                                 .long_name = "raw",
                                 .bool_buffer = &use_raw,
                                 .help_string = "Output raw entropy values"});

  add_argument(&arguments,
               (Arg){.short_name = 'o',
                     .long_name = "out",
                     .string_buffer = out_path,
                     .help_string = "Output file path (default: stdout)"});

  add_argument(&arguments, (Arg){.short_name = POSITIONAL_ARG,
                                 .string_buffer = file_path,
                                 .help_string = "Input FASTA file"});

  process_args(filtered_argc, filtered_argv, &arguments);

  if (show_version) {
    if (strcmp(GIT_COMMIT, "unknown") == 0)
      fprintf(stdout, "loco %s\n", VERSION);
    else
      fprintf(stdout, "loco %s (%s)\n", VERSION, GIT_COMMIT);
    free((void *)filtered_argv);
    return 0;
  }

  if (show_help) {
    print_help(argv[0], &arguments);
    fprintf(stderr,
            "\nChannel options (at least one required):\n"
            "  --ww <kmer> <kmer_units> <threshold> [spacing=N]\n"
            "                                     "
            "Kmer entropy channel (repeatable)\n"
            "  --wh <path> <kmer_units> <threshold> [spacing=N]\n"
            "                                     "
            "Embedding entropy channel (repeatable)\n"
            "  --config <file>                    "
            "Load settings and channels from YAML config\n"
            "                                     "
            "(use 'fit: true' in config to calibrate thresholds)\n"
            "  --ratios <r1,r2,...>               "
            "Letter frequencies for fit mode\n"
            "                                     "
            "DNA order: A,C,G,T\n"
            "                                     "
            "AA order:  A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y\n"
            "                                     "
            "(default: 60%% AT for DNA, Swiss-Prot bg for AA)\n");
    free((void *)filtered_argv);
    return 0;
  }

  /*
   * Merge config file settings with CLI.
   * CLI flags override config values (only when actually provided).
   * CLI channels are appended to config channels.
   */
  if (have_config) {
    /* Find arg entries to check usage_count */
    Arg *arg_aa = NULL, *arg_bed = NULL, *arg_x = NULL;
    Arg *arg_raw = NULL, *arg_out = NULL;
    for (int i = 0; i < arguments.length; i++) {
      Arg *a = &arguments.arguments[i];
      if (a->long_name && strcmp(a->long_name, "aa") == 0)
        arg_aa = a;
      else if (a->long_name && strcmp(a->long_name, "bed") == 0)
        arg_bed = a;
      else if (a->short_name == 'x' && !a->long_name)
        arg_x = a;
      else if (a->long_name && strcmp(a->long_name, "raw") == 0)
        arg_raw = a;
      else if (a->long_name && strcmp(a->long_name, "out") == 0)
        arg_out = a;
    }

    /* Apply config booleans only if CLI did not provide them */
    if (arg_aa && arg_aa->usage_count == 0 && cfg.aa >= 0)
      use_aa = (_Bool)cfg.aa;
    if (arg_bed && arg_bed->usage_count == 0 && cfg.bed >= 0)
      use_bed = (_Bool)cfg.bed;
    if (arg_x && arg_x->usage_count == 0 && cfg.x_mask >= 0)
      use_x = (_Bool)cfg.x_mask;
    if (arg_raw && arg_raw->usage_count == 0 && cfg.raw >= 0)
      use_raw = (_Bool)cfg.raw;
    if (arg_out && arg_out->usage_count == 0 && cfg.out[0] != '\0')
      strncpy(out_path, cfg.out, sizeof(out_path) - 1);

    /* Prepend config channels before CLI channels */
    int total_ww = cfg.ww_count + ww_count;
    int total_wh = cfg.wh_count + wh_count;
    if (total_ww + total_wh > MAX_CHANNELS) {
      fprintf(stderr, "Error: too many channels (max %d)\n", MAX_CHANNELS);
      free((void *)filtered_argv);
      return 1;
    }

    /* Shift CLI ww entries to make room for config entries */
    if (cfg.ww_count > 0) {
      memmove(&ww_entries[cfg.ww_count], &ww_entries[0],
              sizeof(WwEntry) * ww_count);
      memcpy(&ww_entries[0], cfg.ww_entries, sizeof(WwEntry) * cfg.ww_count);
      ww_count = total_ww;
    }

    /* Shift CLI wh entries to make room for config entries */
    if (cfg.wh_count > 0) {
      memmove(&wh_entries[cfg.wh_count], &wh_entries[0],
              sizeof(WhEntry) * wh_count);
      memcpy(&wh_entries[0], cfg.wh_entries, sizeof(WhEntry) * cfg.wh_count);
      wh_count = total_wh;
    }

    /* Pick up metachannels from config */
    meta_entries = cfg.meta_entries;
    meta_count = cfg.meta_count;
  }

  /* Validate metachannels and named channels */
  {
    int has_named = 0;
    for (int i = 0; i < ww_count && !has_named; i++)
      has_named = ww_entries[i].name[0] != '\0';
    for (int i = 0; i < wh_count && !has_named; i++)
      has_named = wh_entries[i].name[0] != '\0';
    if (meta_count > 0 || has_named) {
      if (validate_metachannels(meta_entries, meta_count,
                                ww_entries, ww_count,
                                wh_entries, wh_count) != 0) {
        free((void *)filtered_argv);
        return 1;
      }
      if (meta_count > 1) {
        sort_metachannels(meta_entries, meta_count,
                          ww_entries, ww_count,
                          wh_entries, wh_count);
      }
    }
  }

  /* fit mode: triggered by 'fit: true' in config */
  if (have_config && cfg.fit == 1) {
    if (ww_count == 0 && wh_count == 0) {
      fprintf(stderr,
              "Error: fit mode requires at least one ww or wh channel\n");
      free((void *)filtered_argv);
      return 1;
    }

    FILE *fit_out = stdout;
    if (out_path[0] != '\0') {
      fit_out = fopen(out_path, "w");
      if (!fit_out) {
        fprintf(stderr, "Error: could not open output file: %s\n", out_path);
        free((void *)filtered_argv);
        return 1;
      }
    }

    const char *fit_ratios = ratios_str[0] ? ratios_str : NULL;
    const char *fit_fasta = file_path[0] ? file_path : NULL;

    if (fit_ratios && fit_fasta) {
      fprintf(stderr, "Error: --ratios and input file are mutually "
                      "exclusive in --fit mode\n");
      if (fit_out != stdout)
        fclose(fit_out);
      free((void *)filtered_argv);
      return 1;
    }

    int rc = fit_channels(fit_out, use_aa, ww_entries, ww_count,
                          wh_entries, wh_count, fit_ratios, fit_fasta,
                          meta_entries, meta_count);

    if (fit_out != stdout)
      fclose(fit_out);
    free((void *)filtered_argv);
    return rc == 0 ? 0 : 1;
  }

  if (file_path[0] == '\0') {
    fprintf(stderr, "Error: no input file specified\n");
    print_help(argv[0], &arguments);
    free((void *)filtered_argv);
    return 1;
  }

  if (ww_count == 0 && wh_count == 0) {
    fprintf(stderr,
            "Error: at least one --ww or --wh channel is required\n");
    free((void *)filtered_argv);
    return 1;
  }

  FILE *out_file = stdout;
  if (out_path[0] != '\0') {
    out_file = fopen(out_path, "w");
    if (!out_file) {
      fprintf(stderr, "Error: could not open output file: %s\n", out_path);
      free((void *)filtered_argv);
      return 1;
    }
  }

  init_lookup_tables();

  /* Build lookup tables for --ww entries */
  Entropy **ww_tables = NULL;
  if (ww_count > 0) {
    ww_tables = (Entropy **)malloc(sizeof(Entropy *) * ww_count);
    for (int w = 0; w < ww_count; w++) {
      int num_kmer_samples;
      if (ww_entries[w].spacing > 0) {
        num_kmer_samples = ww_entries[w].kmer_units;
      } else {
        int ws = ww_entries[w].kmer_units * ww_entries[w].kmer_size;
        num_kmer_samples = ws - (ww_entries[w].kmer_size - 1);
      }
      ww_tables[w] = (Entropy *)malloc(sizeof(Entropy) * (num_kmer_samples + 1));
      fill_count_to_entropy_lookup_table_kmer(ww_tables[w],
                                              (uint64_t)num_kmer_samples);
    }
  }

  /* Load --wh bin files and build bio_to_alpha mappings */
  if (wh_count > 0) {
    BioLetter *lookup = use_aa ? char_to_amino : char_to_dna;
    for (int w = 0; w < wh_count; w++) {
      if (kmer_hash_load(&wh_entries[w].hash, wh_entries[w].path) != 0) {
        fprintf(stderr, "Error: could not load bin file: %s\n",
                wh_entries[w].path);
        free((void *)filtered_argv);
        return 1;
      }
      for (int i = 0; i < 21; i++)
        wh_entries[w].bio_to_alpha[i] = 0;
      const char *alpha = wh_entries[w].hash.alphabet;
      for (uint32_t i = 0; i < wh_entries[w].hash.alphabet_size; i++) {
        BioLetter bl = lookup[(unsigned char)alpha[i]];
        wh_entries[w].bio_to_alpha[bl] = (int)i;
      }
    }
  }

  int total_channels = ww_count + wh_count;

  BioSequence *sequences =
      use_aa ? read_fasta_aa(file_path) : read_fasta_dna(file_path);

  BioSequence *seq = sequences;
  while (seq) {
    uint64_t len = seq->sequence_length;

    Entropy **channels =
        (Entropy **)malloc(sizeof(Entropy *) * total_channels);
    float *thresholds = (float *)malloc(sizeof(float) * total_channels);
    int ch = 0;

    /* Compute --ww channels */
    for (int w = 0; w < ww_count; w++) {
      channels[ch] = (Entropy *)malloc(sizeof(Entropy) * len);
      thresholds[ch] = ww_entries[w].threshold;
      int ws;
      if (ww_entries[w].spacing > 0) {
        int num_states = ww_entries[w].spacing + 1;
        ws = (ww_entries[w].kmer_units - 1) * num_states +
             ww_entries[w].kmer_size;
      } else {
        ws = ww_entries[w].kmer_units * ww_entries[w].kmer_size;
      }
      fill_complexity_buffer_for_seq_kmer(
          seq->sequence, len, channels[ch],
          (uint64_t)ws,
          (uint32_t)ww_entries[w].kmer_size, ww_tables[w],
          (uint32_t)ww_entries[w].spacing);
      ch++;
    }

    /* Compute --wh channels */
    for (int w = 0; w < wh_count; w++) {
      channels[ch] = (Entropy *)malloc(sizeof(Entropy) * len);
      thresholds[ch] = wh_entries[w].threshold;
      int ws;
      if (wh_entries[w].spacing > 0) {
        int num_states = wh_entries[w].spacing + 1;
        ws = (wh_entries[w].kmer_units - 1) * num_states +
             (int)wh_entries[w].hash.kmer_size;
      } else {
        ws = wh_entries[w].kmer_units * (int)wh_entries[w].hash.kmer_size;
      }
      fill_complexity_buffer_for_seq_hash(
          seq->sequence, len, channels[ch],
          (uint64_t)ws, &wh_entries[w].hash,
          wh_entries[w].bio_to_alpha,
          (uint32_t)wh_entries[w].spacing);
      ch++;
    }

    if (use_raw) {
      output_raw_entropies(out_file, seq, channels, total_channels);
    } else {
      _Bool *mask = (_Bool *)calloc(len, sizeof(_Bool));

      /* OR unnamed channels (backward compatible) */
      for (int c = 0; c < ww_count; c++) {
        if (ww_entries[c].name[0] != '\0')
          continue;
        for (uint64_t i = 0; i < len; i++)
          mask[i] |= channels[c][i] < thresholds[c];
      }
      for (int c = 0; c < wh_count; c++) {
        if (wh_entries[c].name[0] != '\0')
          continue;
        int ci = ww_count + c;
        for (uint64_t i = 0; i < len; i++)
          mask[i] |= channels[ci][i] < thresholds[ci];
      }

      /* Bulk-compute metachannel masks (topo-sorted, each computed once) */
      if (meta_count > 0) {
        _Bool **meta_masks = compute_meta_masks(
            meta_entries, meta_count,
            ww_entries, ww_count, wh_entries, wh_count,
            channels, thresholds, len);
        if (meta_masks) {
          for (int m = 0; m < meta_count; m++) {
            for (uint64_t i = 0; i < len; i++)
              mask[i] |= meta_masks[m][i];
            free(meta_masks[m]);
          }
          free(meta_masks);
        }
      }

      if (use_bed) {
        output_bed_from_mask(out_file, seq, mask);
      } else {
        output_fasta_masked(out_file, seq, mask, use_x, use_aa);
      }
      free(mask);
    }

    for (int c = 0; c < total_channels; c++)
      free(channels[c]);
    free(channels);
    free(thresholds);

    seq = (BioSequence *)seq->next_sequence;
  }

  /* Cleanup */
  if (ww_tables) {
    for (int w = 0; w < ww_count; w++)
      free(ww_tables[w]);
    free(ww_tables);
  }
  for (int w = 0; w < wh_count; w++)
    kmer_hash_free(&wh_entries[w].hash);

  if (out_file != stdout)
    fclose(out_file);

  free((void *)filtered_argv);
  return 0;
}
