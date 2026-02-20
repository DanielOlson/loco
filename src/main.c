#include "arg.h"
#include "io.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main(int argc, const char **argv) {
  int w1, w2, w3, w4;
  _Bool use_aa = 0;
  _Bool use_bed = 0;
  _Bool use_x = 0;
  _Bool use_raw = 0;
  _Bool show_help = 0;
  float thresh_1, thresh_2, thresh_3, thresh_4;
  char file_path[1024] = {0};
  char out_path[1024] = {0};

  ArgumentList arguments = create_argument_list();

  add_argument(&arguments, (Arg){.short_name = 'h',
                                 .long_name = "help",
                                 .bool_buffer = &show_help,
                                 .help_string = "Show this help message"});

  add_argument(&arguments, (Arg){.long_name = "w1",
                                 .int_buffer = &w1,
                                 .default_value = "20",
                                 .help_string = "1-mer window size"});

  add_argument(&arguments, (Arg){.long_name = "w2",
                                 .int_buffer = &w2,
                                 .default_value = "30",
                                 .help_string = "2-mer window size"});

  add_argument(&arguments, (Arg){.long_name = "w3",
                                 .int_buffer = &w3,
                                 .default_value = "40",
                                 .help_string = "3-mer window size"});

  add_argument(&arguments, (Arg){.long_name = "w4",
                                 .int_buffer = &w4,
                                 .default_value = "50",
                                 .help_string = "4-mer window size"});

  add_argument(&arguments,
               (Arg){.long_name = "aa",
                     .bool_buffer = &use_aa,
                     .help_string = "Amino acid mode (default: DNA)"});

  add_argument(&arguments, (Arg){.short_name = '1',
                                 .float_buffer = &thresh_1,
                                 .default_value = "1.00",
                                 .help_string = "1-mer entropy threshold"});

  add_argument(&arguments, (Arg){.short_name = '2',
                                 .float_buffer = &thresh_2,
                                 .default_value = "2.05",
                                 .help_string = "2-mer entropy threshold"});

  add_argument(&arguments, (Arg){.short_name = '3',
                                 .float_buffer = &thresh_3,
                                 .default_value = "2.90",
                                 .help_string = "3-mer entropy threshold"});

  add_argument(&arguments, (Arg){.short_name = '4',
                                 .float_buffer = &thresh_4,
                                 .default_value = "3.47",
                                 .help_string = "4-mer entropy threshold"});

  add_argument(&arguments,
               (Arg){.short_name = 'b',
                     .long_name = "bed",
                     .bool_buffer = &use_bed,
                     .help_string = "Output BED format (default: FASTA)"});

  add_argument(&arguments,
               (Arg){.short_name = 'x',
                     .bool_buffer = &use_x,
                     .help_string = "X-masking (default: lowercase)"});

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

  process_args(argc, argv, &arguments);

  if (show_help) {
    print_help(argv[0], &arguments);
    return 0;
  }

  if (file_path[0] == '\0') {
    fprintf(stderr, "Error: no input file specified\n");
    print_help(argv[0], &arguments);
    return 1;
  }

  FILE *out_file = stdout;
  if (out_path[0] != '\0') {
    out_file = fopen(out_path, "w");
    if (!out_file) {
      fprintf(stderr, "Error: could not open output file: %s\n", out_path);
      return 1;
    }
  }

  init_lookup_tables();

  // Build entropy lookup tables for all 4 k-mer sizes
  Entropy *ent_table_1 = (Entropy *)malloc(sizeof(Entropy) * (w1 + 1));
  Entropy *ent_table_2 = (Entropy *)malloc(sizeof(Entropy) * (w2 + 1));
  Entropy *ent_table_3 = (Entropy *)malloc(sizeof(Entropy) * (w3 + 1));
  Entropy *ent_table_4 = (Entropy *)malloc(sizeof(Entropy) * (w4 + 1));
  fill_count_to_entropy_lookup_table(ent_table_1, w1);
  fill_count_to_entropy_lookup_table_2mer(ent_table_2, w2);
  fill_count_to_entropy_lookup_table_3mer(ent_table_3, w3);
  fill_count_to_entropy_lookup_table_4mer(ent_table_4, w4);

  BioSequence *sequences =
      use_aa ? read_fasta_aa(file_path) : read_fasta_dna(file_path);

  BioSequence *seq = sequences;
  while (seq) {
    uint64_t len = seq->sequence_length;
    Entropy *ent_1 = (Entropy *)malloc(sizeof(Entropy) * len);
    Entropy *ent_2 = (Entropy *)malloc(sizeof(Entropy) * len);
    Entropy *ent_3 = (Entropy *)malloc(sizeof(Entropy) * len);
    Entropy *ent_4 = (Entropy *)malloc(sizeof(Entropy) * len);

    fill_complexity_buffer_for_seq(seq->sequence, len, ent_1, w1, ent_table_1);
    fill_complexity_buffer_for_seq_2mer(seq->sequence, len, ent_2, w2,
                                        ent_table_2);
    fill_complexity_buffer_for_seq_3mer(seq->sequence, len, ent_3, w3,
                                        ent_table_3);
    fill_complexity_buffer_for_seq_4mer(seq->sequence, len, ent_4, w4,
                                        ent_table_4);

    if (use_raw) {
      output_raw_entropies(out_file, seq, ent_1, ent_2, ent_3, ent_4);
    } else {
      _Bool *mask = (_Bool *)calloc(len, sizeof(_Bool));
      for (uint64_t i = 0; i < len; i++) {
        if (ent_1[i] < thresh_1 || ent_2[i] < thresh_2 || ent_3[i] < thresh_3 ||
            ent_4[i] < thresh_4) {
          mask[i] = 1;
        }
      }

      if (use_bed) {
        output_bed_from_mask(out_file, seq, mask);
      } else {
        output_fasta_masked(out_file, seq, mask, use_x);
      }
      free(mask);
    }

    free(ent_1);
    free(ent_2);
    free(ent_3);
    free(ent_4);
    seq = (BioSequence *)seq->next_sequence;
  }

  free(ent_table_1);
  free(ent_table_2);
  free(ent_table_3);
  free(ent_table_4);

  if (out_file != stdout) {
    fclose(out_file);
  }

  return 0;
}
