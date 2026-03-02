#include "io.h"
#include <fcntl.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>

static void append_node(BioSequence **head, BioSequence **tail, char *name,
                        BioLetter *seq, char *raw, uint64_t len,
                        char *block, uint64_t block_len) {
  BioSequence *node = (BioSequence *)malloc(sizeof(BioSequence));
  node->sequence_name = name;
  node->sequence = seq;
  node->raw_sequence = raw;
  node->raw_block = (char *)malloc(block_len);
  memcpy(node->raw_block, block, block_len);
  node->sequence_length = len;
  node->raw_block_length = block_len;
  node->next_sequence = NULL;
  if (*tail) {
    (*tail)->next_sequence = (struct BioSequence *)node;
  } else {
    *head = node;
  }
  *tail = node;
}

static BioSequence *read_fasta(char *file_path, BioLetter *lookup) {
  int fd = open(file_path, O_RDONLY);
  if (fd < 0) {
    fprintf(stderr, "Could not open file: %s\n", file_path);
    return NULL;
  }

  struct stat st;
  fstat(fd, &st);
  size_t file_size = (size_t)st.st_size;
  if (file_size == 0) {
    close(fd);
    return NULL;
  }

  const char *data =
      (const char *)mmap(NULL, file_size, PROT_READ, MAP_PRIVATE, fd, 0);
  close(fd);
  if (data == MAP_FAILED) {
    fprintf(stderr, "Could not mmap file: %s\n", file_path);
    return NULL;
  }

  // file_size is an upper bound on total sequence characters
  BioLetter *seq_pool = (BioLetter *)malloc(sizeof(BioLetter) * file_size);
  char *raw_pool = (char *)malloc(file_size);

  BioSequence *head = NULL;
  BioSequence *tail = NULL;

  BioLetter *seq_buf = NULL;
  char *raw_buf = NULL;
  uint64_t seq_len = 0;
  uint64_t pool_offset = 0;
  char *name = NULL;
  const char *block_start = NULL;

  const char *p = data;
  const char *end = data + file_size;

  while (p < end) {
    if (*p == '>') {
      // Save previous sequence
      if (name) {
        uint64_t block_len = (uint64_t)(p - block_start);
        // Trim trailing newlines from block
        while (block_len > 0 && (block_start[block_len - 1] == '\n' ||
                                  block_start[block_len - 1] == '\r')) {
          block_len--;
        }
        append_node(&head, &tail, name, seq_buf, raw_buf, seq_len,
                    (char *)block_start, block_len);
        pool_offset += seq_len;
      }

      // Parse header line
      p++;
      const char *line_start = p;
      while (p < end && *p != '\n' && *p != '\r') {
        p++;
      }
      size_t name_len = (size_t)(p - line_start);
      name = (char *)malloc(name_len + 1);
      memcpy(name, line_start, name_len);
      name[name_len] = '\0';

      seq_len = 0;
      seq_buf = seq_pool + pool_offset;
      raw_buf = raw_pool + pool_offset;

      // Skip newline
      if (p < end && *p == '\r')
        p++;
      if (p < end && *p == '\n')
        p++;
      block_start = p;
    } else if (*p == '\n' || *p == '\r') {
      p++;
    } else if (name) {
      // Sequence data - scan to end of line
      while (p < end && *p != '\n' && *p != '\r') {
        raw_buf[seq_len] = *p;
        seq_buf[seq_len] = lookup[(unsigned char)*p];
        seq_len++;
        p++;
      }
    } else {
      p++;
    }
  }

  // Save last sequence
  if (name) {
    uint64_t block_len = (uint64_t)(end - block_start);
    while (block_len > 0 && (block_start[block_len - 1] == '\n' ||
                              block_start[block_len - 1] == '\r')) {
      block_len--;
    }
    append_node(&head, &tail, name, seq_buf, raw_buf, seq_len,
                (char *)block_start, block_len);
  }

  munmap((void *)data, file_size);
  return head;
}

void output_BED_for_entropies(FILE *out_file, BioSequence *sequence,
                              Entropy *complexity_entropies, Entropy cutoff) {
  uint64_t region_start = 0;
  int in_region = 0;

  for (uint64_t i = 0; i < sequence->sequence_length; i++) {
    if (complexity_entropies[i] > cutoff) {
      if (!in_region) {
        region_start = i;
        in_region = 1;
      }
    } else {
      if (in_region) {
        fprintf(out_file, "%s\t%llu\t%llu\n", sequence->sequence_name,
                (unsigned long long)region_start, (unsigned long long)i);
        in_region = 0;
      }
    }
  }

  if (in_region) {
    fprintf(out_file, "%s\t%llu\t%llu\n", sequence->sequence_name,
            (unsigned long long)region_start,
            (unsigned long long)sequence->sequence_length);
  }
}

void output_entropies(FILE *out_file, BioSequence *sequence,
                      Entropy *complexity_entropies) {
  fprintf(out_file, ">%s\n", sequence->sequence_name);
  for (uint64_t i = 0; i < sequence->sequence_length; i++) {
    fprintf(out_file, "%llu: %f\n", i, complexity_entropies[i]);
  }
}

void output_bed_from_mask(FILE *out_file, BioSequence *sequence,
                          _Bool *mask) {
  /* Extract sequence identifier (first whitespace-delimited token) */
  const char *name = sequence->sequence_name;
  int name_len = 0;
  while (name[name_len] && name[name_len] != ' ' && name[name_len] != '\t')
    name_len++;

  uint64_t region_start = 0;
  int in_region = 0;

  for (uint64_t i = 0; i < sequence->sequence_length; i++) {
    if (mask[i]) {
      if (!in_region) {
        region_start = i;
        in_region = 1;
      }
    } else {
      if (in_region) {
        fprintf(out_file, "%.*s\t%llu\t%llu\n", name_len, name,
                (unsigned long long)region_start, (unsigned long long)i);
        in_region = 0;
      }
    }
  }

  if (in_region) {
    fprintf(out_file, "%.*s\t%llu\t%llu\n", name_len, name,
            (unsigned long long)region_start,
            (unsigned long long)sequence->sequence_length);
  }
}

void output_fasta_masked(FILE *out_file, BioSequence *sequence, _Bool *mask,
                         _Bool use_x, _Bool use_aa) {
  fprintf(out_file, ">%s\n", sequence->sequence_name);
  uint64_t seq_idx = 0;
  for (uint64_t i = 0; i < sequence->raw_block_length; i++) {
    char c = sequence->raw_block[i];
    if (c == '\n' || c == '\r') {
      fputc(c, out_file);
    } else {
      if (seq_idx < sequence->sequence_length && mask[seq_idx]) {
        if (use_x) {
          c = use_aa ? 'X' : 'N';
        } else {
          if (c >= 'A' && c <= 'Z') {
            c = c + 32;
          }
        }
      }
      fputc(c, out_file);
      seq_idx++;
    }
  }
  fputc('\n', out_file);
}

void output_raw_entropies(FILE *out_file, BioSequence *sequence,
                          Entropy **ent_channels, int num_channels) {
  fprintf(out_file, ">%s\n", sequence->sequence_name);
  for (uint64_t i = 0; i < sequence->sequence_length; i++) {
    for (int c = 0; c < num_channels; c++) {
      if (c > 0)
        fputc(' ', out_file);
      fprintf(out_file, "%f", ent_channels[c][i]);
    }
    fputc('\n', out_file);
  }
}

BioSequence *read_fasta_aa(char *file_path) {
  return read_fasta(file_path, char_to_amino);
}

BioSequence *read_fasta_dna(char *file_path) {
  return read_fasta(file_path, char_to_dna);
}
