#ifndef ARG_H
#define ARG_H

#define POSITIONAL_ARG '\0'

typedef struct {
  char short_name;
  char *long_name; // optional long name
  char *help_string;

  void (*callback)(char *);

  // Only one of these should be set
  _Bool *bool_buffer; // Note, positional arguments cannot be bools
  int *int_buffer;
  float *float_buffer;
  char *string_buffer;

  char *default_value;

  int usage_count;
} Arg;

typedef struct {
  Arg *arguments;
  int length;
  int capacity;
} ArgumentList;

ArgumentList create_argument_list(void);
void add_argument(ArgumentList *arg_list, Arg argument);
void process_args(int argc, const char **argv, ArgumentList *arg_list);
void print_help(const char *program_name, ArgumentList *arg_list);

#endif
