#include "arg.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

ArgumentList create_argument_list(void) {
  ArgumentList new_list;
  new_list.arguments = (Arg *)malloc(sizeof(Arg) * 512);
  new_list.capacity = 512;
  new_list.length = 0;
  return new_list;
}

void add_argument(ArgumentList *arg_list, Arg argument) {
  if (arg_list->capacity == arg_list->length) {
    arg_list->arguments =
        realloc(arg_list->arguments, sizeof(Arg) * arg_list->capacity * 2);
    arg_list->capacity = arg_list->capacity * 2;
  }

  arg_list->arguments[arg_list->length] = argument;
  arg_list->length++;
}

static void apply_value(Arg *arg, const char *value) {
  arg->usage_count++;
  if (arg->bool_buffer) {
    *arg->bool_buffer = 1;
  } else if (arg->int_buffer) {
    *arg->int_buffer = atoi(value);
  } else if (arg->float_buffer) {
    *arg->float_buffer = (float)atof(value);
  } else if (arg->string_buffer) {
    strcpy(arg->string_buffer, value);
  }
  if (arg->callback) {
    arg->callback((char *)value);
  }
}

static Arg *find_short(ArgumentList *arg_list, char name) {
  for (int j = 0; j < arg_list->length; j++) {
    if (arg_list->arguments[j].short_name == name) {
      return &arg_list->arguments[j];
    }
  }
  return NULL;
}

static Arg *find_long(ArgumentList *arg_list, const char *name) {
  for (int j = 0; j < arg_list->length; j++) {
    if (arg_list->arguments[j].long_name &&
        strcmp(arg_list->arguments[j].long_name, name) == 0) {
      return &arg_list->arguments[j];
    }
  }
  return NULL;
}

void process_args(int argc, const char **argv, ArgumentList *arg_list) {
  // Apply defaults first
  for (int j = 0; j < arg_list->length; j++) {
    Arg *arg = &arg_list->arguments[j];
    if (arg->default_value) {
      apply_value(arg, arg->default_value);
      arg->usage_count = 0;
    }
  }

  int positional = 0;
  for (int i = 1; i < argc; ++i) {
    if (argv[i][0] == '-' && argv[i][1] == '-') {
      // Long argument: --name or --name=value
      const char *name = argv[i] + 2;
      const char *eq = strchr(name, '=');
      Arg *arg;
      if (eq) {
        char buf[256];
        int len = (int)(eq - name);
        memcpy(buf, name, len);
        buf[len] = '\0';
        arg = find_long(arg_list, buf);
        if (!arg) {
          fprintf(stderr, "Unknown argument: --%s\n", buf);
          exit(1);
        }
        apply_value(arg, eq + 1);
      } else {
        arg = find_long(arg_list, name);
        if (!arg) {
          fprintf(stderr, "Unknown argument: --%s\n", name);
          exit(1);
        }
        if (arg->bool_buffer) {
          apply_value(arg, NULL);
        } else if (i + 1 < argc) {
          apply_value(arg, argv[++i]);
        } else {
          fprintf(stderr, "Missing value for --%s\n", name);
          exit(1);
        }
      }
    } else if (argv[i][0] == '-') {
      // Short argument: -w value or -w=value
      char name = argv[i][1];
      Arg *arg = find_short(arg_list, name);
      if (!arg) {
        fprintf(stderr, "Unknown argument: -%c\n", name);
        exit(1);
      }
      if (argv[i][2] == '=') {
        apply_value(arg, argv[i] + 3);
      } else if (arg->bool_buffer) {
        apply_value(arg, NULL);
      } else if (i + 1 < argc) {
        apply_value(arg, argv[++i]);
      } else {
        fprintf(stderr, "Missing value for -%c\n", name);
        exit(1);
      }
    } else {
      // Positional argument
      Arg *arg = NULL;
      for (int j = 0, p = 0; j < arg_list->length; j++) {
        if (arg_list->arguments[j].short_name == POSITIONAL_ARG &&
            !arg_list->arguments[j].long_name) {
          if (p == positional) {
            arg = &arg_list->arguments[j];
            break;
          }
          p++;
        }
      }
      if (!arg) {
        fprintf(stderr, "Unexpected positional argument: %s\n", argv[i]);
        exit(1);
      }
      apply_value(arg, argv[i]);
      positional++;
    }
  }
}

void print_help(const char *program_name, ArgumentList *arg_list) {
  fprintf(stderr, "Usage: %s [options] <input.fa>\n\nOptions:\n", program_name);
  for (int i = 0; i < arg_list->length; i++) {
    Arg *arg = &arg_list->arguments[i];
    if (arg->short_name == POSITIONAL_ARG && !arg->long_name) {
      continue;
    }
    fprintf(stderr, "  ");
    if (arg->short_name) {
      fprintf(stderr, "-%c", arg->short_name);
      if (arg->long_name) {
        fprintf(stderr, ", ");
      }
    }
    if (arg->long_name) {
      fprintf(stderr, "--%s", arg->long_name);
    }
    if (!arg->bool_buffer) {
      fprintf(stderr, " <value>");
    }
    if (arg->help_string) {
      fprintf(stderr, "\t%s", arg->help_string);
    }
    if (arg->default_value) {
      fprintf(stderr, " (default: %s)", arg->default_value);
    }
    fprintf(stderr, "\n");
  }
}
