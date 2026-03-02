CC      = gcc
CFLAGS  = -std=c17 -Wall -Wextra -pedantic
TARGET  = loco
VERSION = 0.1.0

SRCDIR  = src
BUILDDIR = build
SRCS    = $(wildcard $(SRCDIR)/*.c)
OBJS    = $(patsubst $(SRCDIR)/%.c,$(BUILDDIR)/%.o,$(SRCS))

GIT_COMMIT := $(shell git rev-parse --short HEAD 2>/dev/null || echo "unknown")
CFLAGS += -DVERSION=\"$(VERSION)\" -DGIT_COMMIT=\"$(GIT_COMMIT)\"

.PHONY: all debug clean test

all: CFLAGS += -O3
all: $(TARGET)

debug: CFLAGS += -O0 -g -DDEBUG
debug: $(TARGET)

$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) -o $@ $^ -lm

$(BUILDDIR)/%.o: $(SRCDIR)/%.c | $(BUILDDIR)
	$(CC) $(CFLAGS) -c -o $@ $<

$(BUILDDIR):
	mkdir -p $(BUILDDIR)

clean:
	rm -rf $(BUILDDIR) $(TARGET)

test: CFLAGS += -O0 -g -DDEBUG
test: $(TARGET)
	./$(TARGET)
