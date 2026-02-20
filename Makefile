CC      = gcc
CFLAGS  = -std=c17 -Wall -Wextra -pedantic
TARGET  = loco

SRCDIR  = src
BUILDDIR = build
SRCS    = $(wildcard $(SRCDIR)/*.c)
OBJS    = $(patsubst $(SRCDIR)/%.c,$(BUILDDIR)/%.o,$(SRCS))

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
