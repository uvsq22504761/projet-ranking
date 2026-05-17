CC = gcc
CFLAGS = -O2
LDFLAGS = -lm

SRC = pagerank_projet.c
SRC_O = pagerank_original.c
TARGET = build/pagerank
TARGET_O = build/pagerank_o

.PHONY: all clean

all: build $(TARGET) $(TARGET_O)

build:
	mkdir -p build

$(TARGET): $(SRC)
	$(CC) $(CFLAGS) -o $@ $< $(LDFLAGS)
$(TARGET_O): $(SRC_O)
	$(CC) $(CFLAGS) -o $@ $< $(LDFLAGS)

clean:
	rm -rf build