CC = gcc
CFLAGS = -O2
LDFLAGS = -lm

SRC = pagerank_projet.c
TARGET = build/pagerank

.PHONY: all clean

all: build $(TARGET)

build:
	mkdir -p build

$(TARGET): $(SRC)
	$(CC) $(CFLAGS) -o $@ $< $(LDFLAGS)

clean:
	rm -rf build