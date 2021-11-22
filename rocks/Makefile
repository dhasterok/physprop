# Makefile for rockprop
CC=gcc
CFLAGS=-g -Wall

SOURCES = physprop.c \
	rockprop.c

all: rockprop

physprop.o: physprop.c
	${CC} ${CFLAGS} -c -o $@ $^

rockprop: rockprop.c physprop.o
	${CC} ${CFLAGS} -o $@ $^ -lm

clean:
	rm rockprop *.o
