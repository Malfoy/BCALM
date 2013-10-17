CC=g++ 
CFLAGS=-W -Wall -O4 -std=c++0x -march=native
LDFLAGS=


ifeq ($(gprof),1)
CFLAGS=-std=c++0x -pg -march=native -O4
LDFLAGS=-pg 
endif

ifeq ($(valgrind),1)
CFLAGS=-std=c++0x -O4 -g
LDFLAGS=-g
endif



EXEC=bcalm

all: $(EXEC)

bcalm: main.o lm.o ograph.o debug.o
	$(CC) -o $@ $^ $(LDFLAGS)

debug.o: debug.cpp ograph.h
	$(CC) -o $@ -c $< $(CFLAGS)

main.o: main.cpp lm.h ograph.h debug.h
	$(CC) -o $@ -c $< $(CFLAGS)

ograph.o: ograph.cpp
	$(CC) -o $@ -c $< $(CFLAGS)

lm.o: lm.cpp ograph.h
	$(CC) -o $@ -c $< $(CFLAGS)

script: scheck.o ograph.o
	$(CC) -o $@ $^ $(LDFLAGS)

scheck.o:scheck.cpp ograph.h
	$(CC) -o $@ -c $< $(CFLAGS)

clean:
	rm -rf *.o
	rm -rf $(EXEC)
	rm -rf *.dot
	rm -rf randomgenome


rebuild: clean $(EXEC)


