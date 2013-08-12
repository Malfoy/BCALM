CC=g++
CFLAGS=-W -Wall -O3 -std=c++0x
LDFLAGS=
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
	

rebuild: clean $(EXEC)
	

