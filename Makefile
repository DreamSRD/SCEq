HEADER = -I/usr/local/include
LIBB = -L/usr/local/lib
LIBRA = -lm -lfftw3
SOURCES = SCEq.c
CFLAGS =

all:
	gcc -std=c99 $(CFLAGS) $(SOURCES) $(HEADER) $(LIBB) $(LIBRA) -o SCEq
clean:
	rm -rf *.o SCEq
