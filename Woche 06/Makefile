# Common definitions
CC = mpicc

# Compiler flags, paths and libraries
CFLAGS = -std=c11 -pedantic -Wall -Wextra -O2 -ggdb
LFLAGS = $(CFLAGS)
LIBS   = -lm

# Targets ...
all: timempi timempi2

timempi: timempi.c Makefile
	$(CC) $(LFLAGS) -o $@ timempi.c $(LIBS)

timempi2: timempi2.c Makefile
	$(CC) $(LFLAGS) -o $@ timempi2.c $(LIBS)

clean:
	$(RM) *.o *~ timempi2 timempi
