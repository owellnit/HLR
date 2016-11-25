CC     = mpicc
CFLAGS = -std=c11 -pedantic -Wall -Wextra -O2 -g -ggdb
LFLAGS = $(CFLAGS)
LIBS   = -lm

all: timempi timempi2

timempi: timempi.c Makefile
	$(CC) $(LFLAGS) -o $@ timempi.c $(LIBS)

timempi2: timempi2.c Makefile
	$(CC) $(LFLAGS) -o $@ timempi2.c $(LIBS)

clean:
	$(RM) *.o *~ timempi2 timempi
