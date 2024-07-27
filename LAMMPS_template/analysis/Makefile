#	Generic Makefile - Probably requires GnuMake
#	Erik Luijten, original May 1998
#	$Id$

CC	= icc

#CFLAGS	= -Wall -ansi -pedantic

#OPT	= -g
OPT	= -O2 -s

LDFLAGS = -lm

OBJECTS  = alloc2d.o alloc3d.o matrix70.o

PROGRAMS = matrix70

default: matrix70

all: $(PROGRAMS)

%.o : %.c
	$(CC) $(CFLAGS) $(OPT) -c $<

matrix70: $(OBJECTS)
	$(CC) -o matrix70 $(CFLAGS) $(OPT) $(OBJECTS) $(LDFLAGS)

### NO CHANGES BELOW THIS LINE ###

.PHONY: cleanall clean remake

cleanall: clean
	@rm -f $(PROGRAMS)

remake: cleanall default

