# $Id$
# Makefile for GTB-developed C/C++ language software

#DEST_ROOT=/usr/local/gtb

EXES=bardcnv
OBJS=	\
		main.o params.o readin.o viterbi.o probabilities.o printout.o baumwelch.o forward.o backward.o plots.o precision.o memory.o

SRCS= $(OBJS:.o=.c)
CC = gcc
CFLAGS = -m64
#CFLAGS = -omp -g 
CFLAGS = -O3
LIB = -lm

# make all of these targets recursively through the sub-directories named above
all: $(EXES)

$(EXES): $(OBJS) $(INCS)
	$(CC) $(CFLAGS) -o $(EXES) $(OBJS) $(LIB)

clean realclean:
	rm -f $(EXES) *.o

test: all
	@echo "No tests defined (yet) for bardcnv."

