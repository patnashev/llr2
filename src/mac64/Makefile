# Makefile for macOS llr64

CC = gcc

LFLAGS =
LIBELF = ../gwnum/mac64/gwnum.a
LIBS2  = $(LIBELF) -lm -lpthread

COMPOBJS = lprime.o
EXE2 = llr2

#########################################################################

all: $(EXE2)

$(EXE2): $(COMPOBJS) $(LIBELF)
	$(CC) $(LFLAGS) -o $(EXE2) $(COMPOBJS) $(LIBS2)

$(COMPOBJS):
	make -f make64

clean:
	rm -f $(EXE2) $(COMPOBJS)

