# Makefile for Linux llr64 or sllr64

CC = gcc

LFLAGS =
LIBELF = ../gwnum/linux64/gwnum.a
LIBS2  = $(LIBELF) ../gwnum/linux64/gwnum.ld -lm  -lpthread -lstdc++ -static

COMPOBJS = lprime.o
EXE2 = llr2

# insert proper dummy file to get _GWDATA data segment aligned
# on a 32 byte boundary.  Use proper dummyNN.o to do the alignment

DUMMY2   = ../gwnum/linux64/gwnum.ld

#########################################################################

all: $(EXE2)

$(EXE2): $(COMPOBJS) $(LIBELF)
	$(CC) $(LFLAGS) -o $(EXE2) $(COMPOBJS) $(DUMMY2) $(LIBS2)

$(COMPOBJS):
	make -f make64

clean:
	rm -f $(EXE2) $(COMPOBJS)
