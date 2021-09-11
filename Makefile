SRC = ./src
CC = gcc
override CFLAGS += -O3 -mavx -msse2

all:	eds_search clean0

eds_search: eds_search.o translator.o functions.o bndm.o bndm_eds_mp.o protein_table.o sa.o $(SRC)/globals.h
	$(CC) $(CFLAGS) -o eds_search eds_search.o translator.o functions.o bndm.o sa.o bndm_eds_mp.o protein_table.o

eds_search.o: $(SRC)/eds_search.c
	$(CC) $(CFLAGS) -c $(SRC)/eds_search.c

translator.o: $(SRC)/translator.c
	$(CC) $(CFLAGS) -c $(SRC)/translator.c

functions.o: $(SRC)/functions.c
	$(CC) $(CFLAGS) -c $(SRC)/functions.c

bndm.o: $(SRC)/bndm.c
	$(CC) $(CFLAGS) -c $(SRC)/bndm.c

bndm_eds_mp.o: $(SRC)/bndm_eds_mp.c
	$(CC) $(CFLAGS) -c $(SRC)/bndm_eds_mp.c

protein_table.o: $(SRC)/protein_table.c
	$(CC) $(CFLAGS) -c $(SRC)/protein_table.c

sa.o: $(SRC)/sa.c
	$(CC) $(CFLAGS) -c $(SRC)/sa.c

clean0:
	rm -f *.o
