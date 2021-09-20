#define _GNU_SOURCE

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>

#include "globals.h"
#include "functions.h"
#include "bndm.h"
#include "bndm_eds_mp.h"
#include "bndm_aa.h"
#include "sa.h"
#include "translator.h"
#include "protein_table.h"

//Globals
unsigned char *readBuffer;
unsigned int fSize;
unsigned int rbPointer;

unsigned char *writeBuffer;
unsigned int wbPointer;
unsigned int aPointer;

int main(int argc, char *argv[]) {
    rbPointer = 0;
    wbPointer = 0;
    aPointer = 0;

    int LOOPS = atoi(argv[2]);
    size_t pattLen = (unsigned int) atoi(argv[3]), patt2Len = 0;
    const char *algorithm = argv[4];
    unsigned char patterns[MAX_DNA_PATTERNS][MAX_PATTERN_LENGTH + 1];
    memset(patterns, 0, MAX_DNA_PATTERNS*(MAX_PATTERN_LENGTH + 1));

    if (argc >= 6) {
        pattLen = strnlen(argv[5], MAX_PATTERN_LENGTH + 1);
        memcpy(patterns[0], argv[5], pattLen);
    }
    if (argc >= 7) {
        patt2Len = strnlen(argv[6], MAX_PATTERN_LENGTH + 1);
        strncpy(patterns[1], argv[6], pattLen);
    }
    if (pattLen && patt2Len && pattLen!=patt2Len) {
        fprintf(stderr, "Two patterns supplied, but length is not matching!");
        return EXIT_FAILURE;
    }

    readInputFile(argv[1], &readBuffer, &fSize);
    writeBuffer = malloc((fSize + 1000000)*sizeof(unsigned char *));
    //printf("fSize = %d\n", fSize);
    translate();
    writeOutputFile("text/test.out", &writeBuffer, wbPointer);

    init_IUPAC_SYMBOLS_TO_BASES();
    init_AA_TO_COMPR_IUPAC_SYMBOLS();

    srand(time(NULL));

    if (strcasecmp(algorithm, "sa")==0) {
        SA_run(patterns[0], pattLen, LOOPS);
    } else if (strcasecmp(algorithm, "bndm")==0) {
        bndm_eds_run(patterns[0], pattLen, LOOPS);
    } else if (strcasecmp(algorithm, "bndm-mp")==0) {
        bndm_eds_mp_run(patterns[0], pattLen, LOOPS);
    } else if (strcasecmp(algorithm, "bndm-aa")==0) {
        bndm_eds_aa_run(patterns[0], pattLen, LOOPS);
    } else {
        fprintf(stderr, "Invalid algorithm option!\n");
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}
