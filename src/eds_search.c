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

int main(int argc, char *argv[]) {

    int LOOPS = atoi(argv[2]);
    size_t pattLen = (unsigned int) atoi(argv[3]), patt2Len = 0;
    const char *algorithm = argv[4];
    unsigned char patterns[MAX_DNA_PATTERNS][MAX_PATTERN_LENGTH + 1];
    memset(patterns, 0, MAX_DNA_PATTERNS*(MAX_PATTERN_LENGTH + 1));

    if (argc >= 6 && *argv[5]) {
        pattLen = strnlen(argv[5], MAX_PATTERN_LENGTH + 1);
        memcpy(patterns[0], argv[5], pattLen);
    }
    if (argc >= 7 && *argv[6]) {
        patt2Len = strnlen(argv[6], MAX_PATTERN_LENGTH + 1);
        strncpy(patterns[1], argv[6], pattLen);
    }
    if (pattLen && patt2Len && pattLen!=patt2Len) {
        fprintf(stderr, "Two patterns supplied, but length is not matching!\n");
        return EXIT_FAILURE;
    }

    int eds_size, teds_size;
    unsigned char *eds = readInputFile(argv[1], &eds_size);
    unsigned char *teds = translate(eds, eds_size, &teds_size);
    free(eds);

    writeOutputFile("text/test.out", &teds, teds_size);

    init_IUPAC_SYMBOLS_TO_BASES();
    init_AA_TO_COMPR_IUPAC_SYMBOLS();

    srand(time(NULL));

    if (strcasecmp(algorithm, "sa")==0) {
        SA_run(teds, teds_size, patterns[0], pattLen, LOOPS);
    } else if (strcasecmp(algorithm, "bndm")==0) {
        bndm_eds_run(teds, teds_size, patterns[0], pattLen, LOOPS);
    } else if (strcasecmp(algorithm, "bndm-mp")==0) {
        bndm_eds_mp_run(teds, teds_size, patterns[0], pattLen, LOOPS);
    } else if (strcasecmp(algorithm, "bndm-aa")==0) {
        bndm_eds_aa_run(teds, teds_size, patterns[0], patterns[1], pattLen, LOOPS);
    } else {
        fprintf(stderr, "Invalid algorithm option!\n");
        free(teds);
        return EXIT_FAILURE;
    }
    free(teds);
    return EXIT_SUCCESS;
}
