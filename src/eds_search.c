#define _GNU_SOURCE

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <sys/resource.h>
#include "globals.h"
#include "functions.h"
#include "bndm.h"
#include "bndm_eds_mp.h"
#include "bndm_aa.h"
#include "sa.h"
#include "translator.h"
#include "protein_table.h"


//Globals
unsigned char* readBuffer;
unsigned int fSize;
unsigned int rbPointer;

unsigned char* writeBuffer;
unsigned int wbPointer;
unsigned int aPointer;

int main(int argc, char * argv[])
{
	rbPointer = 0;
	wbPointer = 0;
	aPointer = 0;

	int LOOPS = atoi(argv[2]);
	size_t pattLen = (unsigned int)atoi(argv[3]), patt2Len = 0;
	unsigned char algorithm = (unsigned char)atoi(argv[4]);
    unsigned char patterns[MAX_DNA_PATTERNS][MAX_PATTERN_LENGTH+1];
    memset(patterns, 0, MAX_DNA_PATTERNS * (MAX_PATTERN_LENGTH + 1));

    if (argc >= 6){
        pattLen = strnlen(argv[5], MAX_PATTERN_LENGTH+1);
        memcpy(patterns[0], argv[5], pattLen);
    }
    if (argc >= 7) {
        patt2Len = strnlen(argv[6], MAX_PATTERN_LENGTH+1);
        strncpy(patterns[1], argv[6], pattLen);
    }
    if (pattLen && patt2Len && pattLen != patt2Len) {
        fprintf(stderr, "Two patterns supplied, but length is not matching!");
        exit(1);
    }

	readInputFile(argv[1], &readBuffer, &fSize);
	writeBuffer = malloc((fSize + 1000000) *sizeof(unsigned char*));
	//printf("fSize = %d\n", fSize);
	translate();
	writeOutputFile("text/test.out", &writeBuffer, wbPointer);

    init_IUPAC_SYMBOLS_TO_BASES();
    init_AA_TO_COMPR_IUPAC_SYMBOLS();

//	srand(time(NULL));
    srand(123);

	switch (algorithm) {

	case 1:
        SA_run(patterns[0], pattLen, LOOPS);
		break;

	case 2:
        bndm_eds_run(patterns[0], pattLen, LOOPS);
		break;

	case 3:
        bndm_eds_mp_run(patterns[0], pattLen, LOOPS);
		break;

	case 4:
        bndm_eds_aa_run(patterns[0], pattLen, LOOPS);
		break;

	default:
		fprintf(stderr,"Invalid algorithm option!\n");
		return EXIT_FAILURE;
	}

	return 0;
}
