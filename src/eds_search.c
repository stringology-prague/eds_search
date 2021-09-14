#define _GNU_SOURCE

#include <stdlib.h>
#include <stdio.h>
#include <time.h>
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

const size_t MAX_DNA_PATTERNS = 2;

int main(int argc, char * argv[])
{
	rbPointer = 0;
	wbPointer = 0;
	aPointer = 0;

	unsigned int LOOPS = (unsigned int)atoi(argv[2]);
	size_t pattLen = (unsigned int)atoi(argv[3]);
	unsigned char algoritm = (unsigned char)atoi(argv[4]);
    unsigned char patterns[MAX_DNA_PATTERNS][MAX_PATTERN_LENGTH];
    unsigned char *pattern0 = (unsigned char*)calloc(MAX_PATTERN_LENGTH, sizeof(char));
    unsigned char *pattern1 = (unsigned char*)calloc(MAX_PATTERN_LENGTH, sizeof(char));

    if (argc >= 6){
        pattLen = strnlen(argv[5], MAX_PATTERN_LENGTH);
        strncpy(pattern0, argv[5], pattLen);
    }
    if (argc >= 7) {
        if (pattLen != strnlen(argv[6], MAX_PATTERN_LENGTH)) {
            fprintf(stderr, "Two patterns supplied, but length is not matching!");
            exit(1);
        }
        strncpy(pattern1, argv[6], pattLen);
    }

	readInputFile(argv[1], &readBuffer, &fSize);
	writeBuffer = malloc((fSize + 1000000) *sizeof(unsigned char*));
	//printf("fSize = %d\n", fSize);
	translate();
	writeOutputFile("text/test.out", &writeBuffer, wbPointer);

    init_IUPAC_SYMBOLS_TO_BASES();
    init_AA_TO_COMPR_IUPAC_SYMBOLS();
	
	struct rusage ruse, ruse1, ruse2;
	double ssec1, ssec2, usec1, usec2;
//	srand(time(NULL));
    srand(123);

	switch (algoritm) {

	case 1:
		goto BNDM;
		break;

	case 2:
		goto SA;
		break;

	case 3:
		goto BNDM_EDS_MP;
		break;

	case 4:
		goto BNDM_EDS_AA;
		break;

	default:
		printf("Invalid algorithm option!\n");
		return EXIT_FAILURE;
	}

BNDM_EDS_AA:;

//    if (*pattern0 == 0) {
//        fprintf(stderr, "BNDM-EDS-MP requires for AA pattern to be specified as argument!");
//        free(pattern0);
//        free(pattern1);
//        return 1;
//    }
//    int num_dna_patterns = translate_aa_pattern(pattern0, pattLen, patterns, MAX_DNA_PATTERNS, MAX_PATTERN_LENGTH);
//    free(pattern0);
//    free(pattern1);
//    int dna_pattern_length = pattLen * 3;
//    if (num_dna_patterns == 0){
//        fprintf(stderr, "BNDM-EDS-MP failed to generate any DNA patterns!");
//        return 1;
//    }
//
//    pattern0 = patterns[0];
//    if (num_dna_patterns == 2){
//        pattern1 = patterns[1];
//    } else {
//        pattern1 = patterns[0];
//    }

    // TODO remove - temporary for testing
    if (*pattern0 == 0 || *pattern1 == 0) {
        fprintf(stderr, "BNDM-EDS-AA requires two patterns consisting of IUPAC degenerate base symbols to be specified as argument!");
        free(pattern0);
        free(pattern1);
        return 1;
    }

    printf("BNDM-EDS-AA Pattern0: %s\n", pattern0);
    printf("BNDM-EDS-AA Pattern1: %s\n", pattern1);

	getrusage(RUSAGE_SELF, &ruse);
	ssec1 = (double)(ruse.ru_stime.tv_sec * 1000000 + ruse.ru_stime.tv_usec);
	usec1 = (double)(ruse.ru_utime.tv_sec * 1000000 + ruse.ru_utime.tv_usec);
	getrusage(RUSAGE_SELF, &ruse1);

	for (int i = 0; i < LOOPS; i++)
	{
		aPointer = 0;
		bndm_eds_aa_search(pattern0, pattern1, pattLen);
	}

	getrusage(RUSAGE_SELF, &ruse);
	ssec2 = (double)(ruse.ru_stime.tv_sec * 1000000 + ruse.ru_stime.tv_usec);
	usec2 = (double)(ruse.ru_utime.tv_sec * 1000000 + ruse.ru_utime.tv_usec);

	printf("User time:\t%f s\n", (usec2 - usec1) / (double)1000000);
	printf("System time:\t%f s\n", (ssec2 - ssec1) / (double)1000000);
	printf("Total time:\t%f s\n", ((usec2 + ssec2) - (usec1 + ssec1)) / (double)1000000);

	return 0;

BNDM_EDS_MP:;

    if (*pattern0 == 0) {
        fprintf(stderr, "BNDM-EDS-MP requires for AA pattern to be specified as argument!");
        free(pattern0);
        free(pattern1);
        return 1;
    }
    int num_dna_patterns = translate_aa_pattern(pattern0, pattLen, patterns, MAX_DNA_PATTERNS, MAX_PATTERN_LENGTH);
    free(pattern0);
    free(pattern1);
    int dna_pattern_length = pattLen * 3;
    if (num_dna_patterns == 0){
        fprintf(stderr, "BNDM-EDS-MP failed to generate any DNA patterns!");
        return 1;
    }

    pattern0 = patterns[0];
    if (num_dna_patterns == 2){
        pattern1 = patterns[1];
    } else {
        pattern1 = patterns[0];
    }
    printf("BNDM-EDS-MP Pattern0: %s\n", pattern0);
    printf("BNDM-EDS-MP Pattern1: %s\n", pattern1);

	getrusage(RUSAGE_SELF, &ruse);
	ssec1 = (double)(ruse.ru_stime.tv_sec * 1000000 + ruse.ru_stime.tv_usec);
	usec1 = (double)(ruse.ru_utime.tv_sec * 1000000 + ruse.ru_utime.tv_usec);
	getrusage(RUSAGE_SELF, &ruse1);

	for (int i = 0; i < LOOPS; i++)
	{
		aPointer = 0;
		bndm_eds_mp_search(pattern0, pattern1, dna_pattern_length);
	}

	getrusage(RUSAGE_SELF, &ruse);
	ssec2 = (double)(ruse.ru_stime.tv_sec * 1000000 + ruse.ru_stime.tv_usec);
	usec2 = (double)(ruse.ru_utime.tv_sec * 1000000 + ruse.ru_utime.tv_usec);

	printf("User time:\t%f s\n", (usec2 - usec1) / (double)1000000);
	printf("System time:\t%f s\n", (ssec2 - ssec1) / (double)1000000);
	printf("Total time:\t%f s\n", ((usec2 + ssec2) - (usec1 + ssec1)) / (double)1000000);

	return 0;

BNDM:;

	getrusage(RUSAGE_SELF, &ruse);
	ssec1 = (double)(ruse.ru_stime.tv_sec * 1000000 + ruse.ru_stime.tv_usec);
	usec1 = (double)(ruse.ru_utime.tv_sec * 1000000 + ruse.ru_utime.tv_usec);
	getrusage(RUSAGE_SELF, &ruse1);

	for (int i = 0; i < LOOPS; i++)
	{
        if (*pattern0 == 0){
            randomSelectPattern(pattern0, pattLen, readBuffer, fSize);
        }
		printf("Pattern: %s\n", pattern0);
		aPointer = 0;
		bndm_search(pattern0, pattLen);
	}

	getrusage(RUSAGE_SELF, &ruse);
	ssec2 = (double)(ruse.ru_stime.tv_sec * 1000000 + ruse.ru_stime.tv_usec);
	usec2 = (double)(ruse.ru_utime.tv_sec * 1000000 + ruse.ru_utime.tv_usec);

	printf("User time:\t%f s\n", (usec2 - usec1) / (double)1000000);
	printf("System time:\t%f s\n", (ssec2 - ssec1) / (double)1000000);
	printf("Total time:\t%f s\n", ((usec2 + ssec2) - (usec1 + ssec1)) / (double)1000000);

    free(pattern0);
    free(pattern1);
	return 0;

SA:;

	getrusage(RUSAGE_SELF, &ruse);
	ssec1 = (double)(ruse.ru_stime.tv_sec * 1000000 + ruse.ru_stime.tv_usec);
	usec1 = (double)(ruse.ru_utime.tv_sec * 1000000 + ruse.ru_utime.tv_usec);
	getrusage(RUSAGE_SELF, &ruse1);

	for (int i = 0; i < LOOPS; i++)
	{
        if (*pattern0 == 0){
            randomSelectPattern(pattern0, pattLen, readBuffer, fSize);
        }
		printf("Pattern: %s\n", pattern0);
		aPointer = 0;
		SA_search(pattern0, pattLen);
	}

	getrusage(RUSAGE_SELF, &ruse);
	ssec2 = (double)(ruse.ru_stime.tv_sec * 1000000 + ruse.ru_stime.tv_usec);
	usec2 = (double)(ruse.ru_utime.tv_sec * 1000000 + ruse.ru_utime.tv_usec);

	printf("User time:\t%f s\n", (usec2 - usec1) / (double)1000000);
	printf("System time:\t%f s\n", (ssec2 - ssec1) / (double)1000000);
	printf("Total time:\t%f s\n", ((usec2 + ssec2) - (usec1 + ssec1)) / (double)1000000);

    free(pattern0);
    free(pattern1);
	return 0;
}
