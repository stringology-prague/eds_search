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
#include "sa.h"
#include "translator.h"



//Globals
unsigned char* readBuffer;
unsigned int fSize;
unsigned int rbPointer;

unsigned char* writeBuffer;
unsigned int wbPointer;
unsigned int aPointer;

const size_t MAX_PATTERN_LENGTH = sizeof(int) * 8;

int main(int argc, char * argv[])
{
	rbPointer = 0;
	wbPointer = 0;
	aPointer = 0;

	unsigned int LOOPS = (unsigned int)atoi(argv[2]);
	size_t pattLen = (unsigned int)atoi(argv[3]);
	unsigned char algoritm = (unsigned char)atoi(argv[4]);
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

	default:
		printf("Invalid algorithm option!\n");
		return EXIT_FAILURE;
	}

BNDM_EDS_MP:;

	getrusage(RUSAGE_SELF, &ruse);
	ssec1 = (double)(ruse.ru_stime.tv_sec * 1000000 + ruse.ru_stime.tv_usec);
	usec1 = (double)(ruse.ru_utime.tv_sec * 1000000 + ruse.ru_utime.tv_usec);
	getrusage(RUSAGE_SELF, &ruse1);

	for (int i = 0; i < LOOPS; i++)
	{
        if (*pattern0 == 0){
            randomSelectPattern(pattern0, pattLen, readBuffer, fSize);
        }
        if (*pattern1 == 0){
            randomSelectPattern(pattern1, pattLen, readBuffer, fSize);
        }
		printf("Pattern0: %s\n", pattern0);
		printf("Pattern1: %s\n", pattern1);
		aPointer = 0;
		bndm_eds_mp_search(pattern0, pattern1, pattLen);
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
