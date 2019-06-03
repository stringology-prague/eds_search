#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <string.h>
#include <sys/resource.h>
#include "globals.h"
#include "functions.h"
#include "bndm.h"
#include "sa.h"
#include "translator.h"

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

	unsigned int LOOPS = (unsigned int)atoi(argv[2]);
	unsigned int pattLen = (unsigned int)atoi(argv[3]);
	unsigned char algoritm = (unsigned char)atoi(argv[4]);

	readInputFile(argv[1], &readBuffer, &fSize);
	writeBuffer = malloc((fSize + 1000000) *sizeof(unsigned char*));
	//printf("fSize = %d\n", fSize);
	translate();
	writeOutputFile("text/test.out", &writeBuffer, wbPointer);
		
	unsigned char* pattern;
	
	struct rusage ruse, ruse1, ruse2;
	double ssec1, ssec2, usec1, usec2;
	srand(time(NULL));

	switch (algoritm) {

	case 1:
		goto BNDM;
		break;

	case 2:
		goto SA;
		break;

	default:
		printf("Invalid algorithm option!\n");
		return EXIT_FAILURE;
	}


BNDM:;
	
	getrusage(RUSAGE_SELF, &ruse);
	ssec1 = (double)(ruse.ru_stime.tv_sec * 1000000 + ruse.ru_stime.tv_usec);
	usec1 = (double)(ruse.ru_utime.tv_sec * 1000000 + ruse.ru_utime.tv_usec);
	getrusage(RUSAGE_SELF, &ruse1);

	for (int i = 0; i < LOOPS; i++)
	{
		randomSelectPattern(&pattern, pattLen, readBuffer, fSize);
		printf("Pattern: %s\n", pattern);
		aPointer = 0;
		bndm_search(pattern, pattLen);
		free(pattern);
	}

	getrusage(RUSAGE_SELF, &ruse);
	ssec2 = (double)(ruse.ru_stime.tv_sec * 1000000 + ruse.ru_stime.tv_usec);
	usec2 = (double)(ruse.ru_utime.tv_sec * 1000000 + ruse.ru_utime.tv_usec);

	printf("User time:\t%f s\n", (usec2 - usec1) / (double)1000000);
	printf("System time:\t%f s\n", (ssec2 - ssec1) / (double)1000000);
	printf("Total time:\t%f s\n", ((usec2 + ssec2) - (usec1 + ssec1)) / (double)1000000);
	return 0;

SA:;

	getrusage(RUSAGE_SELF, &ruse);
	ssec1 = (double)(ruse.ru_stime.tv_sec * 1000000 + ruse.ru_stime.tv_usec);
	usec1 = (double)(ruse.ru_utime.tv_sec * 1000000 + ruse.ru_utime.tv_usec);
	getrusage(RUSAGE_SELF, &ruse1);

	for (int i = 0; i < LOOPS; i++)
	{
		randomSelectPattern(&pattern, pattLen, readBuffer, fSize);
		printf("Pattern: %s\n", pattern);
		aPointer = 0;
		SA_search(pattern, pattLen);
		free(pattern);
	}

	getrusage(RUSAGE_SELF, &ruse);
	ssec2 = (double)(ruse.ru_stime.tv_sec * 1000000 + ruse.ru_stime.tv_usec);
	usec2 = (double)(ruse.ru_utime.tv_sec * 1000000 + ruse.ru_utime.tv_usec);

	printf("User time:\t%f s\n", (usec2 - usec1) / (double)1000000);
	printf("System time:\t%f s\n", (ssec2 - ssec1) / (double)1000000);
	printf("Total time:\t%f s\n", ((usec2 + ssec2) - (usec1 + ssec1)) / (double)1000000);


	return 0;
}
