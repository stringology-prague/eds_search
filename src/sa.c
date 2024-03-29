#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <string.h>
#include <sys/resource.h>

#include "globals.h"
#include "functions.h"

#include "sa.h"

int SA_run(const unsigned char *teds,
           const size_t len,
           const unsigned char *pattern,
           const size_t m,
           const int loops)
{
    struct rusage ruse, ruse1, ruse2;
    double ssec1, ssec2, usec1, usec2;
    int matches;
    unsigned char rand_pattern[MAX_PATTERN_LENGTH];
    memcpy(rand_pattern, pattern, m);

    getrusage(RUSAGE_SELF, &ruse);
    ssec1 = (double)(ruse.ru_stime.tv_sec * 1000000 + ruse.ru_stime.tv_usec);
    usec1 = (double)(ruse.ru_utime.tv_sec * 1000000 + ruse.ru_utime.tv_usec);
    getrusage(RUSAGE_SELF, &ruse1);

    for (int i = 0; i < loops; i++)
    {
        if (m == 0){
            randomSelectPattern(rand_pattern, m, teds, len);
        }
        printf("Pattern: %s\n", rand_pattern);
        matches = SA_search(teds, len, rand_pattern, m);
    }

    getrusage(RUSAGE_SELF, &ruse);
    ssec2 = (double)(ruse.ru_stime.tv_sec * 1000000 + ruse.ru_stime.tv_usec);
    usec2 = (double)(ruse.ru_utime.tv_sec * 1000000 + ruse.ru_utime.tv_usec);

    printf("User time:\t%f s\n", (usec2 - usec1) / (double)1000000);
    printf("System time:\t%f s\n", (ssec2 - ssec1) / (double)1000000);
    printf("Total time:\t%f s\n", ((usec2 + ssec2) - (usec1 + ssec1)) / (double)1000000);

    return matches;
}

int SA_search(const unsigned char *teds,
              const size_t len,
              unsigned char *x,
              unsigned int m) {
	unsigned int S[SIGMA];
	int i, j, F, D, R1, R2, last, count;

	/*SA Preprocessing*/
	for (i = 0; i < SIGMA; ++i) S[i] = 0;
	for (i = 0, j = 1; i < m; ++i, j <<= 1)
		S[x[i]] |= j;

	F = 1; F <<= m - 1;
	/* Searching */

	unsigned int elementNum;
	unsigned int elementLength;
	unsigned int elementStart, elementEnd;
	unsigned int segmentCounter = 0;
	unsigned int elementCounter = 0;

	count = 0;
	R1 = R2 = 0;

    int aPointer = 0;
	while (aPointer < len)
	{
		//Processing the beginning of a segment
        aPointer += byteDecodeInt(teds + aPointer, &elementNum);
		segmentCounter++;
		elementCounter += elementNum;
			
		//Process one element of the segment.
		for (unsigned int k = 0; k < elementNum; k++)
		{
            aPointer += byteDecodeInt(teds + aPointer, &elementLength);
			elementStart = aPointer;
			elementEnd = elementStart + elementLength;
			//printf("bndm_search 1: k = %d, elementNum = %d, elementLength = %d, elementStart = %d, elementEnd = %d\n", k, elementNum, elementLenght, elementStart, elementEnd);

			//Perform SA search at the beginning of the element.
			j = elementStart;
			D = R1;
			while (j < elementEnd) {
				D = ((D << 1) | 1) & S[teds[j]];



				if (D & F) { count++; printf("SA HIT: j = %d, c = %c, D = %x, S = %x, elementStart = %d, m = %d\n", j, teds[j], D, S[teds[j]], elementStart, m); }
				j++;
			}
			
			R2 |= D; //Merge the SA registers of single elements.
			aPointer += elementLength;
		}
		R1 = R2;
		R2 = 0;
	}

	printf("count = %d, segments = %d, elements = %d\n", count, segmentCounter, elementCounter);

	return count;
}