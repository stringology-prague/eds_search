#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <string.h>
#include <sys/resource.h>

#include "globals.h"
#include "functions.h"

#include "bndm.h"

int bndm_eds_run(const unsigned char *teds,
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

    printf("BNDM-EDS\n");
    printf("Pattern:\t\"%.*s\"\n", (int) m, pattern);

    getrusage(RUSAGE_SELF, &ruse);
    ssec1 = (double)(ruse.ru_stime.tv_sec * 1000000 + ruse.ru_stime.tv_usec);
    usec1 = (double)(ruse.ru_utime.tv_sec * 1000000 + ruse.ru_utime.tv_usec);
    getrusage(RUSAGE_SELF, &ruse1);

    for (int i = 0; i < loops; i++)
    {
        if (m == 0){
            randomSelectPattern(rand_pattern, m, teds, len);
        }
        DEBUG_PRINT("Pattern: %.*s\n", (int)m, rand_pattern);
        matches = bndm_search(teds, len, rand_pattern, m);
    }

    getrusage(RUSAGE_SELF, &ruse);
    ssec2 = (double)(ruse.ru_stime.tv_sec * 1000000 + ruse.ru_stime.tv_usec);
    usec2 = (double)(ruse.ru_utime.tv_sec * 1000000 + ruse.ru_utime.tv_usec);

    printf("Matches:\t%d\n", matches);
    printf("User time:\t%f s\n", (usec2 - usec1) / (double)1000000);
    printf("System time:\t%f s\n", (ssec2 - ssec1) / (double)1000000);
    printf("Total time:\t%f s\n", ((usec2 + ssec2) - (usec1 + ssec1)) / (double)1000000);

    return matches;
}

int bndm_search(const unsigned char *teds,
                const size_t len,
                unsigned char *x,
                unsigned int m) {
	unsigned int B[SIGMA];
	unsigned int S[SIGMA];
	unsigned int* R = (unsigned int*)malloc(m*sizeof(unsigned int)); //Stores starting register values for change from BNDM to SA. Basically one-set-bit on different positions 0 to m-1.
	int i, j, F, D, D2, last, count, R1, R2;

	/* BADM Preprocessing */
	for (i = 0; i<SIGMA; i++) B[i] = 0;
	F = 1;
	
	for (i = m - 1; i >= 0; i--) {
		//printf("bndm_search: 3, i = %d\n", i);
		B[x[i]] |= F;
		R[i] = F;
		F <<= 1;
	}
	F >>= 1;
	/*SA Preprocessing*/
	for (i = 0; i < SIGMA; ++i) S[i] = 0;
		for (i = 0, j = 1; i < m; ++i, j <<= 1)
			S[x[i]] |= j;
	
	/* Searching */

	unsigned int elementNum;
	unsigned int elementLength;
	unsigned int elementStart, elementEnd;
	unsigned int segmentCounter = 0;
	unsigned int elementCounter = 0;
		
	count = 0;
	R1 = R2 = D2 = 0;
	
    int aPointer = 0;
	while (aPointer < len)
	{
		//Processing the beginning of a segment
        aPointer += byteDecodeInt(teds + aPointer, &elementNum);
		segmentCounter++;
		elementCounter += elementNum;
		R1 = R2; R2 = 0;

		//Process one element of the segment.
		for (unsigned int k = 0; k < elementNum; k++)
		{
			//elementLenght = (unsigned char)teds[aPointer++]; elementLenght <<= 8;
			//elementLenght |= (unsigned char)teds[aPointer++];
            aPointer += byteDecodeInt(teds + aPointer, &elementLength);
			elementStart = aPointer;
			elementEnd = elementStart + elementLength;
			//printf("bndm_search 1: k = %d, elementNum = %d, elementLength = %d, elementStart = %d, elementEnd = %d\n", k, elementNum, elementLenght, elementStart, elementEnd);

			//Perform SA search at the beginning of the element.
			j = elementStart;
			D = R1;
			while (j < (elementStart + m) && j < elementEnd) {
				D = ((D << 1) | 1) & S[teds[j]];

				

				if (D & F) { count++;DEBUG_PRINT("SA HIT: j = %d, c = %c, D = %x, S = %x, elementStart = %d, m = %d, R1 = %x\n", j, teds[j], D, S[teds[j]], elementStart, m, R1);}
				j++;
			}
			if (elementLength < m)
				goto element_end;

			//Perform BNDM search.
			j = elementStart;
			while (j + m - 1 < elementEnd)
			{	
				D2 = 0;
				last = m;
				i = m - 1;
				D = ~0;
				while (i >= 0 && D != 0) {
					D &= B[teds[j + i]];

					

					i--;
					if (D != 0  && (D & F) != 0) {
						if (i >= 0) { last = i + 1; D2 |= R[last]; }
						else { count++; DEBUG_PRINT("BNDM HIT: j = %d, i = %d, c = %c, D = %x, B = %x\n", j, i, teds[j + i], D, B[teds[j + i]]);}
					}
					D <<= 1;
				}
				j += last;
			}
			//printf("bndm_search 4: j = %d, last = %d, F = %x\n", j, last, F);


			//Perform SA search at the end of the element.
			
			j += m - last; //Moving j pointer to the initial position for SA.
			D = D2; //Setting the initial value to SA register.
			


			if (j >= 6136576 && j <= 6136584)
                DEBUG_PRINT("BNDM: j = %d, D = %x, last = %d, elementEnd = %d, D2 = %x, R[last] = %x\n", j, D, last, elementEnd, D2, R[last]);
			
			while(j < elementEnd) {
				D = ((D << 1) | 1) & S[teds[j]];
				//printf("bndm_search 5: j = %d, c = %c, D = %x, S = %x\n", j, teds[j], D, S[teds[j]]);

				if (D & F)
				{
					count++;//This cannot happen...
                    DEBUG_PRINT("SA2 HIT: j = %d, D = %x\n",j,D);
				}
				j++;
			}

			element_end:
			R2 |= D; //Merge the SA registers of single elements.

			//Set the pointer to the beginning of the next element.
			aPointer += elementLength;
		}
	}

    DEBUG_PRINT("count = %d, segments = %d, elements = %d\n", count, segmentCounter, elementCounter);

	return count;
}
