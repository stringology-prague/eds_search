#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <string.h>
#include <sys/resource.h>
#include "globals.h"
#include "functions.h"

int SA_search(unsigned char *x, unsigned int m) {
	unsigned int S[SIGMA];
	int i, j, F, D, R1, R2, last, count;

	/*SA Preprocessing*/
	for (i = 0; i < SIGMA; ++i) S[i] = 0;
	for (i = 0, j = 1; i < m; ++i, j <<= 1)
		S[x[i]] |= j;

	F = 1; F <<= m - 1;
	/* Searching */

	unsigned char elementNum;
	unsigned int elementLength;
	unsigned int elementStart, elementEnd;
	unsigned int segmentCounter = 0;
	unsigned int elementCounter = 0;

	count = 0;
	R1 = R2 = 0;

	while (aPointer < wbPointer)
	{
		//Processing the beginning of a segment
		elementNum = (unsigned char)writeBuffer[aPointer++];
		segmentCounter++;
		elementCounter += elementNum;
			
		//Process one element of the segment.
		for (unsigned char k = 0; k < elementNum; k++)
		{
			elementLength = byteDecodeInt();
			elementStart = aPointer;
			elementEnd = elementStart + elementLength;
			//printf("bndm_search 1: k = %d, elementNum = %d, elementLength = %d, elementStart = %d, elementEnd = %d\n", k, elementNum, elementLenght, elementStart, elementEnd);

			//Perform SA search at the beginning of the element.
			j = elementStart;
			D = R1;
			while (j < elementEnd) {
				D = ((D << 1) | 1) & S[writeBuffer[j]];



				if (D & F) { count++; printf("SA HIT: j = %d, c = %c, D = %x, S = %x, elementStart = %d, m = %d\n", j, writeBuffer[j], D, S[writeBuffer[j]], elementStart, m); }
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