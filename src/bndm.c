#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <string.h>
#include <sys/resource.h>
#include "globals.h"
#include "functions.h"

int bndm_search(unsigned char *x, unsigned int m) {
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

	unsigned char elementNum;
	unsigned int elementLength;
	unsigned int elementStart, elementEnd;
	unsigned int segmentCounter = 0;
	unsigned int elementCounter = 0;
		
	count = 0;
	R1 = R2 = D2 = 0;
	

	while (aPointer < wbPointer)
	{
		//Processing the beginning of a segment
		elementNum = (unsigned char)writeBuffer[aPointer++];
		segmentCounter++;
		elementCounter += elementNum;
		R1 = R2; R2 = 0;

		//Process one element of the segment.
		for (unsigned char k = 0; k < elementNum; k++)
		{
			//elementLenght = (unsigned char)writeBuffer[aPointer++]; elementLenght <<= 8;
			//elementLenght |= (unsigned char)writeBuffer[aPointer++];
			elementLength = byteDecodeInt();
			elementStart = aPointer;
			elementEnd = elementStart + elementLength;
			//printf("bndm_search 1: k = %d, elementNum = %d, elementLength = %d, elementStart = %d, elementEnd = %d\n", k, elementNum, elementLenght, elementStart, elementEnd);

			//Perform SA search at the beginning of the element.
			j = elementStart;
			D = R1;
			while (j < (elementStart + m) && j < elementEnd) {
				D = ((D << 1) | 1) & S[writeBuffer[j]];

				

				if (D & F) { count++; printf("SA HIT: j = %d, c = %c, D = %x, S = %x, elementStart = %d, m = %d, R1 = %x\n", j, writeBuffer[j], D, S[writeBuffer[j]], elementStart, m, R1);}
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
					D &= B[writeBuffer[j + i]];

					

					i--;
					if (D != 0  && (D & F) != 0) {
						if (i >= 0) { last = i + 1; D2 |= R[last]; }
						else { count++; printf("BNDM HIT: j = %d, i = %d, c = %c, D = %x, B = %x\n", j, i, writeBuffer[j + i], D, B[writeBuffer[j + i]]);}
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
				printf("BNDM: j = %d, D = %x, last = %d, elementEnd = %d, D2 = %x, R[last] = %x\n", j, D, last, elementEnd, D2, R[last]);
			
			while(j < elementEnd) {
				D = ((D << 1) | 1) & S[writeBuffer[j]];
				//printf("bndm_search 5: j = %d, c = %c, D = %x, S = %x\n", j, writeBuffer[j], D, S[writeBuffer[j]]);

				if (D & F)
				{
					count++;//This cannot happen...
					printf("SA2 HIT: j = %d, D = %x\n",j,D);
				}
				j++;
			}

			element_end:
			R2 |= D; //Merge the SA registers of single elements.

			//Set the pointer to the beginning of the next element.
			aPointer += elementLength;
		}
	}

	printf("count = %d, segments = %d, elements = %d\n", count, segmentCounter, elementCounter);

	return count;
}
