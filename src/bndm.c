#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <string.h>
#include <sys/resource.h>
#include "globals.h"
#include "functions.h"

int bndm_search(unsigned char *pattern0, unsigned char *pattern1, unsigned int m) {
	unsigned int B[2][SIGMA];
	unsigned int S[2][SIGMA];
    //Stores starting register values for change from BNDM to SA. Basically one-set-bit on different positions 0 to m-1.
	unsigned int* R = (unsigned int*)malloc(m*sizeof(unsigned int));
	int i, j, F, D[2], D2[2], last, last_candidate[2], count, R1[2], R2[2];

	/* BADM Preprocessing */
	for (i = 0; i<SIGMA; i++){
        B[0][i] = 0;
        B[1][i] = 0;
    }
	F = 1;
	
	for (i = m - 1; i >= 0; i--) {
		//printf("bndm_search: 3, i = %d\n", i);
		B[0][pattern0[i]] |= F;
		B[1][pattern1[i]] |= F; // TODO use second pattern
		R[i] = F;
		F <<= 1;
	}
	F >>= 1;
	/*SA Preprocessing*/
	for (i = 0; i < SIGMA; ++i){
        S[0][i] = 0;
        S[1][i] = 0;
    }
    for (i = 0, j = 1; i < m; ++i, j <<= 1){
        S[0][pattern0[i]] |= j;
        S[1][pattern1[i]] |= j; // TODO use second pattern
    }

	/* Searching */

	unsigned char elementNum;
	unsigned int elementLength;
	unsigned int elementStart, elementEnd;
	unsigned int segmentCounter = 0;
	unsigned int elementCounter = 0;
		
	count = 0;
	R1[0] = R2[0] = D2[0] = 0;
	R1[1] = R2[1] = D2[1] = 0;

	while (aPointer < wbPointer)
	{
		//Processing the beginning of a segment
		elementNum = (unsigned char)writeBuffer[aPointer++];
		segmentCounter++;
		elementCounter += elementNum;
		R1[0] = R2[0];
		R1[1] = R2[2];
        R2[0] = 0;
        R2[1] = 0;

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
			D[0] = R1[0];
			D[1] = R1[0];
			while (j < (elementStart + m) && j < elementEnd) {
				D[0] = ((D[0] << 1) | 1) & S[0][writeBuffer[j]];
				D[1] = ((D[1] << 1) | 1) & S[1][writeBuffer[j]];

				if (D[0] & F) {
                    count++;
                    printf("SA HIT: j = %d, c = %c, D = %x, S = %x, elementStart = %d, m = %d, R1 = %x\n", j, writeBuffer[j], D[0], S[0][writeBuffer[j]], elementStart, m, R1[0]);
                }
                if (D[1] & F) {
                    count++;
                    printf("SA HIT: j = %d, c = %c, D = %x, S = %x, elementStart = %d, m = %d, R1 = %x\n", j, writeBuffer[j], D[1], S[1][writeBuffer[j]], elementStart, m, R1[0]);
                }
				j++;
			}
			if (elementLength < m)
				goto element_end;

			//Perform BNDM search.
			j = elementStart;
			while (j + m - 1 < elementEnd)
			{	
				D2[0] = 0;
				D2[1] = 0;
				last = m;
                last_candidate[0] = last;
                last_candidate[1] = last;
				i = m - 1;
				D[0] = ~0;
				D[1] = ~0;
				while (i >= 0 && (D[0] != 0 || D[1] != 0)) {
					D[0] &= B[0][writeBuffer[j + i]];
					D[1] &= B[1][writeBuffer[j + i]];
					i--;

                    int update_last = 0;
					if (D[0] != 0  && (D[0] & F) != 0) {
						if (i >= 0) { // Pure pattern prefix
                            last_candidate[0] = i + 1;
                            update_last = 1;

                        }
						else { // Match
                            count++;
                            printf("BNDM HIT: j = %d, i = %d, c = %c, D = %x, B = %x\n", j, i, writeBuffer[j + i], D[0], B[0][writeBuffer[j + i]]);
                        }
					}
					if (D[1] != 0  && (D[1] & F) != 0) {
						if (i >= 0) {
                            last_candidate[1] = i + 1;
                            update_last = 1;
                        }
						else {
                            count++;
                            printf("BNDM HIT: j = %d, i = %d, c = %c, D = %x, B = %x\n", j, i, writeBuffer[j + i], D[1], B[0][writeBuffer[j + i]]);
                        }
					}
                    if (update_last != 0){
                        last = min(last_candidate[0], last_candidate[1]);
                        D2[0] |= R[last];
                        D2[1] |= R[last];
                    }
					D[0] <<= 1;
					D[1] <<= 1;
				}

				j += last;
			}
			//printf("bndm_search 4: j = %d, last = %d, F = %x\n", j, last, F);

			//Perform SA search at the end of the element.
			j += m - last; //Moving j pointer to the initial position for SA.
			D[0] = D2[0]; //Setting the initial value to SA register.
			D[1] = D2[1]; //Setting the initial value to SA register.

//			if (j >= 6136576 && j <= 6136584)
//				printf("BNDM: j = %d, D = %x, last = %d, elementEnd = %d, D2 = %x, R[last] = %x\n", j, D, last, elementEnd, D2, R[last]);
			
			while(j < elementEnd) {
				D[0] = ((D[0] << 1) | 1) & S[0][writeBuffer[j]];
				D[1] = ((D[1] << 1) | 1) & S[1][writeBuffer[j]];
				//printf("bndm_search 5: j = %d, c = %c, D = %x, S = %x\n", j, writeBuffer[j], D, S[writeBuffer[j]]);

				if (D[0] & F)
				{
					count++;//This cannot happen...
					printf("SA2 HIT: j = %d, D = %x\n",j,D[0]);
				}
				if (D[1] & F)
				{
					count++;//This cannot happen...
					printf("SA2 HIT: j = %d, D = %x\n",j,D[1]);
				}
				j++;
			}

			element_end:
			R2[0] |= D[0]; //Merge the SA registers of single elements.
			R2[1] |= D[1]; //Merge the SA registers of single elements.

			//Set the pointer to the beginning of the next element.
			aPointer += elementLength;
		}
	}

	printf("count = %d, segments = %d, elements = %d\n", count, segmentCounter, elementCounter);

	return count;
}
