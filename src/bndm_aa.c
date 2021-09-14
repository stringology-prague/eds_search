#include <assert.h>
#include <stdlib.h>
#include <stdio.h>

#include "globals.h"
#include "functions.h"
#include "protein_table.h"

int bndm_eds_aa_search(unsigned char *pattern0, unsigned char *pattern1, unsigned int m) {
    assert(m % 3 == 0);
    unsigned int B[2][SIGMA];
    unsigned int S[2][SIGMA];

    // ShiftAnd Dual-Pattern Mask - used to mask incoming states from the other pattern
    const unsigned int SA_DP_MASK = 0x24924924U;
    // const unsigned long int SA_DP_MASK = 0x4924924924924924UL;

    // Stores starting register values for change from BNDM to SA. Basically one-set-bit on different positions 0 to m-1.
    unsigned int* R[2];
    R[0] = (unsigned int*)malloc(m*sizeof(unsigned int));
    R[1] = (unsigned int*)malloc(m*sizeof(unsigned int));
    int i, j, F, D[2], D2[2], last, last_candidate[2], count, R1[2], R2[2];

    /* BADM Preprocessing */
    for (i = 0; i<SIGMA; i++){
        B[0][i] = 0;
        B[1][i] = 0;
    }
    F = 1;
    for (i = m - 1; i >= 0; i--) {

        for (int b = 0; b < BASES; b++){
            B[0][IUPAC_SYMBOLS_TO_BASES[pattern0[i]][b]] |= F;
            B[1][IUPAC_SYMBOLS_TO_BASES[pattern1[i]][b]] |= F;
        }
        R[0][i] = F;
        R[1][i] = F;
        F <<= 1;
    }
    F >>= 1;

    /*SA Preprocessing*/
    for (i = 0; i < SIGMA; ++i){
        S[0][i] = 0;
        S[1][i] = 0;
    }
    for (i = 0, j = 1; i < m; ++i, j <<= 1){
        for (int b = 0; b < BASES; b++){
            S[0][IUPAC_SYMBOLS_TO_BASES[pattern0[i]][b]] |= j;
            S[1][IUPAC_SYMBOLS_TO_BASES[pattern1[i]][b]] |= j;
        }
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
        R1[1] = R2[1];
        R2[0] = 0;
        R2[1] = 0;

        //Process one element of the segment.
        for (unsigned char k = 0; k < elementNum; k++)
        {
            elementLength = byteDecodeInt();
            elementStart = aPointer;
            elementEnd = elementStart + elementLength;
            DEBUG_PRINT("  ELEMENT segment = %d, element = %d/%d, start = %d, end = %d, len = %d, elementCounter = %d\n",
                        segmentCounter, k, elementNum, elementStart, elementEnd, elementLength, elementCounter);

            //Perform SA search at the beginning of the element.
            j = elementStart;
            D[0] = R1[0];
            D[1] = R1[1];
            while (j < (elementStart + m) && j < elementEnd) {
                char curr_symbol = writeBuffer[j];
                int DC[2]; // Vector states taken from the other pattern (for every 3 characters only)
                DC[0] = D[1] & SA_DP_MASK;
                DC[1] = D[0] & SA_DP_MASK;
                D[0] = (((D[0] | DC[0]) << 1) | 1) & S[0][curr_symbol];
                D[1] = (((D[1] | DC[1]) << 1) | 1) & S[1][curr_symbol];
                DEBUG_PRINT("    SA1 j=%d, curr_symbol=%c, D[0]=0x%x, D[1]=0x%x\n", j, curr_symbol, D[0], D[1]);

                if (D[0] & F) {
                    count++;
                    DEBUG_PRINT("      ");
                    printf("SA HIT (pattern0): j = %d, c = %c, D[0] = 0x%x, S = 0x%x, elementStart = %d, m = %d, R1 = 0x%x\n", j, curr_symbol, D[0], S[0][curr_symbol], elementStart, m, R1[0]);
                }
                if (D[1] & F) {
                    count++;
                    DEBUG_PRINT("      ");
                    printf("SA HIT (pattern1): j = %d, c = %c, D[1] = 0x%x, S = 0x%x, elementStart = %d, m = %d, R1 = 0x%x\n", j, curr_symbol, D[1], S[1][curr_symbol], elementStart, m, R1[1]);
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
                DEBUG_PRINT("    BNDM j=%d, i=%d\n", j, i);
                while (i >= 0 && (D[0] != 0 || D[1] != 0)) {
                    char curr_symbol = writeBuffer[j + i];

                    if ((m - 1 - i) % 3 == 0) { // Every 3 symbols, we can use a state from either pattern
                        D[0] = (D[0] | D[1]) & B[0][curr_symbol];
                        D[1] = (D[0] | D[1]) & B[1][curr_symbol];
                    } else {
                        D[0] = D[0] & B[0][curr_symbol];
                        D[1] = D[1] & B[1][curr_symbol];
                    }

                    DEBUG_PRINT("      BNDM j=%d, i=%d, j+i=%d, curr_symbol=%c, D[0]=0x%x, D[1]=0x%x\n", j, i, j+i, curr_symbol, D[0], D[1]);
                    i--;

                    if (D[0] != 0  && (D[0] & F) != 0) {
                        if (i >= 0) { // Pure pattern prefix
                            last_candidate[0] = i + 1;
                            D2[0] |= R[0][last_candidate[0]];
                        }
                        else { // Match
                            count++;
                            DEBUG_PRINT("        ");
                            printf("BNDM HIT (pattern0): j = %d, i = %d, c = %c, D = 0x%x, B = 0x%x\n", j, i, curr_symbol, D[0], B[0][curr_symbol]);
                        }
                    }
                    if (D[1] != 0  && (D[1] & F) != 0) {
                        if (i >= 0) {
                            last_candidate[1] = i + 1;
                            D2[1] |= R[1][last_candidate[1]];
                        }
                        else {
                            count++;
                            DEBUG_PRINT("        ");
                            printf("BNDM HIT (pattern1): j = %d, i = %d, c = %c, D = 0x%x, B = 0x%x\n", j, i, curr_symbol, D[1], B[1][curr_symbol]);
                        }
                    }
                    D[0] <<= 1;
                    D[1] <<= 1;
                }
                // TODO Consier if `last` also needs to be reduced to
                last = min(last_candidate[0], last_candidate[1]);
                j += last;
            }
            //Perform SA search at the end of the element.
            DEBUG_PRINT("  BNDM->SA2 j=%d last=%d, jnext = %d, D[0]=0x%x, D[1]=0x%x\n", j, last, j + m - last, D2[0], D2[1]);
            j += m - last; //Moving j pointer to the initial position for SA.
            D[0] = D2[0]; //Setting the initial value to SA register.
            D[1] = D2[1]; //Setting the initial value to SA register.
            while(j < elementEnd) {
                char curr_symbol = writeBuffer[j];
                int DC[2]; // Vector states taken from the other pattern (for every 3 characters only)
                DC[0] = D[1] & SA_DP_MASK;
                DC[1] = D[0] & SA_DP_MASK;
                D[0] = (((D[0] | DC[0]) << 1) | 1) & S[0][curr_symbol];
                D[1] = (((D[1] | DC[1]) << 1) | 1) & S[1][curr_symbol];
                DEBUG_PRINT("    SA2 j=%d, curr_symbol=%c, D[0]=0x%x, D[1]=0x%x\n", j, curr_symbol, D[1], D[2]);

                if (D[0] & F)
                {
                    count++;//This cannot happen...
                    DEBUG_PRINT("      SA2 HIT (pattern0): j = %d, D[0] = 0x%x\n",j,D[0]);
                }
                if (D[1] & F)
                {
                    count++;//This cannot happen...
                    DEBUG_PRINT("      SA2 HIT (pattern1): j = %d, D[1] = 0x%x\n",j,D[1]);
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

    free(R[0]);
    free(R[1]);

    return count;
}
