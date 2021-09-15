#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "bndm_aa.h"
#include "globals.h"
#include "functions.h"
#include "protein_table.h"


int bndm_eds_aa_search(const unsigned char* text,
                       const size_t len,
                       const unsigned char *pattern,
                       const size_t m)
{
    assert(m <= MAX_AA_PATTERN_LENGTH);

    DEBUG_PRINT("BNDM-EDS-AA len=%lu, pattern=%.*s, m=%lu\n", len, (int)m, pattern, m);

    unsigned char patterns[2][MAX_PATTERN_LENGTH];
    if (translate_aa_iupac(pattern, m, patterns) == 0){
        fprintf(stderr, "BNDM-EDS-AA failed to generate any DNA patterns!\n");
        return -1;
    }

    DEBUG_PRINT("  BNDM-EDS-AA Pattern0: %s\n", patterns[0]);
    DEBUG_PRINT("  BNDM-EDS-AA Pattern1: %s\n", patterns[1]);

    return bndm_eds_iupac_search(text, len, patterns[0], patterns[1], m*3);
}

int bndm_eds_iupac_search(const unsigned char* text,
                          const size_t len,
                          const unsigned char *pattern0,
                          const unsigned char *pattern1,
                          const size_t m)
{
    assert(m % 3 == 0);
    assert(m <= MAX_PATTERN_LENGTH);
    assert(strnlen(pattern0, m) == strnlen(pattern1, m));

    int S[2][SIGMA], B[2][SIGMA]; // Preprocessed matching vectors for both patterns, S for ShiftAND, B for BNDM
    int F; // Preprocessed masking vector to verify full matches or prefix matches from matching vectors D[*]

    const unsigned int STATE_VECTOR_MERGE_MASK = 0x24924924U; // ShiftAnd Dual-Pattern Mask
    // const unsigned long int STATE_VECTOR_MERGE_MASK = 0x4924924924924924UL; // For 64bit matching

    int D[2]; // Matching vectors used for both SA and BNDM
    int D2[2]; // Matching vector used to save state from BNDM prefix match to be used later in SA matching
    int DC[2]; // Matching states taken from the other pattern (masked by STATE_VECTOR_MERGE_MASK for every 3 chars)
    int R1[2], R2[2]; // Matching vectors used to store state between segments

    /* BNDM Preprocessing */
    for (int i = 0; i<SIGMA; i++){
        B[0][i] = 0;
        B[1][i] = 0;
    }
    F = 1;
    for (int i = m - 1; i >= 0; i--) {

        for (int b = 0; b < BASES; b++){
            B[0][IUPAC_SYMBOLS_TO_BASES[pattern0[i]][b]] |= F;
            B[1][IUPAC_SYMBOLS_TO_BASES[pattern1[i]][b]] |= F;
        }
        F <<= 1;
    }
    F >>= 1;

    /* SA Preprocessing */
    for (int i = 0; i < SIGMA; ++i){
        S[0][i] = 0;
        S[1][i] = 0;
    }
    for (int i = 0, j = 1; i < m; ++i, j <<= 1){
        for (int b = 0; b < BASES; b++){
            S[0][IUPAC_SYMBOLS_TO_BASES[pattern0[i]][b]] |= j;
            S[1][IUPAC_SYMBOLS_TO_BASES[pattern1[i]][b]] |= j;
        }
    }

    /* Searching */
    R1[0] = R1[1] = R2[0] = R2[1] = D2[0] = D2[1] = 0;
    unsigned int matches = 0, segmentCounter = 0, elementCounter = 0;

    while (aPointer < len)
    {
        // Processing the beginning of a segment
        unsigned char segmentSize = (unsigned char)text[aPointer++];
        segmentCounter++;
        elementCounter += segmentSize;
        R1[0] = R2[0];
        R1[1] = R2[1];
        R2[0] = 0;
        R2[1] = 0;
        DEBUG_PRINT(" SEGMENT segment=%d, elements=%d, totalElements=%d, aPointer=%d\n",
                    segmentCounter, segmentSize, elementCounter, aPointer);

        // Process one element of the segment.
        for (unsigned char k = 0; k < segmentSize; k++)
        {
            unsigned int elementLength = byteDecodeInt();
            unsigned int elementStart = aPointer;
            unsigned int elementEnd = elementStart + elementLength;
            DEBUG_PRINT("  ELEMENT segment = %d, element = %d/%d, start = %d, end = %d, len = %d, elementCounter = %d\n",
                        segmentCounter, k, segmentSize, elementStart, elementEnd, elementLength, elementCounter);

            // Perform SA search at the beginning of the element.
            D[0] = R1[0];
            D[1] = R1[1];
            for (int j = elementStart; j < (elementStart + m) && j < elementEnd; j++) {
                char curr_symbol = text[j];
                DC[0] = D[1] & STATE_VECTOR_MERGE_MASK;
                DC[1] = D[0] & STATE_VECTOR_MERGE_MASK;
                D[0] = (((D[0] | DC[0]) << 1) | 1) & S[0][curr_symbol];
                D[1] = (((D[1] | DC[1]) << 1) | 1) & S[1][curr_symbol];
                DEBUG_PRINT("    SA1 j=%d, curr_symbol=%c, D[0]=0x%x, D[1]=0x%x\n", j, curr_symbol, D[0], D[1]);

                if (D[0] & F | D[1] & F) {
                    matches++;
                    DEBUG_PRINT("      MATCH (p0): j=%d, c=%c, D[0]=0x%x, S=0x%x, elementStart=%d, m=%lu, R1=0x%x\n",
                                j, curr_symbol, D[0], S[0][curr_symbol], elementStart, m, R1[0]);
                    DEBUG_PRINT("      MATCH (p1): j=%d, c=%c, D[1]=0x%x, S=0x%x, elementStart=%d, m=%lu, R1=0x%x\n",
                                j, curr_symbol, D[1], S[1][curr_symbol], elementStart, m, R1[1]);
                }
            }
            if (elementLength <= m) {
                R2[0] |= D[0];
                R2[1] |= D[1];
            }

            // Perform BNDM search.
            int j, last;
            for (j = elementStart + 1; j + m - 1 < elementEnd; j += last)
            {
                D2[0] = D2[1] = 0;
                D[0] = D[1] = ~0;
                last = m;
                DEBUG_PRINT("    BNDM j=%d, D2[0]=0x%x, D2[1]=0x%x, last=%d, text=%.*s\n",
                            j, D2[0], D2[1], last, (int)m, text+j);
                for (int i = m -1; i >= 0 && (D[0] != 0 || D[1] != 0); --i) {
                    char curr_symbol = text[j + i];
                    D[0] = D[0] & B[0][curr_symbol];
                    D[1] = D[1] & B[1][curr_symbol];

                    DEBUG_PRINT("      BNDM j=%d, i=%d, j+i=%d, curr_symbol=%c, (m-1-i)mod3=%lu, "
                                "D[0]=0x%x, DC[0]=0x%x,B[0][curr_symbol]=0x%x, "
                                "D[1]=0x%x, DC[1]=0x%x, B[1][curr_symbol]=0x%x, D[1] | DC[1] = %x\n",
                                j, i, j+i,curr_symbol, (m - 1 - i) % 3, D[0], DC[0], B[0][curr_symbol], D[1], DC[1],
                                B[1][curr_symbol], D[1] | DC[1]);

                    if ((D[0] & F) != 0 && i > 0) {
                        last = i;
                        D2[0] |= 1 << (m - 1 - i);
                        DEBUG_PRINT("        BNDM Prefix: j=%d, i=%d, curr_symbol=%c, last=%d, D2[0]=0x%x\n",
                                    j, i, curr_symbol, last, D2[0]);
                    }
                    if ((D[1] & F) != 0 && i > 0) {
                        last = i;
                        D2[1] |= 1 << (m - 1 - i);
                        DEBUG_PRINT("        BNDM Prefix: j=%d, i=%d, curr_symbol=%c, last=%d, D2[1]=0x%x\n",
                                    j, i, curr_symbol, last, D2[1]);
                    }
                    if (((D[0] | D[1]) & F) != 0 && i == 0){
                        matches++;
                        DEBUG_PRINT("         MATCH (p0): j=%d, i=%d, c=%c, D=0x%x, B=0x%x\n",
                                    j, i, curr_symbol, D[0], B[0][curr_symbol]);
                        DEBUG_PRINT("         MATCH (p1): j=%d, i=%d, c=%c, D=0x%x, B=0x%x\n",
                                    j, i, curr_symbol, D[1], B[1][curr_symbol]);
                    }
                    DC[0] = D[1] & STATE_VECTOR_MERGE_MASK;
                    DC[1] = D[0] & STATE_VECTOR_MERGE_MASK;
                    D[0] = (D[0] | DC[0]) << 1;
                    D[1] = (D[1] | DC[1]) << 1;
                }
                DEBUG_PRINT("    BNDM last=%d, m=%lu\n", last, m);
            }

            // Perform SA search at the end of the element.
            D[0] = D2[0]; // Setting the initial value to SA register.
            D[1] = D2[1]; // Setting the initial value to SA register.
            DEBUG_PRINT("  BNDM->SA2 j=%d, D2[0]=0x%x, D2[1]=0x%x\n", j, D2[0], D2[1]);
            for(j = j + m - last; j < elementEnd; j++) {
                char curr_symbol = text[j];
                DC[0] = D[1] & STATE_VECTOR_MERGE_MASK;
                DC[1] = D[0] & STATE_VECTOR_MERGE_MASK;
                D[0] = (((D[0] | DC[0]) << 1) | 1) & S[0][curr_symbol];
                D[1] = (((D[1] | DC[1]) << 1) | 1) & S[1][curr_symbol];
                DEBUG_PRINT("    SA2 j=%d, curr_symbol=%c, D[0]=0x%x, D[1]=0x%x\n", j, curr_symbol, D[0], D[1]);
            }

            // Merge the SA registers of individual elements into a cumulative register for the whole segment
            R2[0] |= D[0];
            R2[1] |= D[1];

            // Set the pointer to the beginning of the next element.
            aPointer += elementLength;
        }
    }

    printf("BNDM-EDS-AA Matches=%d, segments=%d, elements=%d\n", matches, segmentCounter, elementCounter);

    return matches;
}
