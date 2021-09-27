#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <sys/resource.h>

#include "globals.h"
#include "functions.h"
#include "protein_table.h"

#include "bndm_aa.h"

int bndm_eds_aa_run(const unsigned char *pattern0,
                    const unsigned char *pattern1,
                    size_t m,
                    const int loops) {
    struct rusage ruse, ruse1;
    double ssec1, ssec2, usec1, usec2;
    int matches;

    if (m==0) {
        fprintf(stderr, "BNDM-EDS-AA requires a pattern!\n");
        return -1;
    }
    int pattern1_len = 0;
    if (pattern1) {
        pattern1_len = strnlen((const char *) pattern1, MAX_PATTERN_LENGTH + 1);
    }
    if (!pattern1_len && m > MAX_AA_PATTERN_LENGTH) {
        fprintf(stderr, "BNDM-EDS-AA requires AA pattern of maximum length %lu!\n", MAX_AA_PATTERN_LENGTH);
        return -1;
    }
    if (pattern1_len && m > MAX_PATTERN_LENGTH) {
        fprintf(stderr, "BNDM-EDS-AA requires IUPAC pattern of maximum length %lu!\n", MAX_PATTERN_LENGTH);
        return -1;
    }

    unsigned char IUPAC_patterns[2][MAX_PATTERN_LENGTH];
    if (pattern1_len) {
        printf("BNDM-EDS-AA (dual IUPAC pattern - supplied as arguments)\n");
        strncpy((unsigned char *) IUPAC_patterns[0], pattern0, m);
        strncpy((unsigned char *) IUPAC_patterns[1], pattern1, m);
    } else {
        if (translate_aa_iupac(pattern0, m, IUPAC_patterns)==0) {
            fprintf(stderr, "BNDM-EDS-AA failed to generate IUPAC patterns!\n");
            return -1;
        }
        printf("BNDM-EDS-AA (dual IUPAC pattern - generated from AA argument)\n");
        m = 3*m;
    }
    printf("Pattern0:\t\"%.*s\"\n", (int) m, IUPAC_patterns[0]);
    printf("Pattern1:\t\"%.*s\"\n", (int) m, IUPAC_patterns[1]);

    getrusage(RUSAGE_SELF, &ruse);
    ssec1 = (double) (ruse.ru_stime.tv_sec*1000000 + ruse.ru_stime.tv_usec);
    usec1 = (double) (ruse.ru_utime.tv_sec*1000000 + ruse.ru_utime.tv_usec);
    getrusage(RUSAGE_SELF, &ruse1);

    for (int i = 0; i < loops; i++) {
        aPointer = 0;
        matches = bndm_eds_iupac_search(writeBuffer, wbPointer, IUPAC_patterns[0], IUPAC_patterns[1], m);
    }

    getrusage(RUSAGE_SELF, &ruse);
    ssec2 = (double) (ruse.ru_stime.tv_sec*1000000 + ruse.ru_stime.tv_usec);
    usec2 = (double) (ruse.ru_utime.tv_sec*1000000 + ruse.ru_utime.tv_usec);

    printf("Matches:\t%d\n", matches);
    printf("User time:\t%f s\n", (usec2 - usec1)/(double) 1000000);
    printf("System time:\t%f s\n", (ssec2 - ssec1)/(double) 1000000);
    printf("Total time:\t%f s\n", ((usec2 + ssec2) - (usec1 + ssec1))/(double) 1000000);

    return matches;
}

int bndm_eds_aa_search(const unsigned char *text,
                       const size_t len,
                       const unsigned char *pattern,
                       size_t m) {
    assert(m <= MAX_AA_PATTERN_LENGTH);

    DEBUG_PRINT("BNDM-EDS-AA len=%lu, pattern=%.*s, m=%lu\n", len, (int) m, pattern, m);

    unsigned char IUPAC_patterns[2][MAX_PATTERN_LENGTH];
    if (translate_aa_iupac(pattern, m, IUPAC_patterns)==0) {
        fprintf(stderr, "BNDM-EDS-AA failed to generate IUPAC patterns!\n");
        return -1;
    }

    return bndm_eds_iupac_search(text, len, IUPAC_patterns[0], IUPAC_patterns[1], m*3);
}

int bndm_eds_iupac_search(const unsigned char *text,
                          const size_t len,
                          const unsigned char *pattern0,
                          const unsigned char *pattern1,
                          const size_t m) {
    assert(m%3==0);
    assert(m <= MAX_PATTERN_LENGTH);
    assert(strnlen(pattern0, m)==strnlen(pattern1, m));

    DEBUG_PRINT("  BNDM-EDS-AA IUPAC Pattern 0: %.*s\n", (int) m, pattern0);
    DEBUG_PRINT("                    Pattern 1: %.*s\n", (int) m, pattern1);

    int S0[SIGMA], S1[SIGMA], B0[SIGMA],
        B1[SIGMA]; // Preprocessed matching vectors for both patterns, S for ShiftAND, B for BNDM
    // int R[MAX_PATTERN_LENGTH];
    int F; // Preprocessed masking vector to verify full matches or prefix matches from matching vectors D[*]

    const unsigned int STATE_VECTOR_MERGE_MASK = 0x24924924U; // ShiftAnd Dual-Pattern Mask
    // const unsigned long int STATE_VECTOR_MERGE_MASK = 0x4924924924924924UL; // For 64bit matching

    int D0, D1; // Matching vectors used for both SA and BNDM
    int D2, D3; // Matching vector used to save state from BNDM prefix match to be used later in SA matching
    int DC0, DC1; // Matching states taken from the other pattern (masked by STATE_VECTOR_MERGE_MASK for every 3 chars)
    int R10, R11, R20, R21; // Matching vectors used to store state between segments

    /* BNDM Preprocessing */
    for (int i = 0; i < SIGMA; i++) {
        B0[i] = 0;
        B1[i] = 0;
    }
    F = 1;
    for (int i = m - 1; i >= 0; i--) {

        for (int b = 0; b < BASES; b++) {
            B0[IUPAC_SYMBOLS_TO_BASES[pattern0[i]][b]] |= F;
            B1[IUPAC_SYMBOLS_TO_BASES[pattern1[i]][b]] |= F;
        }
        // R[i] = F;
        F <<= 1;
    }
    F >>= 1;

    /* SA Preprocessing */
    for (int i = 0; i < SIGMA; ++i) {
        S0[i] = 0;
        S1[i] = 0;
    }
    for (int i = 0, j = 1; i < m; ++i, j <<= 1) {
        for (int b = 0; b < BASES; b++) {
            S0[IUPAC_SYMBOLS_TO_BASES[pattern0[i]][b]] |= j;
            S1[IUPAC_SYMBOLS_TO_BASES[pattern1[i]][b]] |= j;
        }
    }

    /* Searching */
    R10 = R11 = R20 = R21 = D2 = D3 = 0;
    unsigned int matches = 0, segmentCounter = 0, elementCounter = 0;

    while (aPointer < len) {
        // Processing the beginning of a segment
        unsigned char segmentSize = (unsigned char) text[aPointer++];
        segmentCounter++;
        elementCounter += segmentSize;
        R10 = R20;
        R11 = R21;
        R20 = 0;
        R21 = 0;
        DEBUG_PRINT(" SEGMENT segment=%d, elements=%d, totalElements=%d, aPointer=%d\n",
                    segmentCounter, segmentSize, elementCounter, aPointer);

        // Process one element of the segment.
        for (unsigned char k = 0; k < segmentSize; k++) {
            unsigned int elementLength = byteDecodeInt();
            unsigned int elementStart = aPointer;
            unsigned int elementEnd = elementStart + elementLength;
            DEBUG_PRINT(
                "  ELEMENT segment = %d, element = %d/%d, start = %d, end = %d, len = %d, elementCounter = %d\n",
                segmentCounter, k, segmentSize, elementStart, elementEnd, elementLength, elementCounter);

            // Perform SA search at the beginning of the element.
            D0 = R10;
            D1 = R11;
            for (int j = elementStart; j < (elementStart + m) && j < elementEnd; j++) {
                char curr_symbol = text[j];
                DC0 = D1 & STATE_VECTOR_MERGE_MASK;
                DC1 = D0 & STATE_VECTOR_MERGE_MASK;
                D0 = (((D0 | DC0) << 1) | 1) & S0[curr_symbol];
                D1 = (((D1 | DC1) << 1) | 1) & S1[curr_symbol];
                DEBUG_PRINT("    SA1 j=%d, curr_symbol=%c, D0=0x%x, D1=0x%x\n", j, curr_symbol, D0, D1);

                if (D0 & F | D1 & F) {
                    matches++;
                    DEBUG_PRINT("      MATCH (p0): j=%d, c=%c, D0=0x%x, S=0x%x, elementStart=%d, m=%lu, R1=0x%x\n",
                                j, curr_symbol, D0, S0[curr_symbol], elementStart, m, R10);
                    DEBUG_PRINT("      MATCH (p1): j=%d, c=%c, D1=0x%x, S=0x%x, elementStart=%d, m=%lu, R1=0x%x\n",
                                j, curr_symbol, D1, S1[curr_symbol], elementStart, m, R11);
                }
            }
            if (elementLength > m) {
                // Perform BNDM search.
                int j, last;
                for (j = elementStart + 1; j + m - 1 < elementEnd; j += last) {
                    D2 = D3 = 0;
                    D0 = D1 = ~0;
                    last = m;
                    DEBUG_PRINT("    BNDM j=%d, D2=0x%x, D3=0x%x, last=%d, text=%.*s\n",
                                j, D2, D3, last, (int) m, text + j);
                    for (int i = m - 1; i >= 0 && (D0!=0 || D1!=0); --i) {
                        char curr_symbol = text[j + i];
                        D0 = D0 & B0[curr_symbol];
                        D1 = D1 & B1[curr_symbol];

                        DEBUG_PRINT("      BNDM j=%d, i=%d, j+i=%d, curr_symbol=%c, (m-1-i)mod3=%lu, "
                                    "D0=0x%x, DC0=0x%x,B0[curr_symbol]=0x%x, "
                                    "D1=0x%x, DC1=0x%x, B1[curr_symbol]=0x%x, D1 | DC1 = %x\n",
                                    j, i, j + i, curr_symbol, (m - 1 - i)%3, D0, DC0, B0[curr_symbol], D1,
                                    DC1,
                                    B1[curr_symbol], D1 | DC1);

                        if ((D0 & F)!=0 && i > 0) {
                            last = i;
                            D2 |= 1 << (m - 1 - i); // D2 |= 1 << R[last];
                            DEBUG_PRINT("        BNDM Prefix: j=%d, i=%d, curr_symbol=%c, last=%d, D2=0x%x\n",
                                        j, i, curr_symbol, last, D2);
                        }
                        if ((D1 & F)!=0 && i > 0) {
                            last = i;
                            D3 |= 1 << (m - 1 - i); // D3 |= 1 << R[last];
                            DEBUG_PRINT("        BNDM Prefix: j=%d, i=%d, curr_symbol=%c, last=%d, D3=0x%x\n",
                                        j, i, curr_symbol, last, D3);
                        }
                        if (((D0 | D1) & F)!=0 && i==0) {
                            matches++;
                            DEBUG_PRINT("         MATCH (p0): j=%d, i=%d, c=%c, D=0x%x, B=0x%x\n",
                                        j, i, curr_symbol, D0, B0[curr_symbol]);
                            DEBUG_PRINT("         MATCH (p1): j=%d, i=%d, c=%c, D=0x%x, B=0x%x\n",
                                        j, i, curr_symbol, D1, B1[curr_symbol]);
                        }
                        DC0 = D1 & STATE_VECTOR_MERGE_MASK;
                        DC1 = D0 & STATE_VECTOR_MERGE_MASK;
                        D0 = (D0 | DC0) << 1;
                        D1 = (D1 | DC1) << 1;
                    }
                    DEBUG_PRINT("    BNDM last=%d, m=%lu\n", last, m);
                }

                // Perform SA search at the end of the element.
                D0 = D2; // Setting the initial value to SA register.
                D1 = D3; // Setting the initial value to SA register.
                DEBUG_PRINT("  BNDM->SA2 j=%d, D2=0x%x, D3=0x%x\n", j, D2, D3);

                for (j = j + m - last; j < elementEnd; j++) {
                    char curr_symbol = text[j];
                    DC0 = D1 & STATE_VECTOR_MERGE_MASK;
                    DC1 = D0 & STATE_VECTOR_MERGE_MASK;
                    D0 = (((D0 | DC0) << 1) | 1) & S0[curr_symbol];
                    D1 = (((D1 | DC1) << 1) | 1) & S1[curr_symbol];
                    DEBUG_PRINT("    SA2 j=%d, curr_symbol=%c, D0=0x%x, D1=0x%x\n", j, curr_symbol, D0, D1);
                }
            }

            // Merge the SA registers of individual elements into a cumulative register for the whole segment
            R20 |= D0;
            R21 |= D1;

            // Set the pointer to the beginning of the next element.
            aPointer += elementLength;
        }
    }

    DEBUG_PRINT("BNDM-EDS-AA Matches=%d, segments=%d, elements=%d\n", matches, segmentCounter, elementCounter);

    return matches;
}
