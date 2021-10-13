#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <sys/resource.h>

#include "globals.h"
#include "functions.h"
#include "protein_table.h"

#include "bndm_eds_mp.h"

// If USE_STATIC_NUM_PATTERNS is defined, all the pattern loops will use a constant number of patterns set to
// MAX_DNA_PATTERNS, even if the actual runtime number of patterns is lower. This is to allow the compiler optimize
// the code using loop unrolling, which is not possible for runtime-dependent number of patterns.

#define USE_STATIC_NUM_PATTERNS

#ifdef USE_STATIC_NUM_PATTERNS
#define NUM_PATTERNS (MAX_DNA_PATTERNS)
#else
#define NUM_PATTERNS (num_patterns)
#endif

int bndm_eds_mp_run(const unsigned char *teds,
                    const size_t len,
                    const unsigned char *pattern,
                    const size_t m,
                    const int loops) {
    struct rusage ruse, ruse1, ruse2;
    double ssec1, ssec2, usec1, usec2;
    int matches;

    if (m==0) {
        fprintf(stderr, "BNDM-EDS-MP requires for AA pattern to be specified as argument!\n");
        return 1;
    }
    unsigned char IUPAC_patterns[MAX_DNA_PATTERNS][MAX_PATTERN_LENGTH];
    size_t num_iupac_patterns = translate_aa_iupac_all_combinations(pattern, m, IUPAC_patterns,
                                                                    MAX_DNA_PATTERNS, MAX_PATTERN_LENGTH);
    if (num_iupac_patterns==0) {
        fprintf(stderr, "BNDM-EDS-MP failed to generate any DNA patterns!\n");
        return 1;
    }
    for (size_t i = num_iupac_patterns; i < MAX_DNA_PATTERNS; i++) {
        memset(IUPAC_patterns[i], 0, sizeof(unsigned char)*MAX_PATTERN_LENGTH);
    }

    printf("BNDM-EDS-MP (Multi-Patterns generated from AA)\n");
    for (size_t i = 0; i < num_iupac_patterns; i++) {
        printf("Pattern%ld:\t\"%.*s\"\n", i, (int) m*3, IUPAC_patterns[i]);
    }

    getrusage(RUSAGE_SELF, &ruse);
    ssec1 = (double) (ruse.ru_stime.tv_sec*1000000 + ruse.ru_stime.tv_usec);
    usec1 = (double) (ruse.ru_utime.tv_sec*1000000 + ruse.ru_utime.tv_usec);
    getrusage(RUSAGE_SELF, &ruse1);

    for (int i = 0; i < loops; i++) {
        matches = bndm_eds_mp_search(teds, len, IUPAC_patterns, num_iupac_patterns, m*3);
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

int bndm_eds_mp_search(const unsigned char *text,
                       const size_t len,
                       unsigned char patterns[MAX_DNA_PATTERNS][MAX_PATTERN_LENGTH],
                       unsigned int num_patterns,
                       unsigned int m) {
    assert(m%3==0);
    assert(m <= MAX_PATTERN_LENGTH);
    assert(num_patterns <= MAX_DNA_PATTERNS);

    DEBUG_PRINT("  BNDM-EDS-AA IUPAC Pattern 0: %.*s\n", (int) m, pattern0);
    DEBUG_PRINT("                    Pattern 1: %.*s\n", (int) m, pattern1);

    int S[MAX_DNA_PATTERNS][SIGMA],
        B[MAX_DNA_PATTERNS][SIGMA]; // Preprocessed matching vectors for both patterns, S for ShiftAND, B for BNDM
    // int R[MAX_PATTERN_LENGTH];
    int F; // Preprocessed masking vector to verify full matches or prefix matches from matching vectors D[*]

    const unsigned int STATE_VECTOR_MERGE_MASK = 0x24924924U; // ShiftAnd Dual-Pattern Mask
    // const unsigned long int STATE_VECTOR_MERGE_MASK = 0x4924924924924924UL; // For 64bit matching

    int D[MAX_DNA_PATTERNS]; // Matching vectors used for both SA and BNDM
    int D2[MAX_DNA_PATTERNS
    ]; // Matching vector used to save state from BNDM prefix match to be used later in SA matching
    int R1[MAX_DNA_PATTERNS], R2[MAX_DNA_PATTERNS]; // Matching vectors used to store state between segments

    /* BNDM Preprocessing */
    for (int i = 0; i < SIGMA; i++) {
        for (unsigned int p = 0; p < NUM_PATTERNS; p++) {
            B[p][i] = 0;
        }
    }
    F = 1;
    for (int i = m - 1; i >= 0; i--) {
        for (int b = 0; b < BASES; b++) {
            for (unsigned int p = 0; p < NUM_PATTERNS; p++) {
                B[p][IUPAC_SYMBOLS_TO_BASES[patterns[p][i]][b]] |= F;
            }
        }
        F <<= 1;
    }
    F >>= 1;

    /* SA Preprocessing */
    for (int i = 0; i < SIGMA; ++i) {
        for (unsigned int p = 0; p < NUM_PATTERNS; p++) {
            S[p][i] = 0;
        }
    }
    for (int i = 0, j = 1; i < m; ++i, j <<= 1) {
        for (int b = 0; b < BASES; b++) {
            for (unsigned int p = 0; p < NUM_PATTERNS; p++) {
                S[p][IUPAC_SYMBOLS_TO_BASES[patterns[p][i]][b]] |= j;
            }
        }
    }

    /* Searching */
    for (unsigned int p = 0; p < NUM_PATTERNS; p++) {
        R1[p] = R2[p] = D2[p] = 0;
    }
    unsigned int matches = 0, segmentCounter = 0, elementCounter = 0;

    int aPointer = 0;
    while (aPointer < len) {
        // Processing the beginning of a segment
        unsigned int segmentSize;
        aPointer += byteDecodeInt(text + aPointer, &segmentSize);
        segmentCounter++;
        elementCounter += segmentSize;
        for (unsigned int p = 0; p < NUM_PATTERNS; p++) {
            R1[p] = R2[p];
            R2[p] = 0;
        }
        DEBUG_PRINT(" SEGMENT segment=%d, elements=%d, totalElements=%d, aPointer=%d\n",
                    segmentCounter, segmentSize, elementCounter, aPointer);

        // Process one element of the segment.
        for (unsigned int k = 0; k < segmentSize; k++) {
            unsigned int elementLength;
            aPointer += byteDecodeInt(text + aPointer, &elementLength);
            unsigned int elementStart = aPointer;
            unsigned int elementEnd = elementStart + elementLength;
            DEBUG_PRINT(
                "  ELEMENT segment = %d, element = %d/%d, start = %d, end = %d, len = %d, elementCounter = %d\n",
                segmentCounter, k, segmentSize, elementStart, elementEnd, elementLength, elementCounter);

            // Perform SA search at the beginning of the element.
            for (unsigned int p = 0; p < NUM_PATTERNS; p++) {
                D[p] = R1[p];
            }
            for (int j = elementStart; j < (elementStart + m) && j < elementEnd; j++) {
                char curr_symbol = text[j];
                for (unsigned int p = 0; p < NUM_PATTERNS; p++) {
                    D[p] = ((D[p] << 1) | 1) & S[p][curr_symbol];
                }
                DEBUG_PRINT("    SA1 j=%d, curr_symbol=%c, D0=0x%x, D1=0x%x\n", j, curr_symbol, D0, D1);

                int match = 0;
                for (unsigned int p = 0; p < NUM_PATTERNS; p++) {
                    match |= (D[p] & F);
                }
                if (match) {
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
                    for (unsigned int p = 0; p < NUM_PATTERNS; p++) {
                        D2[p] = 0;
                        D[p] = ~0;
                    }
                    int D_OR = ~0;
                    last = m;
                    DEBUG_PRINT("    BNDM j=%d, D2=0x%x, D3=0x%x, last=%d, text=%.*s\n",
                                j, D2, D3, last, (int) m, text + j);
                    for (int i = m - 1; i >= 0 && D_OR; --i) {
                        char curr_symbol = text[j + i];
                        for (unsigned int p = 0; p < NUM_PATTERNS; p++) {
                            D[p] = D[p] & B[p][curr_symbol];
                        }

                        DEBUG_PRINT("      BNDM j=%d, i=%d, j+i=%d, curr_symbol=%c, (m-1-i)mod3=%lu, "
                                    "D0=0x%x, DC0=0x%x,B0[curr_symbol]=0x%x, "
                                    "D1=0x%x, DC1=0x%x, B1[curr_symbol]=0x%x, D1 | DC1 = %x\n",
                                    j, i, j + i, curr_symbol, (m - 1 - i)%3, D0, DC0, B0[curr_symbol], D1,
                                    DC1,
                                    B1[curr_symbol], D1 | DC1);

                        for (unsigned int p = 0; p < NUM_PATTERNS; p++) {
                            if ((D[p] & F)!=0 && i > 0) {
                                last = i;
                                D2[p] |= 1 << (m - 1 - i); // D2 |= 1 << R[last];
                                DEBUG_PRINT("        BNDM Prefix: j=%d, i=%d, curr_symbol=%c, last=%d, D2=0x%x\n",
                                            j, i, curr_symbol, last, D2);
                            }
                        }
                        int D_OR = 0;
                        for (unsigned int p = 0; p < NUM_PATTERNS; p++) {
                            D_OR |= D[p];
                        }
                        if ((D_OR & F)!=0 && i==0) {
                            matches++;
                            DEBUG_PRINT("         MATCH (p0): j=%d, i=%d, c=%c, D=0x%x, B=0x%x\n",
                                        j, i, curr_symbol, D0, B0[curr_symbol]);
                            DEBUG_PRINT("         MATCH (p1): j=%d, i=%d, c=%c, D=0x%x, B=0x%x\n",
                                        j, i, curr_symbol, D1, B1[curr_symbol]);
                        }
                        for (unsigned int p = 0; p < NUM_PATTERNS; p++) {
                            D[p] = D[p] << 1;
                        }
                        D_OR = D_OR << 1;
                    }
                    DEBUG_PRINT("    BNDM last=%d, m=%lu\n", last, m);
                }

                // Perform SA search at the end of the element.
                for (unsigned int p = 0; p < NUM_PATTERNS; p++) {
                    D[p] = D2[p]; // Setting the initial value to SA register.
                }
                DEBUG_PRINT("  BNDM->SA2 j=%d, D2=0x%x, D3=0x%x\n", j, D2, D3);

                for (j = j + m - last; j < elementEnd; j++) {
                    char curr_symbol = text[j];
                    for (unsigned int p = 0; p < NUM_PATTERNS; p++) {
                        D[p] = ((D[p] << 1) | 1) & S[p][curr_symbol];
                    }
                    DEBUG_PRINT("    SA2 j=%d, curr_symbol=%c, D0=0x%x, D1=0x%x\n", j, curr_symbol, D0, D1);
                }
            }

            // Merge the SA registers of individual elements into a cumulative register for the whole segment
            for (unsigned int p = 0; p < NUM_PATTERNS; p++) {
                R2[p] |= D[p];
            }

            // Set the pointer to the beginning of the next element.
            aPointer += elementLength;
        }
    }

    DEBUG_PRINT("BNDM-EDS-AA Matches=%d, segments=%d, elements=%d\n", matches, segmentCounter, elementCounter);

    return matches;
}
