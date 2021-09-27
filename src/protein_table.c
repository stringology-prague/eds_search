#include <string.h>
#include <stdio.h>

#include "protein_table.h"

unsigned char IUPAC_SYMBOLS_TO_BASES[SIGMA][BASES];

void init_IUPAC_SYMBOLS_TO_BASES() {
    memset(IUPAC_SYMBOLS_TO_BASES, 0, SIGMA*BASES*sizeof(unsigned char));

    IUPAC_SYMBOLS_TO_BASES['A'][0] = 'A';
    IUPAC_SYMBOLS_TO_BASES['C'][0] = 'C';
    IUPAC_SYMBOLS_TO_BASES['G'][0] = 'G';
    IUPAC_SYMBOLS_TO_BASES['T'][0] = 'T';
    IUPAC_SYMBOLS_TO_BASES['U'][0] = 'U';

    IUPAC_SYMBOLS_TO_BASES['W'][0] = 'A';
    IUPAC_SYMBOLS_TO_BASES['S'][0] = 'C';
    IUPAC_SYMBOLS_TO_BASES['M'][0] = 'A';
    IUPAC_SYMBOLS_TO_BASES['K'][0] = 'G';
    IUPAC_SYMBOLS_TO_BASES['R'][0] = 'A';
    IUPAC_SYMBOLS_TO_BASES['Y'][0] = 'C';

    IUPAC_SYMBOLS_TO_BASES['W'][1] = 'T';
    IUPAC_SYMBOLS_TO_BASES['S'][1] = 'G';
    IUPAC_SYMBOLS_TO_BASES['M'][1] = 'C';
    IUPAC_SYMBOLS_TO_BASES['K'][1] = 'T';
    IUPAC_SYMBOLS_TO_BASES['R'][1] = 'G';
    IUPAC_SYMBOLS_TO_BASES['Y'][1] = 'T';

    IUPAC_SYMBOLS_TO_BASES['B'][0] = 'C';
    IUPAC_SYMBOLS_TO_BASES['D'][0] = 'A';
    IUPAC_SYMBOLS_TO_BASES['H'][0] = 'A';
    IUPAC_SYMBOLS_TO_BASES['V'][0] = 'A';

    IUPAC_SYMBOLS_TO_BASES['B'][1] = 'G';
    IUPAC_SYMBOLS_TO_BASES['D'][1] = 'G';
    IUPAC_SYMBOLS_TO_BASES['H'][1] = 'C';
    IUPAC_SYMBOLS_TO_BASES['V'][1] = 'C';

    IUPAC_SYMBOLS_TO_BASES['B'][2] = 'T';
    IUPAC_SYMBOLS_TO_BASES['D'][2] = 'T';
    IUPAC_SYMBOLS_TO_BASES['H'][2] = 'T';
    IUPAC_SYMBOLS_TO_BASES['V'][2] = 'G';

    IUPAC_SYMBOLS_TO_BASES['N'][0] = 'A';
    IUPAC_SYMBOLS_TO_BASES['N'][1] = 'C';
    IUPAC_SYMBOLS_TO_BASES['N'][2] = 'G';
    IUPAC_SYMBOLS_TO_BASES['N'][3] = 'T';

    return;
}

unsigned char AA_TO_COMPR_IUPAC_SYMBOLS[SIGMA][COMPR_IUPAC_MAX_VARIATIONS][COMPR_IUPAC_LEN];

void init_AA_TO_COMPR_IUPAC_SYMBOLS() {
    memset(AA_TO_COMPR_IUPAC_SYMBOLS, 0,
           SIGMA*COMPR_IUPAC_MAX_VARIATIONS*COMPR_IUPAC_LEN*sizeof(unsigned char));

    strncpy((char *) AA_TO_COMPR_IUPAC_SYMBOLS['A'][0], "GCN", COMPR_IUPAC_LEN);
    strncpy((char *) AA_TO_COMPR_IUPAC_SYMBOLS['I'][0], "ATH", COMPR_IUPAC_LEN);
    strncpy((char *) AA_TO_COMPR_IUPAC_SYMBOLS['R'][0], "CGN", COMPR_IUPAC_LEN);
    strncpy((char *) AA_TO_COMPR_IUPAC_SYMBOLS['R'][1], "AGR", COMPR_IUPAC_LEN);
//    strncpy((char *) AA_TO_COMPR_IUPAC_SYMBOLS['R'][0], "CGY", COMPR_IUPAC_LEN);
//    strncpy((char *) AA_TO_COMPR_IUPAC_SYMBOLS['R'][1], "MGR", COMPR_IUPAC_LEN);
    strncpy((char *) AA_TO_COMPR_IUPAC_SYMBOLS['L'][0], "CTN", COMPR_IUPAC_LEN);
    strncpy((char *) AA_TO_COMPR_IUPAC_SYMBOLS['L'][1], "TTR", COMPR_IUPAC_LEN);
//    strncpy((char *) AA_TO_COMPR_IUPAC_SYMBOLS['L'][0], "CTY", COMPR_IUPAC_LEN);
//    strncpy((char *) AA_TO_COMPR_IUPAC_SYMBOLS['L'][1], "YTR", COMPR_IUPAC_LEN);
    strncpy((char *) AA_TO_COMPR_IUPAC_SYMBOLS['N'][0], "AAY", COMPR_IUPAC_LEN);
    strncpy((char *) AA_TO_COMPR_IUPAC_SYMBOLS['K'][0], "AAR", COMPR_IUPAC_LEN);
    strncpy((char *) AA_TO_COMPR_IUPAC_SYMBOLS['D'][0], "GAY", COMPR_IUPAC_LEN);
    strncpy((char *) AA_TO_COMPR_IUPAC_SYMBOLS['M'][0], "ATG", COMPR_IUPAC_LEN);
    strncpy((char *) AA_TO_COMPR_IUPAC_SYMBOLS['B'][0], "RAY", COMPR_IUPAC_LEN);
    strncpy((char *) AA_TO_COMPR_IUPAC_SYMBOLS['F'][0], "TTY", COMPR_IUPAC_LEN);
    strncpy((char *) AA_TO_COMPR_IUPAC_SYMBOLS['C'][0], "TGY", COMPR_IUPAC_LEN);
    strncpy((char *) AA_TO_COMPR_IUPAC_SYMBOLS['P'][0], "CCN", COMPR_IUPAC_LEN);
    strncpy((char *) AA_TO_COMPR_IUPAC_SYMBOLS['Q'][0], "CAR", COMPR_IUPAC_LEN);
    strncpy((char *) AA_TO_COMPR_IUPAC_SYMBOLS['S'][0], "TCN", COMPR_IUPAC_LEN);
    strncpy((char *) AA_TO_COMPR_IUPAC_SYMBOLS['S'][1], "AGY", COMPR_IUPAC_LEN);
    strncpy((char *) AA_TO_COMPR_IUPAC_SYMBOLS['E'][0], "GAR", COMPR_IUPAC_LEN);
    strncpy((char *) AA_TO_COMPR_IUPAC_SYMBOLS['T'][0], "ACN", COMPR_IUPAC_LEN);
    strncpy((char *) AA_TO_COMPR_IUPAC_SYMBOLS['Z'][0], "SAR", COMPR_IUPAC_LEN);
    strncpy((char *) AA_TO_COMPR_IUPAC_SYMBOLS['W'][0], "TGG", COMPR_IUPAC_LEN);
    strncpy((char *) AA_TO_COMPR_IUPAC_SYMBOLS['G'][0], "GGN", COMPR_IUPAC_LEN);
    strncpy((char *) AA_TO_COMPR_IUPAC_SYMBOLS['Y'][0], "TAY", COMPR_IUPAC_LEN);
    strncpy((char *) AA_TO_COMPR_IUPAC_SYMBOLS['H'][0], "CAY", COMPR_IUPAC_LEN);
    strncpy((char *) AA_TO_COMPR_IUPAC_SYMBOLS['V'][0], "GTN", COMPR_IUPAC_LEN);
    strncpy((char *) AA_TO_COMPR_IUPAC_SYMBOLS[AA_START][0], "ATG", COMPR_IUPAC_LEN);
    strncpy((char *) AA_TO_COMPR_IUPAC_SYMBOLS[AA_STOP][0], "TRA", COMPR_IUPAC_LEN);
    strncpy((char *) AA_TO_COMPR_IUPAC_SYMBOLS[AA_STOP][1], "TAR", COMPR_IUPAC_LEN);
}

size_t translate_aa_iupac_all_combinations(const unsigned char *aa_pattern,
                                           const size_t aa_pattern_size,
                                           unsigned char dna_patterns[][MAX_PATTERN_LENGTH],
                                           const size_t max_dna_patterns,
                                           const size_t max_dna_pattern_length) {
    DEBUG_PRINT("aa_pattern: \"%1.*s\", aa_pattern_size=%lu, max_dna_patterns=%lu, max_dna_pattern_length=%lu\n",
           (int) aa_pattern_size,
           aa_pattern,
           aa_pattern_size,
           max_dna_patterns,
           max_dna_pattern_length);
    if (aa_pattern_size*3 > max_dna_pattern_length || aa_pattern_size==0 || aa_pattern==NULL) {
        return 0;
    }

    // Do a preliminary check on the number of DNA sequences to be generated - and fail if too many
    int combinations = 1;
    for (int i = 0; i < aa_pattern_size; i++) {
        if (AA_TO_COMPR_IUPAC_SYMBOLS[aa_pattern[i]][1][0]!=0) {
            combinations <<= 1;
        }
    }
    if (combinations > max_dna_patterns) {
        fprintf(stderr,
                "ERROR: AA pattern \"%1.*s\" translated to %d DNA patterns, max allowed combinations are %lu!\n",
                (int) aa_pattern_size,
                aa_pattern,
                combinations,
                max_dna_patterns);
        return 0;
    }

    for (int i = 0; i < max_dna_patterns; i++) {
        memset(dna_patterns[i], 0, max_dna_pattern_length);
    }

    int active_dna_patterns = 1;
    for (int i = 0; i < aa_pattern_size; i++) {
        // Check if this is a symbol that codes into multiple disjoint DNA sequences
        if (AA_TO_COMPR_IUPAC_SYMBOLS[aa_pattern[i]][1][0]!=0) {
            // Duplicate current DNA patterns and append next 3 DNA symbols
            for (int p = 0; p < active_dna_patterns; p++) {
                memcpy(dna_patterns[p + active_dna_patterns], dna_patterns[p], 3*i);
                dna_patterns[p][3*i + 0] = AA_TO_COMPR_IUPAC_SYMBOLS[aa_pattern[i]][0][0];
                dna_patterns[p][3*i + 1] = AA_TO_COMPR_IUPAC_SYMBOLS[aa_pattern[i]][0][1];
                dna_patterns[p][3*i + 2] = AA_TO_COMPR_IUPAC_SYMBOLS[aa_pattern[i]][0][2];
                dna_patterns[p + active_dna_patterns][3*i + 0] = AA_TO_COMPR_IUPAC_SYMBOLS[aa_pattern[i]][1][0];
                dna_patterns[p + active_dna_patterns][3*i + 1] = AA_TO_COMPR_IUPAC_SYMBOLS[aa_pattern[i]][1][1];
                dna_patterns[p + active_dna_patterns][3*i + 2] = AA_TO_COMPR_IUPAC_SYMBOLS[aa_pattern[i]][1][2];
            }
            active_dna_patterns <<= 1;
        } else {
            // Append next 3 DNA symbols
            for (int p = 0; p < active_dna_patterns; p++) {
                dna_patterns[p][3*i + 0] = AA_TO_COMPR_IUPAC_SYMBOLS[aa_pattern[i]][0][0];
                dna_patterns[p][3*i + 1] = AA_TO_COMPR_IUPAC_SYMBOLS[aa_pattern[i]][0][1];
                dna_patterns[p][3*i + 2] = AA_TO_COMPR_IUPAC_SYMBOLS[aa_pattern[i]][0][2];
            }
        }
    }
    for (int p = 0; p < active_dna_patterns; p++) {
        DEBUG_PRINT("dna_patterns[%d]: %1.*s\n", p, (int) aa_pattern_size*3, dna_patterns[p]);
    }
    return active_dna_patterns;
}

size_t translate_aa_iupac(const unsigned char *aa_pattern,
                          const size_t aa_pattern_size,
                          unsigned char dna_patterns[2][MAX_PATTERN_LENGTH]) {
    DEBUG_PRINT("TRANSLATE_AA_IUPAC aa_pattern: \"%1.*s\", aa_pattern_size=%lu\n",
                (int) aa_pattern_size, aa_pattern, aa_pattern_size);
    if (aa_pattern_size > MAX_AA_PATTERN_LENGTH || aa_pattern_size==0 || aa_pattern==NULL) {
        return 0;
    }

    memset(dna_patterns[0], 0, MAX_PATTERN_LENGTH);
    memset(dna_patterns[1], 0, MAX_PATTERN_LENGTH);

    for (int i = 0; i < aa_pattern_size; i++) {
        // Expand AA symbol to IUPAC into the first pattern
        dna_patterns[0][3*i + 0] = AA_TO_COMPR_IUPAC_SYMBOLS[aa_pattern[i]][0][0];
        dna_patterns[0][3*i + 1] = AA_TO_COMPR_IUPAC_SYMBOLS[aa_pattern[i]][0][1];
        dna_patterns[0][3*i + 2] = AA_TO_COMPR_IUPAC_SYMBOLS[aa_pattern[i]][0][2];

        // Check if this is AA symbol that codes into multiple disjoint DNA sequences
        if (AA_TO_COMPR_IUPAC_SYMBOLS[aa_pattern[i]][1][0]!=0) {
            // Fill in second pattern with the second IUPAC codon
            dna_patterns[1][3*i + 0] = AA_TO_COMPR_IUPAC_SYMBOLS[aa_pattern[i]][1][0];
            dna_patterns[1][3*i + 1] = AA_TO_COMPR_IUPAC_SYMBOLS[aa_pattern[i]][1][1];
            dna_patterns[1][3*i + 2] = AA_TO_COMPR_IUPAC_SYMBOLS[aa_pattern[i]][1][2];

        } else {
            // Fill in second pattern with the only IUPAC codon
            dna_patterns[1][3*i + 0] = dna_patterns[0][3*i + 0];
            dna_patterns[1][3*i + 1] = dna_patterns[0][3*i + 1];
            dna_patterns[1][3*i + 2] = dna_patterns[0][3*i + 2];
        }
    }
    DEBUG_PRINT("  dna_patterns[0]: %1.*s\n", (int) aa_pattern_size*3, dna_patterns[0]);
    DEBUG_PRINT("  dna_patterns[1]: %1.*s\n", (int) aa_pattern_size*3, dna_patterns[1]);
    return 2;
}