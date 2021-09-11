#include "globals.h"

#define BASES 4

extern unsigned char IUPAC_SYMBOLS_TO_BASES[SIGMA][BASES];

void init_IUPAC_SYMBOLS_TO_BASES();

#define COMPR_IUPAC_MAX_VARIATIONS 2
#define COMPR_IUPAC_LEN 3
#define AA_START 0x00
#define AA_STOP 0xFF

extern unsigned char AA_TO_COMPR_IUPAC_SYMBOLS[SIGMA][COMPR_IUPAC_MAX_VARIATIONS][COMPR_IUPAC_LEN];

void init_AA_TO_COMPR_IUPAC_SYMBOLS();

size_t translate_aa_pattern(const unsigned char *aa_pattern,
                            size_t aa_pattern_size,
                            unsigned char dna_patterns[][MAX_PATTERN_LENGTH],
                            size_t max_dna_patterns,
                            size_t max_dna_pattern_length);
