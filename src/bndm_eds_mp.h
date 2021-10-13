#include <stddef.h>

int bndm_eds_mp_search(const unsigned char *text,
                       const size_t len,
                       unsigned char patterns[MAX_DNA_PATTERNS][MAX_PATTERN_LENGTH],
                       unsigned int num_patterns,
                       unsigned int m);

int bndm_eds_mp_run(const unsigned char *teds,
                    const size_t len,
                    const unsigned char *aa_pattern,
                    const size_t m,
                    const int loops);