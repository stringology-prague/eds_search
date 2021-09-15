#include <stddef.h>

int bndm_eds_mp_search(unsigned char *pattern0, unsigned char *pattern1, unsigned int m);

int bndm_eds_mp_run(const unsigned char *pattern,
                    const size_t m,
                    const int loops);