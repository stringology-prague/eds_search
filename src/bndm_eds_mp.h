#include <stddef.h>

int bndm_eds_mp_search(const unsigned char *teds,
                       const size_t len,
                       unsigned char *pattern0,
                       unsigned char *pattern1,
                       unsigned int m);

int bndm_eds_mp_run(const unsigned char *teds,
                    const size_t len,
                    const unsigned char *pattern,
                    const size_t m,
                    const int loops);