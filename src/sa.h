#include <stddef.h>

int SA_search(const unsigned char *teds,
              const size_t len,
              unsigned char *x,
              unsigned int m);

int SA_run(const unsigned char *teds,
           const size_t len,
           const unsigned char *pattern,
           const size_t m, const int loops);
