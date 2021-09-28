#include <stdio.h>

unsigned int byteEncodeInt(unsigned char *teds, const unsigned int x_in);
unsigned int byteDecodeInt(const unsigned char *teds, unsigned int *x_out);
int getFileSize(FILE *f);
unsigned char * readInputFile(char *fName, int *size);
void writeOutputFile(char *fName, unsigned char **outBuffer, int outSize);
unsigned int randomSelectPattern(unsigned char* pattern, unsigned short pattLen, const unsigned char* file, unsigned int fileSize);
int min(int a, int b);