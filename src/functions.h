#include <stdio.h>

void byteEncodeInt(unsigned int x);
unsigned int byteDecodeInt();
int getFileSize(FILE *f);
void readInputFile(char *fName, unsigned char **readBuffer, int *size);
void writeOutputFile(char *fName, unsigned char **outBuffer, int outSize);
void printString(unsigned int start, unsigned char length);
unsigned int randomSelectPattern(unsigned char* pattern, unsigned short pattLen, unsigned char* file, unsigned int fileSize);
int min(int a, int b);