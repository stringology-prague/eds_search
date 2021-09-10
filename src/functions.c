#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <string.h>
#include <sys/resource.h>
#include "globals.h"

void byteEncodeInt(unsigned int x)
{
	//1B codeword
	if (x < ENC_1B)
	{
		writeBuffer[wbPointer++] = (unsigned char)x;
		return;
	}

	//2B codeword
	if (x < ENC_2B)
	{
		x -= ENC_1B;
		writeBuffer[wbPointer++] = (unsigned char)((x>>8)+DEC_1B);
		writeBuffer[wbPointer++] = (unsigned char)(x & 0x000000ff);
		return;
	}

	//3B codeword
	if (x < ENC_3B)
	{
		x -= ENC_2B;
		writeBuffer[wbPointer++] = (unsigned char)((x >> 16) + DEC_2B);
		writeBuffer[wbPointer++] = (unsigned char)((x >> 8) & 0x000000ff);
		writeBuffer[wbPointer++] = (unsigned char)(x & 0x000000ff);
		return;
	}

	//4B codeword
	x -= ENC_3B;
	writeBuffer[wbPointer++] = (unsigned char)((x >> 24) + DEC_3B);
	writeBuffer[wbPointer++] = (unsigned char)((x >> 16) & 0x000000ff);
	writeBuffer[wbPointer++] = (unsigned char)((x >> 8) & 0x000000ff);
	writeBuffer[wbPointer++] = (unsigned char)(x & 0x000000ff);
}

unsigned int byteDecodeInt()
{
	unsigned int x;
	unsigned char c = writeBuffer[aPointer++];

	if (c < DEC_1B)
		return (unsigned int)c;

	if (c < DEC_2B)
	{
		x = (unsigned int)(c - DEC_1B);
		x <<= 8;
		x |= (unsigned int)writeBuffer[aPointer++];
		x += ENC_1B;
		return x;
	}

	if (c < DEC_3B)
	{
		x = (unsigned int)(c - DEC_2B);
		x <<= 8;
		x |= (unsigned int)writeBuffer[aPointer++];
		x <<= 8;
		x |= (unsigned int)writeBuffer[aPointer++];
		x += ENC_2B;
		return x;
	}

	x = (unsigned int)(c - DEC_3B);
	x <<= 8;
	x |= (unsigned int)writeBuffer[aPointer++];
	x <<= 8;
	x |= (unsigned int)writeBuffer[aPointer++];
	x <<= 8;
	x |= (unsigned int)writeBuffer[aPointer++];
	x += ENC_3B;
	return x;
}

int getFileSize(FILE *f)
{
	int size;
	fseek(f, 0, SEEK_END);
	size = ftell(f);
	fseek(f, 0, SEEK_SET);
	return size;
}

void readInputFile(char *fName, unsigned char **readBuffer, int *size)
{
    int unused __attribute__((unused));

    FILE *fr = fopen(fName, "r");
    *size = getFileSize(fr);
    *readBuffer = (unsigned char*)malloc((*size) * sizeof(char));
    unused = fread(*readBuffer, 1, *size, fr);
	fclose(fr);

	//printf("readInputFile: %c\n",*readBuffer[0]);
}

void writeOutputFile(char *fName, unsigned char **outBuffer, int outSize)
{
	FILE *fw = fopen(fName, "w");
	fwrite(*outBuffer, sizeof(char), outSize, fw);
	fclose(fw);
}

void printString(unsigned int start, unsigned char length)
{
	printf("start = %d, length = %d", start, length);
	
	for (unsigned int i = start; i < (start + length); i++)
		printf("i = %d: %c\n", i, writeBuffer[i]);

	printf("\n");
}

unsigned int randomSelectPattern(unsigned char** pattern, unsigned int pattLen, unsigned char* file, unsigned int fileSize) {

	unsigned int  start = rand() % fileSize;
	
	(*pattern) = (char*)malloc((pattLen) * sizeof(char));

	unsigned int pPos = 0;
	unsigned fPos = start;
	
	while (pPos < pattLen && fPos < fileSize) {
		if (file[fPos] == 'N' || file[fPos] == ',' || file[fPos] == '{' || file[fPos] == '}') {
			fPos++;
			continue;
		}
		(*pattern)[pPos++] = file[fPos++];
	}
	(*pattern)[pPos] = '\0';
	//memcpy(*pattern, file + start, pattLen);

	return start;
}

int min(int a, int b)
{
    if (a < b)
    {
        return a;
    }
    return b;
}