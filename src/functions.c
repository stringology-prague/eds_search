#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <string.h>
#include <sys/resource.h>
#include "globals.h"

unsigned int byteEncodeInt(unsigned char *writeBuffer, unsigned int x)
{
	unsigned int pos = 0;
	//1B codeword
	if (x < ENC_1B)
	{
		writeBuffer[pos++] = (unsigned char)x;
		return pos;
	}

	//2B codeword
	if (x < ENC_2B)
	{
		x -= ENC_1B;
		writeBuffer[pos++] = (unsigned char)((x>>8)+DEC_1B);
		writeBuffer[pos++] = (unsigned char)(x & 0x000000ff);
		return pos;
	}

	//3B codeword
	if (x < ENC_3B)
	{
		x -= ENC_2B;
		writeBuffer[pos++] = (unsigned char)((x >> 16) + DEC_2B);
		writeBuffer[pos++] = (unsigned char)((x >> 8) & 0x000000ff);
		writeBuffer[pos++] = (unsigned char)(x & 0x000000ff);
		return pos;
	}

	//4B codeword
	x -= ENC_3B;
	writeBuffer[pos++] = (unsigned char)((x >> 24) + DEC_3B);
	writeBuffer[pos++] = (unsigned char)((x >> 16) & 0x000000ff);
	writeBuffer[pos++] = (unsigned char)((x >> 8) & 0x000000ff);
	writeBuffer[pos++] = (unsigned char)(x & 0x000000ff);
	return pos;
}

unsigned int byteDecodeInt(const unsigned char *writeBuffer, unsigned int *x_out)
{
	unsigned int pos = 0;
	unsigned char c = writeBuffer[pos++], x;

	if (c < DEC_1B)
	{
		*x_out = c;
		return pos;
	}

	if (c < DEC_2B)
	{
		x = (unsigned int)(c - DEC_1B);
		x <<= 8;
		x |= (unsigned int)writeBuffer[pos++];
		x += ENC_1B;
		*x_out = x;
		return pos;
	}

	if (c < DEC_3B)
	{
		x = (unsigned int)(c - DEC_2B);
		x <<= 8;
		x |= (unsigned int)writeBuffer[pos++];
		x <<= 8;
		x |= (unsigned int)writeBuffer[pos++];
		x += ENC_2B;
		*x_out = x;
		return pos;
	}

	x = (unsigned int)(c - DEC_3B);
	x <<= 8;
	x |= (unsigned int)writeBuffer[pos++];
	x <<= 8;
	x |= (unsigned int)writeBuffer[pos++];
	x <<= 8;
	x |= (unsigned int)writeBuffer[pos++];
	x += ENC_3B;
	*x_out = x;
	return pos;
}

int getFileSize(FILE *f)
{
	int size;
	fseek(f, 0, SEEK_END);
	size = ftell(f);
	fseek(f, 0, SEEK_SET);
	return size;
}

unsigned char* readInputFile(char *fName, int *size)
{
    int unused __attribute__((unused));

    FILE *fr = fopen(fName, "r");
    *size = getFileSize(fr);
	unsigned char *readBuffer = (unsigned char*)malloc((*size) * sizeof(char));
    unused = fread(readBuffer, 1, *size, fr);
	fclose(fr);
	return readBuffer;
}

void writeOutputFile(char *fName, unsigned char **outBuffer, int outSize)
{
	FILE *fw = fopen(fName, "w");
	fwrite(*outBuffer, sizeof(char), outSize, fw);
	fclose(fw);
}

unsigned int randomSelectPattern(unsigned char* pattern, unsigned int pattLen, unsigned char* file, unsigned int fileSize) {

	unsigned int  start = rand() % fileSize;

	unsigned int pPos = 0;
	unsigned fPos = start;
	
	while (pPos < pattLen && fPos < fileSize) {
		if (file[fPos] == 'N' || file[fPos] == ',' || file[fPos] == '{' || file[fPos] == '}') {
			fPos++;
			continue;
		}
		pattern[pPos++] = file[fPos++];
	}
	pattern[pPos] = '\0';
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