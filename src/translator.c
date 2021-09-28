#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <string.h>
#include <sys/resource.h>
#include "globals.h"
#include "functions.h"

unsigned char* translate(unsigned char *readBuffer, int eds_size,  int *translated_eds_size)
{
	unsigned char c;
	unsigned int rbPointer = 0;
	unsigned int elementNum = 0;
	unsigned int length = 0;
	unsigned int wbPointer = 0;

    unsigned char *writeBuffer = malloc((eds_size + 1000000)*sizeof(unsigned char *));
    unsigned int tmpPointer = 0;
    unsigned char *tmpBuffer = malloc((eds_size + 1000000)*sizeof(unsigned char *));

	while (rbPointer < eds_size)
	{
		c = readBuffer[rbPointer++];
		//if(rbPointer>2814597)
			//printf("c = %c, rbPointer = %d, wbPointer = %d\n",c,rbPointer,wbPointer);
		if (length && (c == '{' || c == '}' || c == ','))
		{
            tmpPointer += byteEncodeInt(tmpBuffer + tmpPointer, length-1);
            memcpy(tmpBuffer + tmpPointer, readBuffer + rbPointer - length, length - 1);
            tmpPointer += length-1;
            elementNum++;
            length = 0;
        }

        if (c == '{' || c == '}')
		{
            wbPointer += byteEncodeInt(writeBuffer + wbPointer, elementNum);
            memcpy(writeBuffer + wbPointer, tmpBuffer, tmpPointer);
            wbPointer += tmpPointer;
            tmpPointer = 0;
			elementNum = 0;
			if (c == '}' && readBuffer[rbPointer] == '{')
				rbPointer++;
		}

		length++;
	}
    if (c == '}'){
        *translated_eds_size = wbPointer;
        free(tmpBuffer);
        return writeBuffer;
    }
    length--;
    wbPointer += byteEncodeInt(writeBuffer + wbPointer, elementNum+1) ;
    wbPointer += byteEncodeInt(writeBuffer + wbPointer, length);
	memcpy(writeBuffer + wbPointer, readBuffer + rbPointer - length, length);
	wbPointer += length;
    *translated_eds_size = wbPointer;
    free(tmpBuffer);
    return writeBuffer;
}
