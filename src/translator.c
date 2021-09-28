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
	unsigned int segmentStart = 0;
	unsigned int rbPointer = 0;
	unsigned char elementNum = 0;
	unsigned int length = 0;
	unsigned int wbPointer = 1;

    unsigned char *writeBuffer = malloc((eds_size + 1000000)*sizeof(unsigned char *));

	while (rbPointer < eds_size)
	{
		c = readBuffer[rbPointer++];
		//if(rbPointer>2814597)
			//printf("c = %c, rbPointer = %d, wbPointer = %d\n",c,rbPointer,wbPointer);
		if (length && (c == '{' || c == '}' || c == ','))
		{
            wbPointer += byteEncodeInt(writeBuffer + wbPointer, length-1);
//            strncpy(tmpBuffer, readBuffer + rbPointer - length, length - 1);
            strncpy(writeBuffer + wbPointer, readBuffer + rbPointer - length, length - 1);
            wbPointer += length-1;
            elementNum++;
            length = 0;
		}

		if (c == '{' || c == '}')
		{
			writeBuffer[segmentStart] = (unsigned char)elementNum;
			segmentStart = wbPointer++;
			elementNum = 0;
			if (c == '}' && readBuffer[rbPointer] == '{')
				rbPointer++;
		}

		length++;
	}
    if (c == '}'){
        *translated_eds_size = wbPointer;
        return writeBuffer;
    }
    length--;
	writeBuffer[segmentStart] = (unsigned char)elementNum+1;
	
	//writeBuffer[wbPointer + 1] = (unsigned char)(length & 0x00ff);
	//writeBuffer[wbPointer] = (unsigned char)(length >> 8); wbPointer += 2;
    wbPointer += byteEncodeInt(writeBuffer + wbPointer, length);
	strncpy(writeBuffer + wbPointer, readBuffer + rbPointer - length, length);
	wbPointer += length;
    *translated_eds_size = wbPointer;
    return writeBuffer;
}
