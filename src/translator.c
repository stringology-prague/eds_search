#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <string.h>
#include <sys/resource.h>
#include "globals.h"
#include "functions.h"

void translate()
{
	unsigned char c, prevC;
	unsigned int segmentStart = 0;
	unsigned char elementNum = 0;
	unsigned int length = 0;
	wbPointer = 1;

	while (rbPointer < fSize)
	{
		c = readBuffer[rbPointer++];
		//if(rbPointer>2814597)
			//printf("c = %c, rbPointer = %d, wbPointer = %d\n",c,rbPointer,wbPointer);
		if (length && (c == '{' || c == '}' || c == ','))
		{
            byteEncodeInt(length-1);
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
		prevC = c;
	}
    if (c == '}'){
        return;
    }
    length--;
	writeBuffer[segmentStart] = (unsigned char)elementNum+1;
	
	//writeBuffer[wbPointer + 1] = (unsigned char)(length & 0x00ff);
	//writeBuffer[wbPointer] = (unsigned char)(length >> 8); wbPointer += 2;
	byteEncodeInt(length);
	strncpy(writeBuffer + wbPointer, readBuffer + rbPointer - length, length);
	wbPointer += length;
}
