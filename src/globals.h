//Definitions
#define SIGMA 256

//Byte encoder limits for single codeword lengths
#define ENC_1B 200
#define ENC_2B 7880
#define ENC_3B 1318600

//Byte decoder limits for single codeword lengths
#define DEC_1B 200
#define DEC_2B 230
#define DEC_3B 250

//Globals
extern unsigned char* readBuffer;
extern unsigned int fSize;
extern unsigned int rbPointer;

extern unsigned char* writeBuffer;
extern unsigned int wbPointer;
extern unsigned int aPointer;