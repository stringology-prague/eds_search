//Definitions
#define SIGMA 256
#define MAX_PATTERN_LENGTH (sizeof(int) * 8)
#define MAX_AA_PATTERN_LENGTH ((sizeof(int) * 8) / 3)
#define MAX_DNA_PATTERNS (2)

//Byte encoder limits for single codeword lengths
#define ENC_1B 200
#define ENC_2B 7880
#define ENC_3B 1318600

//Byte decoder limits for single codeword lengths
#define DEC_1B 200
#define DEC_2B 230
#define DEC_3B 250

//Globals
extern unsigned char *readBuffer;
extern unsigned int fSize;
extern unsigned int rbPointer;

extern unsigned char *writeBuffer;
extern unsigned int wbPointer;
extern unsigned int aPointer;

#ifdef DEBUG
#define DEBUG_PRINT(...) do{ fprintf( stderr, __VA_ARGS__ ); } while( 0 )
#else
#define DEBUG_PRINT(...) do{ } while ( 0 )
#endif
