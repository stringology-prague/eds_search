#include <string.h>

#include "protein_table.h"

unsigned char IUPAC_SYMBOLS_TO_BASES[SIGMA][BASES+1];

void init_IUPAC_SYMBOLS_TO_BASES(){
    memset(IUPAC_SYMBOLS_TO_BASES,0,SIGMA * (BASES+1) * sizeof (unsigned char));

    IUPAC_SYMBOLS_TO_BASES['A'][0] = 'A';
    IUPAC_SYMBOLS_TO_BASES['C'][0] = 'C';
    IUPAC_SYMBOLS_TO_BASES['G'][0] = 'G';
    IUPAC_SYMBOLS_TO_BASES['T'][0] = 'T';
    IUPAC_SYMBOLS_TO_BASES['U'][0] = 'U';

    IUPAC_SYMBOLS_TO_BASES['W'][0] = 'A';
    IUPAC_SYMBOLS_TO_BASES['S'][0] = 'C';
    IUPAC_SYMBOLS_TO_BASES['M'][0] = 'A';
    IUPAC_SYMBOLS_TO_BASES['K'][0] = 'G';
    IUPAC_SYMBOLS_TO_BASES['R'][0] = 'A';
    IUPAC_SYMBOLS_TO_BASES['Y'][0] = 'C';

    IUPAC_SYMBOLS_TO_BASES['W'][1] = 'T';
    IUPAC_SYMBOLS_TO_BASES['S'][1] = 'G';
    IUPAC_SYMBOLS_TO_BASES['M'][1] = 'C';
    IUPAC_SYMBOLS_TO_BASES['K'][1] = 'T';
    IUPAC_SYMBOLS_TO_BASES['R'][1] = 'G';
    IUPAC_SYMBOLS_TO_BASES['Y'][1] = 'T';

    IUPAC_SYMBOLS_TO_BASES['B'][0] = 'C';
    IUPAC_SYMBOLS_TO_BASES['D'][0] = 'A';
    IUPAC_SYMBOLS_TO_BASES['H'][0] = 'A';
    IUPAC_SYMBOLS_TO_BASES['V'][0] = 'A';

    IUPAC_SYMBOLS_TO_BASES['B'][1] = 'G';
    IUPAC_SYMBOLS_TO_BASES['D'][1] = 'G';
    IUPAC_SYMBOLS_TO_BASES['H'][1] = 'C';
    IUPAC_SYMBOLS_TO_BASES['V'][1] = 'C';

    IUPAC_SYMBOLS_TO_BASES['B'][2] = 'T';
    IUPAC_SYMBOLS_TO_BASES['D'][2] = 'T';
    IUPAC_SYMBOLS_TO_BASES['H'][2] = 'T';
    IUPAC_SYMBOLS_TO_BASES['V'][2] = 'G';

    IUPAC_SYMBOLS_TO_BASES['N'][0] = 'A';
    IUPAC_SYMBOLS_TO_BASES['N'][1] = 'C';
    IUPAC_SYMBOLS_TO_BASES['N'][2] = 'G';
    IUPAC_SYMBOLS_TO_BASES['N'][3] = 'T';

    return;
}