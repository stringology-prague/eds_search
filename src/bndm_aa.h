#include <stddef.h>

/**
 * BNDM-EDS matching function that accepts a string of amino acids as a <pattern> and performs exact search on an
 * elastic degenerate string <text> of DNA symbols.
 *
 * The amino acid string <pattern> is internally translated to compressed IUPAC, which results in two sub-patterns,
 * since each amino acid can be represented by at most two IUPAC codons. Special characters can be used on the amino
 * acid pattern to encode START and STOP sequences, see definitions of AA_START and AA_STOP. These two patterns are
 * then preprocessed into masking vectors into DNA alphabet and used for matching.
 *
 * @example Amino acid pattern "CRLV" will internally translate into two compressed IUPAC patterns "TGY CGN CTN GTN"
 * and "TGY AGR TTR GTN" (without spaces).
 *
 * See the table AA_TO_COMPR_IUPAC_SYMBOLS for transation from amino acid pattern into two compressed
 * IUPAC patterns, and table IUPAC_SYMBOLS_TO_BASES for transation from compressed IUPAC to DNA bases.
 *
 * @param text elastic degenerate string of DNA sequences, expected characters: 'A', 'C', 'G', 'T', '{', '}' and ','
 * @param len length of text
 * @param pattern amino acid pattern, expected characters are from "ABCDEFGHIKLMNPQRSTVWYZ", AA_START and AA_STOP
 * @param m length of pattern, maximum length is 10 for 32bit architectures, or 21 for 64bit architectures
 * @return number of matches
 */
int bndm_eds_aa_search(const unsigned char *text,
                       const size_t len,
                       const unsigned char *pattern,
                       const size_t m);

/**
 * BNDM-EDS matching function that accepts a string of DNA codons in compressed IUPAC format <pattern0,pattern1> and
 * performs exact search on an elastic degenerate string <text> of DNA symbols.
 *
 * The two sub-patterns are preprocessed into masking vectors into DNA alphabet and used for matching. Matches and
 * partial matches at every 3 positions of the patterns are merged between the patterns.
 *
 * @param text elastic degenerate string of DNA sequences, expected characters: 'A', 'C', 'G', 'T', '{', '}' and ','
 * @param len length of text
 * @param pattern compressed IUPAC pattern, expected characters are from "ABCDGHKMNRSTUVWY"
 * @param m length of pattern, maximum length is 10 for 32bit architectures, or 21 for 64bit architectures
 * @return number of matches
 */
int bndm_eds_iupac_search(const unsigned char *text,
                          const size_t len,
                          const unsigned char *pattern0,
                          const unsigned char *pattern1,
                          const size_t m);

/**
 * Wrapper function that runs either bndm_eds_aa_search(...) or bndm_eds_iupac_search(...). If the <pattern1> is NULL,
 * then the argument <pattern0> is treated as amino acid string and bndm_eds_aa_search is run.
 *
 * Otherwise the patterns are treated as IUPAC and bndm_eds_iupac_search is run. In this case, both patterns must be
 * of the same length <m>.
 *
 * This wrapper performs time measurements and runs the the underlying search function <loops>-times.
 *
 * @param teds Translated EDS text to search on
 * @param len Length of <teds>
 * @param pattern0 AA (see bndm_eds_aa_search) or IUPAC pattern (see bndm_eds_iupac_search)
 * @param pattern1 NULL or IUPAC pattern (see bndm_eds_iupac_search)
 * @param m length of pattern (for maximum allowed, see bndm_eds_aa_search or bndm_eds_iupac_search)
 * @param loops number of repeated execution of the underlying search algorithm
 * @return number of matches (per single run)
 */
int bndm_eds_aa_run(const unsigned char *teds,
                    const size_t len,
                    const unsigned char *pattern0,
                    const unsigned char *pattern1,
                    const size_t m,
                    const int loops);
