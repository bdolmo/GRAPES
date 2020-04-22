
#define BASE_MASK 0x3 /* binary: 11 */

#include <iostream>
#include <string>
#include <vector>
#include <regex>
#include "dna2bit.h"
#include <cstring>

using namespace std;

enum
{
    BASE_A = 0x0, /* binary: 00 */
    BASE_C = 0x1, /* binary: 01 */
    BASE_G = 0x2, /* binary: 10 */
    BASE_T = 0x3, /* binary: 11 */
};


uint8_t* dna2bit::compressDNA( std::string& dna_str) {

	size_t dna_bytes = (dna_str.length() / 4) + (dna_str.length() % 4 != 0);
 
        uint8_t* m_data = new uint8_t[dna_bytes];
 
        memset(m_data, 0, dna_bytes);

	for(int i = 0; i < dna_str.length(); i++) {
		uint8_t shift = 6 - 2 * (i % 4);
		switch (dna_str[i]) {
          		case 'A':
		                m_data[i / 4] |= BASE_A << shift;
				break;
			case 'C':
		                m_data[i / 4] |= BASE_C << shift;
				break;
			case 'T':
		                m_data[i / 4] |= BASE_T << shift;
				break;
			case 'G':
		                m_data[i / 4] |= BASE_G << shift;
				break;		
		}
                shift = (shift == 0) ? 6 : shift - 2;
	}
	return m_data;
}


string dna2bit::decompressDNA( uint8_t* m_data, int seqLen) {

	string new_dna;
	
	for (size_t i = 0; i < seqLen; ++i) {
        
            uint8_t shift = 6 - 2 * (i % 4);
            uint8_t mask = BASE_MASK << shift;
 
            /* get the i-th DNA base */
            uint8_t base = (m_data[i / 4] & mask) >> shift;

            switch (base) {
                case BASE_A:
		    new_dna = new_dna + "A";
                    break;
                case BASE_C:
		    new_dna = new_dna + "C";
                    break;
                case BASE_G:
		    new_dna = new_dna + "G";
                    break;
                case BASE_T:
		    new_dna = new_dna + "T";
                    break;
           }
        }
 
        new_dna[seqLen] = '\0';
	return new_dna;
}

