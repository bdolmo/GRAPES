
#define BASE_MASK 0x3 /* binary: 11 */

#include <iostream>
#include <string>
#include <vector>
#include <regex>

class dna2bit {

	public:
		dna2bit( std::string& );

		uint8_t*compressDNA( std::string&);
		std::string decompressDNA( uint8_t* , int);		
		
	private:
	    std::string dna_str;
    	    uint8_t* m_data;
};
