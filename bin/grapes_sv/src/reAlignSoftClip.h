#include <iostream>
#include <vector>
#include <string> 
#include <fstream>
#include "vcf_t.h"

class reAlignSoftClip {

	public:
		// Constructor
		reAlignSoftClip(std::string, std::string ,int, int, std::string, std::string, std::string, std::string, std::string, int, int, int);

		void Align(std::vector<vcf_t>&);

	private:
		std::string read;
		std::string chr;
		int pos;
		int end;
		std::string cigar;
		std::string ref_sequence;
		std::string softclip_type;	
		std::string VCF;
		std::string svlength;
		int posA;
		int posB;
		int mapq;

};

