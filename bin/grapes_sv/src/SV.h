#include <iostream>
#include <string>
#include <vector>
#include "sam_t.h"
#include "vcf_t.h"
#include <fstream>

class SV {

	public:

		// Constructor,  <vector_of_reads>,  chromosome, position
		SV( std::vector<sam_t>&, std::string, std::string, int, int, std::string, int, int, std::string, std::string&, std::string, int, double, double);

		bool classify(std::vector<vcf_t>&, std::string&);

	private:
		std::vector<sam_t> reads;
		std::string chrA;
		std::string chrB;
		int posA;
		int posB;
		std::string softclip_type;
		int total_reads;
		int total_assembled;
		std::string bamFile;
		std::string ref_sequence;
		std::string svlength;
		int nDiscordants;
		double pvalue_discordant;
		double kmer_diversity;


};
