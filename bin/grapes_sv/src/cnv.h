#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include "vcf_t.h"

#include "GenomicRange.h"

class Cnv {

	public:
		// Declaring constructor
		Cnv( std::string&, std::string&, int, int, int );

		// Methods	
		void normalize_gc();
		void normalize_mappability();
		std::vector<vcf_t> dumpCNV (std::string, std::string&, std::string&);

		std::vector<GenomicRange> getRatios();
		void segment();
	
	private:
		std::string outDir;
		std::string raw_counts;
		int medianSomaticCounts;
		int medianGerminalCounts;
		int n_states;
		
};
