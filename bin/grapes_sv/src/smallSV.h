#include <iostream>
#include <string> // carreguem la llibreria d'strings
#include <fstream> // llibreria per obrir inputs
#include <vector>
#include <list>
#include <algorithm>
#include <omp.h>

#include "GenomicRange.h"
#include "pe_adjacency.h"
#include "ReadPairInfo.h"
#include "sam_t.h"
#include "vcf_t.h"

 class smallSV {

	// Mètodes
	public:
		// Constructor
		smallSV(std::string&, std::string&, int, int, int);

		void callSmallSV( std::map<string, std::vector<sam_t>>&, std::map<std::string, std::vector<GenomicRange>>&, std::string& fastq, std::vector<vcf_t>&, std::string&, std::string&);
		void callSmallSvExome (std::map<string, vector<sam_t>>&, std::map<std::string, std::vector<GenomicRange>>&, std::string&, std::vector<vcf_t>&, std::string&, std::string&);


	private:
		std::string bamFile;
		std::string reference;
		int minOlapBases;
		int maxMismatches;
		int readLength;
 };
