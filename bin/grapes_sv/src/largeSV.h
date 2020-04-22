#include <iostream>
#include <string> // carreguem la llibreria d'strings
#include <fstream> // llibreria per obrir inputs
#include <vector>
#include <list>
#include <algorithm>

#include "GenomicRange.h"
#include "pe_adjacency.h"
#include "ReadPairInfo.h"
#include "sam_t.h"
#include "vcf_t.h"

extern bool debug;


 class largeSV {

	// Mètodes
	public:
		// Constructor
		largeSV(std::string&, std::string&, int, int);

		void callStructVar( std::map<std::string, std::vector<GenomicRange>>&, std::map<std::string, std::vector<GenomicRange>>&, std::map<string, std::vector<sam_t>>& breakReadClusters, std::vector<vcf_t>&, std::string svtype, std::string);

	private:
		std::string bamFile;
		std::string reference;
		int minOlapBases;
		int maxMismatches;
 };



