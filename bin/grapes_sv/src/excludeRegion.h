#include <string>
#include <vector>
#include <iostream>
#include <map>
#include "pe_adjacency.h"
#include "Utils.cpp"
#include "ReadPairInfo.h"
#include "GenomicRange.h"
class excludeRegion {

	// mètodes
	public:
		// Constructor
		excludeRegion( const std::string&);
		std::map<std::string, std::vector<GenomicRange>> readAndSave();
	
	private:
		std::string regionsExclude;
		std::map<std::string, GenomicRange> mapOfRegions;
};


