#include <iostream>
#include <string>
#include <vector>
#include <map>
#include "excludeRegion.h"
#include "Utils.cpp"
#include "GenomicRange.h"
#include <regex>

std::map<std::string, std::vector<GenomicRange>> excludeRegion::readAndSave() {

	ifstream inputRegionsFile;
	inputRegionsFile.open( regionsExclude );
	std::string line;

	std::map<std::string, std::vector<GenomicRange>> blackRegions;

	if (inputRegionsFile.is_open()) {
    		while ( getline (inputRegionsFile, line) ) {

			std::vector<std::string> tmp = split( line, '\t' );
	
			GenomicRange grange;
			grange.chromosomeA      = tmp[0];
			grange.chromosome_start = atoi(tmp[1].c_str());
			grange.chromosome_end   = atoi(tmp[2].c_str());

			grange.chr_centromere_start = grange.chromosome_start;
			grange.chr_centromere_end   = grange.chromosome_end;

        		std::smatch m;
			std::string tmps = tmp[3];
			std::regex centromere  ("centromere");

			if (grange.chromosome_start-1000000 > 0 && std::regex_search (tmp[3], m, centromere)) {
				grange.chromosome_start     = atoi(tmp[1].c_str())-1000000;
				grange.chromosome_end       = atoi(tmp[2].c_str())+1000000;
				grange.chr_centromere_start = atoi(tmp[1].c_str())-1000000;
				grange.chr_centromere_end   = atoi(tmp[2].c_str())+1000000;
			}  
			blackRegions[grange.chromosomeA].push_back(grange);
		}
	}
	else {
		std::cout << "ERROR: Unable to open " << regionsExclude << "\n";
	}
	
	if (!blackRegions.empty()) {
		std::cout << " INFO: " << regionsExclude << " successfully loaded!" << endl;
	}
	else {
		std::cout << " WARNING: " << regionsExclude << " could not be loaded" << endl;
	}
	
	return blackRegions;
}
