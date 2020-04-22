#include <string>
#include <vector>
#include <iostream>
#include <map>
#include "sam_t.h"
#include "GenomicRange.h"

extern bool debug;

class clusterSC {

	// mètodes
	public:
		// Constructor,  <vector_of_reads>,  chromosome, position
		clusterSC( const std::string&, const std::string&, int, float, int, int);
		void extractSC( std::map<std::string, std::vector<sam_t>>&, std::map<std::string, std::vector<GenomicRange>>&, std::map<std::string, std::vector<GenomicRange>>&, std::string& );
		int getReadLength();
	private:
		std::string bamFile;
		std::string outDir;
		int minSoftLength;
		float softClipLenPerc;
		int minSoftClusterSize;
		int totalSR;
		int readLength;
};
