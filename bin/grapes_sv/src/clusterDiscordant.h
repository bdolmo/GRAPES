#include <string>
#include <vector>
#include <iostream>
#include <map>
#include "sam_t.h"
#include "ReadPairInfo.h"
#include "GenomicRange.h"


class clusterDiscordant {

	// mètodes
	public:
		// Constructor
		clusterDiscordant( const std::string&, const std::string&, const std::string&, int, int, long&, const std::string&, std::string, int&);

		std::tuple<double, double> calculateInsertSizeLimits();
		void extractDiscordant(std::map<std::string, std::vector<GenomicRange>>&);
		void cluster( std::string&, std::map<std::string, std::vector<GenomicRange>>& );

		// getters
		float getMeanCoverage();
		std::vector<int> getSomaticCounts();
		std::vector<int> getGerminalCounts();
		int getTotalFR();
		int getTotalRF();
		int getTotalFF();
		int getTotalRR();
		int getBinSize();
		int getTotalSR();
		int getReadLength();

	private:
		std::string bamFile;
       	std::string outdir;
		std::string sampleName;
		int minClusterSize;
		int numSDs;
		long genomeSize;
		std::string genome;
		std::string getCountsByWindow;
		int maxInsertSize;
		double mean;
		double stdev;
		double upper_limit;
		double mean_coverage;
		int totalDiscordants;
		std::map< std::string, ReadPairInfo>  readPair_map;

		int totalFR = 0;
		int totalRF = 0;
		int totalFF = 0;
		int totalRR = 0;
		int totalSR = 0;

		int binSize;
		std::vector<int> readCounts_somatic;
		std::vector<int> readCounts_germinal;
		int readLength;
};
