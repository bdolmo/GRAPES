#include <string>
#include <iostream>
#include <vector>
#include <map>

#include "sample_metrics_t.h"
#include "sample_metrics_region_t.h"
#include "bed_t.h"

class ExtractDepth {

	// Methods
	public:
		// Constructor
		ExtractDepth( std::vector<bed_t>&, std::string&, std::string&, std::string&,  bool, bool, bool, std::map<std::string, std::vector<sample_metrics_region_t>>&, std::map<std::string,	std::vector<sample_metrics_region_t>>&, std::string  );

		void extract(std::map<std::string, std::vector<sample_metrics_region_t>>&, std::map<std::string, std::vector<sample_metrics_region_t>>& );

		// Getters
		int getTotalReads();
		int getTotalReadsX();
		double getMeanCount();
		double getMeanCountX();
		double getMeanCov();
		double getMeanCovX();
		double getMeanInsertSize();
		double getStdDevInsertSize();
		
	private:
		std::vector<bed_t> RegionsVector;
		std::string sampleNamePath;
		std::string sampleName;
		std::string genome;
		bool is_first;
		bool extractCounts;
		bool extractCoverage;
		std::map<std::string, std::vector<sample_metrics_region_t>> mapSummaryCount;
		std::map<std::string, std::vector<sample_metrics_region_t>> mapSummaryCoverage;
		std::string outdir;

		int totalReads;
		int totalReadsX;
		double meanCount;
		double meanCountX;
		double meanCov;
		double meanCovX;
		double meanIsize;
		double stdDevIsize;
};
