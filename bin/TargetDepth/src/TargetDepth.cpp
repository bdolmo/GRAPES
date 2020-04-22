#include <iostream>
#include <string.h>
#include <fstream>
#include <vector>
#include <map>
#include <sstream>
#include <algorithm>

// Boost library
#include <boost/program_options.hpp>
#include <boost/iterator/transform_iterator.hpp>
#include <boost/algorithm/string/replace.hpp>
#include <boost/algorithm/string.hpp>

#include "SeqLib/RefGenome.h"
#include "SeqLib/BamReader.h"
#include "SeqLib/BamRecord.h"

// htslib
#include "htslib/htslib/faidx.h"
#include "htslib/htslib/kstring.h"
#include "htslib/htslib/khash.h"

#include "ExtractDepth.h"
#include "sample_metrics_t.h"
#include "sample_metrics_region_t.h"
#include "bed_t.h"

using namespace std;
using namespace SeqLib;
namespace po = boost::program_options;


string getFileName(const string& s);
vector<string> split( std::string const& original, char separator );
void bed2Vector( std::string& bed, vector<bed_t>& RegionsVector, RefGenome& ref);
float computeGC ( std::string chr, int start, int end, RefGenome& ref);
std::string IntToString ( int& );
//g++  -I htslib/ -lhts -I/home/bernat/C++/SeqLib/htslib/htslib -L/home/bernat/C++/SeqLib/htslib/ -lhts -I/home/bernat/C++/SeqLib -L/home/bernat/C++/SeqLib/bin -lseqlib -o class -std=c++11

// Declaring Constructor
ExtractDepth::ExtractDepth(std::vector<bed_t>& _RegionsVector, string& _sampleNamePath, string& _sampleName, string& _genome,  bool _is_first, bool _extractCounts, bool _extractCoverage, map<std::string, vector<sample_metrics_region_t>>&  _mapSummaryCount, map<std::string, vector<sample_metrics_region_t>>& _mapSummaryCoverage, std::string _outdir ) {
	RegionsVector = _RegionsVector;
	sampleNamePath = _sampleNamePath;
	sampleName = _sampleName;
	genome = _genome;
	is_first = _is_first;
 	extractCounts = _extractCounts;
	extractCoverage = _extractCoverage;
	mapSummaryCount= _mapSummaryCount;
	mapSummaryCoverage= _mapSummaryCoverage;
	outdir = _outdir;
}

int main (int argc, char* argv[]) {
	
	string bam_file;
	string sample_file;
	string reference;
	string bed;
	string outdir;
	
	bool extractCoverage;
	bool extractCounts;
	
	po::options_description description("\n Program: TargetDepth: A simple tool for extracting read depth from a list of BAM files \n Contact: bdelolmo@gencardio.com \n Usage");
	description.add_options()
	    ("input,i", po::value<std::string>(&bam_file), "BAM file as a single input")
	    ("list,l", po::value<std::string>(&sample_file), "File specifying multiple BAMs (with paths)")
	    ("genome,g", po::value<std::string>(&reference)->required(), "Genome reference FASTA")
	    ("bed,b", po::value<std::string>(&bed)->required(), "BED region file")
	    ("outdir,o", po::value<std::string>(&outdir)->default_value("."), "Output directory")
	    ("depth,d", "Extract per-base coverage")
	    ("counts,c", "Extract toal read count per region")
	    ("help,h", "Display help message\n");

	po::variables_map vm;
	po::positional_options_description input_options;
	po::store(po::parse_command_line(argc, argv, description), vm);
	po::store(po::command_line_parser(argc, argv).options(description).positional(input_options).run(), vm);
	time_t now = time(0);
	// convert now to string form
	char* dt = ctime(&now);
	try {

		if(vm.count("help")){
			std::cout << description;
			return 0;
		}
		if(vm.count("depth")){
			extractCoverage = true;
		}
		else {
			extractCoverage = false;
		}
		if(vm.count("counts")){
			extractCounts = true;
		}
		else {
			extractCounts = false;
		}
	        po::notify(vm);
		if(vm.count("input") && vm.count("list")) {
			cout << endl << " ERROR: choose between -i (for single BAM) or -l (for list of BAMs)" << endl;
			return 0;
		}
		else if(vm.count("input")){

		}
		else if(vm.count("list")){

		}
		else {
			cout << " ERROR: Missing input BAM/s. Please use -i option for single BAM or -l for multiple BAMs" << endl;
		}
		if(vm.count("bed")){

		}
		else {
			cout << " ERROR: Missing a BED regions file" << endl;
		}
		if(vm.count("genome")){

		}
		else {
			cout << " ERROR: Missing Reference Genome file in FASTA format" << endl;
		}

	} catch (po::error& e) {
		std::cout << description;
 		return 1;
	}
	RefGenome ref;
	ref.LoadIndex(reference);

	if (ref.LoadIndex(reference)) {
	}
	else {
		cout << " ERROR: please provide a valid Genome file in FASTA format" << endl;
		return 0;
	} 

	if (extractCoverage == false && extractCounts == false) {
		cout << "\nERROR: please provide -d (coverage per base) and/or -c (read count) depth extraction" << endl << description << endl;
		return 0;
	}
	std::cout << endl << " Executing TargetDepth on " << dt << endl;
	if(extractCoverage == true) {
		std::cout <<   " INFO: Extract coverage: true "  << endl;
	}
	if (extractCounts == true) {
		std::cout <<  " INFO: Extract counts: true "  << endl;
	}
	std::cout << " INFO: Reference Genome: " << reference << endl;
	std::cout << " INFO: BED region file: " << bed << endl;

	// Getting samples from file
	vector<string> samplePaths;

	if (vm.count("list")) {
		std::cout <<  " INFO: list of BAM Files: " << sample_file << endl;
		ifstream sampFile;
	  	sampFile.open (sample_file);
		if (sampFile.is_open()) {
			string line;
			while ( std::getline (sampFile, line)) {
				samplePaths.push_back(line);
			}
		}
	}
	if (vm.count("input")) {
		std::cout << " INFO: single BAM File: " << bam_file << endl;
		samplePaths.push_back(bam_file);
	}

	// Reading BED file and saving info for further accessing
	vector<bed_t> RegionsVector;
	bed2Vector(bed, RegionsVector, ref);
	map<std::string, std::vector<sample_metrics_region_t>> mapSummaryCount;

	map<std::string, std::vector<sample_metrics_region_t>> mapSummaryCoverage;

	string logfile = outdir + "/" + "summary_metrics.log";

	// Saving general metrics (e.g mean coverage, mean counts, etc) 
	vector<sample_metrics_t> GeneralMetrics;

	bool is_first = true;
	vector<string> sampleBasenames;
	string header;

	// Printing general metrics in a log file
	std::ofstream LOG;
	LOG.open (logfile, std::ofstream::out | std::ofstream::app);

	for (auto& f : samplePaths) {

		string sampleName = getFileName(f);

		ExtractDepth e1 (RegionsVector, f, sampleName, reference, is_first, extractCounts, extractCoverage, mapSummaryCount, mapSummaryCoverage, outdir);
		e1.extract(mapSummaryCount, mapSummaryCoverage);
		is_first = false;
		int total_reads  = e1.getTotalReads();
		int total_readsX = e1.getTotalReadsX();
		double mean_coverage = e1.getMeanCov();
		double mean_coverageX= e1.getMeanCovX();
		double mean_count = e1.getMeanCount();
		double mean_countX= e1.getMeanCountX();
		double mean_isize = e1.getMeanInsertSize();
		double sd_isize = e1.getStdDevInsertSize();
	        time_t now = time(0);
	        // convert now to string form
	        char* dt = ctime(&now);
		dt[strlen(dt) - 1] = '\0';

		cout << " INFO: [" << dt << "] " << sampleName << "\t" << total_reads << "\t" << total_readsX << "\t"  << mean_coverage << "\t" << mean_count << "\t" << mean_isize << "\t" << sd_isize << "\t" << mean_coverageX << "\t" << mean_countX << "\n";
		LOG << sampleName << "\t" << total_reads << "\t" << total_readsX  << "\t" << mean_coverage << "\t" << mean_count << "\t" << mean_isize << "\t" << sd_isize << "\t" << mean_coverageX << "\t" << mean_countX << "\n";

	}
  return 0;
}


void bed2Vector( std::string& bed, vector<bed_t>& RegionsVector, RefGenome& ref) {
	ifstream bedFile;
  	bedFile.open (bed);
        if (bedFile.is_open()) {
		string line;
		while ( std::getline (bedFile, line)) {
			vector<string> tmp = split(line, '\t');
			string chr = tmp[0];
			int start  = atoi(tmp[1].c_str());
			int end    = atoi(tmp[2].c_str());

			float gc = computeGC(chr, start, end, ref);
			bed_t region;
			region.chr = chr;
			region.start = start;
			region.end = end;
			region.gc = gc;
			if (tmp.size() > 3) {
				string info = tmp[3];
				region.info = info;
			}
			RegionsVector.push_back(region);
		}
	}
}


vector<string> inline split( std::string const& original, char separator )
{
    std::vector<std::string> results;
    std::string::const_iterator start = original.begin();
    std::string::const_iterator end = original.end();
    string::const_iterator next = std::find( start, end, separator );
    while ( next != end ) {
        results.push_back( string( start, next ) );
        start = next + 1;
        next = std::find( start, end, separator );
    }
    results.push_back( std::string( start, next ) );
    return results;
}

string getFileName(const string& s) {

   char sep = '/';
   size_t i = s.rfind(sep, s.length());
   if (i != string::npos) {
      return(s.substr(i+1, s.length() - i));
   }

   return("");
}
string IntToString (int& a){
	ostringstream temp;
	temp<<a;
	return temp.str();
}


float computeGC ( std::string chr, int start, int end, RefGenome& reference ) {

	std::string sequence = reference.QueryRegion(chr, start, end);
	int count = 0;
	for (auto& ntd : sequence) {
		if (ntd == 'C' || ntd == 'c' || ntd == 'G' || ntd == 'g' ) {
			count++;
		}
	}
	float gc_content = (count/ (float) sequence.length()) * (float) 100;
	return gc_content;
}




