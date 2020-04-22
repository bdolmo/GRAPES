#include <iostream>
#include <string.h>
#include <fstream>
#include <vector>
#include <map>
#include <sstream>
#include <algorithm>
#include <unordered_map>
#include <ctime>
#include <cstdlib>
#include <math.h>
#include <thread>
#include <omp.h>

// Boost library
#include <boost/program_options.hpp>
#include <boost/iterator/transform_iterator.hpp>
#include <boost/algorithm/string/replace.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/filesystem.hpp>
#include <boost/math/distributions.hpp>
#include <boost/progress.hpp>

// smith-waterman SIMD
#include "ssw_cpp.h"

// grapes headers
#include "Assembler.h"
#include "alignContig.h"
#include "Aligner.h"
#include "cnv.h"
#include "excludeRegion.h"
#include "clusterDiscordant.h"
#include "clusterSC.h"
#include "clust_t.h"
#include "dna2bit.h"
#include "GenomicRange.h"
#include "largeSV.h"
#include "pe_adjacency.h"
#include "ReadPairInfo.h"
//#include "reAlignSoftClip.h"
#include "sam_t.h"
#include "smallSV.h"
#include "SV.h"
#include "Trimmer.h"
#include "Utils.cpp"
#include "vcf_t.h"
#include "varCall.h"
#include "hashAligner.h"
#include "seed_t.h"
#include "callSV.h"

#include "SeqLib/RefGenome.h"
#include "SeqLib/BamReader.h"
#include "SeqLib/BWAWrapper.h"
#include "SeqLib/BamRecord.h"
#include "SeqLib/BFC.h"
#include "SeqLib/FermiAssembler.h"


using namespace std;
using namespace SeqLib;
namespace po = boost::program_options;
namespace fs = boost::filesystem;

bool debug = false;

// Declaring prototypes
//vector<int> returnAlignment (const StripedSmithWaterman::Alignment& alignment);
std::pair<int, int> adjustMinimumBreakReads(float&, int&, int&);
void writeRawCalls( std::vector<vcf_t>&, long&, std::string&, std::string&, std::string&, const int&); 

bool is_file_exist(const char *);
void resolveLargeSV ( string&, map<string, vector<GenomicRange>>&, RefGenome&, std::string, std::map<std::string, std::vector<GenomicRange>>&, std::map<std::string, std::vector<GenomicRange>>&, std::map<std::string, std::vector<GenomicRange>>&, std::map<std::string, std::vector<GenomicRange>>&, std::map<std::string, std::vector<GenomicRange>>&, map<string, vector<sam_t>>&, std::vector<vcf_t>&, std::string);
void locateBreakReads( std::string&, std::map<std::string, std::vector<GenomicRange>>&, std::map<std::string, std::vector<GenomicRange>>&, map<string, vector<sam_t>>&, RefGenome&, std::string&, std::string, std::vector<vcf_t>&,  std::string, std::string);
void findSmallSV ( string&, map<string, vector<sam_t>>&, int &, RefGenome&, std::string, std::string, std::vector<vcf_t>&, std::string);
void Analyze ( string&, vector<string>&, string, int, string, int, int, string&, std::string, std::vector<vcf_t>&, std::string, double);
std::tuple<double, double> reciprocalOverlap ( int&, int&, int&, int& );
string fillQualString ( std::string& );

// lexicographical comparison provides strict weak ordering
inline bool comp_gr(const GenomicRange& lhs, const GenomicRange& rhs){
	return std::tie(lhs.chromosomeA, lhs.chromosomeB, lhs.start, lhs.end, lhs.supporting, lhs.t_cov, lhs.softclip_type, lhs.discordant_ratio) < std::tie(rhs.chromosomeA, rhs.chromosomeB, rhs.start, rhs.end, rhs.supporting, rhs.t_cov, rhs.softclip_type, rhs.discordant_ratio);
}
inline bool comp_vcf(const vcf_t& lhs, const vcf_t& rhs){
        return std::tie(lhs.chr, lhs.start, lhs.end, lhs.precision, lhs.svtype, lhs.breakReads, lhs.assembled, lhs.covPvalue, lhs.discordants, lhs.alleleBalance) < std::tie(rhs.chr, rhs.start, rhs.end, rhs.precision, 		rhs.svtype, rhs.breakReads, rhs.assembled, rhs.covPvalue, rhs.discordants, rhs.alleleBalance);
}
inline bool comp_sam(const sam_t& lhs, const sam_t& rhs){
  return std::tie(lhs.align_pos) < std::tie(rhs.align_pos);
}

callSV::callSV(std::string& _bamFile , std::string& _reference) {
	bamFile = _bamFile;
	reference = _reference;
}

hashAligner::hashAligner(int _kSize, int _minSeeds, int _maxMismatches) {
	kSize    = _kSize;
	minSeeds  = _minSeeds;
	maxMismatches = _maxMismatches;
}

// Inicialitzem el constructor per la classe Aligner
Aligner::Aligner(const string& _query, const string& _reference, int _max_mismatches) {
	query          = _query;
	reference      = _reference;
	max_mismatches = _max_mismatches;
}

// Inicialitzem el constructor per la classe Aligner
alignContig::alignContig( std::vector<std::string>& _contigs, std::string& _reference, std::string& _chr, int& _genomicStart, int& _genomicEnd, int& _readLength, std::string& _genome) {
	contigs        = _contigs;
	reference      = _reference;
	chr            = _chr;
	genomicStart   = _genomicStart;
	genomicEnd     = _genomicEnd;
	readLength     = _readLength;
	genome	       = _genome;
}

Assembler::Assembler( std::vector<std::string>& _reads, int _minOlapLength, int _maxMismatches, int _maxGaps, int _kSize) {
	reads         = _reads;
	minOlapLength = _minOlapLength;
	maxMismatches = _maxMismatches;
	maxGaps       = _maxGaps;
	kSize         = _kSize;
}

excludeRegion::excludeRegion ( const std::string& _regionsExclude) {
	regionsExclude = _regionsExclude;
} 

dna2bit::dna2bit( std::string& _dna_str ) {
	dna_str = _dna_str;
}

clusterDiscordant::clusterDiscordant( const std::string& _bamFile, const std::string& _outdir, const std::string& _sampleName, int _minClusterSize, int _numSDs, long& _genomeSize, const std::string& _genome, std::string _getCountsByWindow, int& _maxInsertSize ) {
	bamFile        = _bamFile;
	outdir         = _outdir;
	sampleName     = _sampleName;
	minClusterSize = _minClusterSize;
	numSDs         = _numSDs;
	genomeSize     = _genomeSize;
	genome         = _genome;
	getCountsByWindow = _getCountsByWindow;
	maxInsertSize  = _maxInsertSize;
}

Cnv::Cnv( std::string& _outDir, std::string& _raw_counts, int _medianSomaticCounts, int _medianGerminalCounts, int _n_states ) {
	outDir                = _outDir,
	raw_counts            = _raw_counts;
	medianSomaticCounts   = _medianSomaticCounts;
	medianGerminalCounts  = _medianGerminalCounts;
	n_states              = _n_states;
}

clusterSC::clusterSC( const string& _bamFile, const string& _outDir, int _minSoftLength, float _softClipLenPerc, int _minSoftClusterSize, int _totalSR) {
	bamFile           = _bamFile;
	outDir            = _outDir;
	minSoftLength     = _minSoftLength;
	softClipLenPerc   = _softClipLenPerc;
	minSoftClusterSize= _minSoftClusterSize;
	totalSR           = _totalSR;
}

// Inicialitzem el constructor per la classe largeSV
largeSV::largeSV(std::string& _bamFile, std::string& _reference, int _minOlapBases, int _maxMismatches) {
	bamFile   = _bamFile;
	reference = _reference;
	minOlapBases = _minOlapBases;
	maxMismatches= _maxMismatches;
}

// Inicialitzem el constructor per la classe Trimmer
Trimmer::Trimmer(string _input_seq, string _qualities) {
	input_sequence = _input_seq;
	qualities      = _qualities;	
}

// Inicialitzem el constructor per la classe smallSV
smallSV::smallSV(std::string& _bamFile, std::string& _reference, int _minOlapBases, int _maxMismatches, int _readLength) {
	bamFile   = _bamFile;
	reference = _reference;
	minOlapBases = _minOlapBases;
	maxMismatches= _maxMismatches;
	readLength   = _readLength;
}

// Inicialitzem el constructor per la classe SV
SV::SV(vector<sam_t>& _reads, string _chrA, string _chrB, int _posA, int _posB, string _soft_type, int _total_reads, int _total_assembled, string _bam, string& _ref_sequence, std::string _svlength, int _nDiscordants, double _pvalue_discordant, double _kmer_diversity) {
	reads             = _reads;
	chrA              = _chrA;
	chrB              = _chrB;
	posA              = _posA;
	posB              = _posB;
	softclip_type     = _soft_type;
	total_reads       = _total_reads;
	total_assembled   = _total_assembled;
	bamFile           = _bam;
	ref_sequence      = _ref_sequence;
	svlength          = _svlength;
	//vcf_v             = _vcf_v;
	nDiscordants      = _nDiscordants;
	pvalue_discordant = _pvalue_discordant;
	kmer_diversity    = _kmer_diversity;
}

// Inicialitzem el constructor per la classe reAlignSoftClip
/*reAlignSoftClip::reAlignSoftClip(std::string _read, string _chrom, int _pos, int _end,  std::string _cigar, std::string _ref_sequence, std::string _soft_type, std::string _VCF, std::string _svlength, int _posA, int _posB, int _mapq) {
	read = _read;
	chr  = _chrom;
	pos   = _pos;
	end = _end;
	cigar = _cigar;
	ref_sequence = _ref_sequence;
	softclip_type = _soft_type;
	VCF = _VCF;
        svlength = _svlength;
	posA = _posA;
	posB = _posB;
	mapq = _mapq;
}*/

varCall::varCall( std::string _bamFile, std::string _genome, int _minMapQ, int _minBaseQ, int _minSNV, int _minCov, float _minHomRatio ) {
	bamFile    = _bamFile;
	genome     = _genome;
	minMapQ    = _minMapQ;
	minBaseQ   = _minBaseQ;
	minSNV     = _minSNV;
	minCov     = _minCov;
	minHomRatio= _minHomRatio;
}

int totalFR, totalRF, totalRR, totalFF, totalSR, readLength, minOlapBases, maxMismatches;
std::vector<vcf_t> vcf_v;

//####### MAIN #######
int main (int argc, char* argv[]) {

	string reference, bamFile, output, outdir, regionsExclude;
	string wes = "off";
	string wgs = "on";
	string cnv = "on";
	string fs =  "on";
	string fl =  "on";
	string fc =  "on";
	string exclude;

	int nDiscordants;
	int nBreakReads;
	int nSds;
	int total_seq_trimmed = 0;
	int minSize;
	int maxSize;
	float softClipLenPerc;
	int threads;
	int minOlapPerc;
	int maxMismatch;
	std::map<std::string, std::vector<sam_t>> breakReadClusters;

	std::map<std::string, std::vector<GenomicRange>> FR_ranges;
	std::map<std::string, std::vector<GenomicRange>> RF_ranges;
	std::map<std::string, std::vector<GenomicRange>> FF_ranges;
	std::map<std::string, std::vector<GenomicRange>> RR_ranges;
	std::map<std::string, std::vector<GenomicRange>> INT_ranges;
	double meanDiscordants;

	char* dt;
	time_t now;	
	po::options_description description("\n Program: Grapes SV Version 0.9.3 \n Contact: bdelolmo@gencardio.com \n Usage");

	description.add_options()
	    ("bam,b", po::value<std::string>(&bamFile)->required(), "BAM file")
	    ("genome,g", po::value<std::string>(&reference)->required(), "Genome reference FASTA")
	    ("outname,n", po::value<std::string>(&output)->required(), "Output file name")
	    ("outdir,o", po::value<std::string>(&outdir)->required(), "Output directory")
	    ("wgs", po::value<std::string>(&wgs)->default_value("on"), "Whole-genome analysis")
	    ("wes", po::value<std::string>(&wes)->default_value("off"), "Whole-exome/Targeted panel analysis")
	    ("clust,c", po::value<int>(&nDiscordants)->default_value(2), "Number of discordant pairs to create a cluster")
	    ("sds,s", po::value<int>(&nSds)->default_value(3), "Number of standard deviations from the man insert size")
	    ("breads,r", po::value<int>(&nBreakReads)->default_value(2), "Number of breakreads")
	    ("sclen,l", po::value<float>(&softClipLenPerc)->default_value(10), "Minimum sofclipped bases (%) considered")
	    ("minOlap,a", po::value<int>(&minOlapPerc)->default_value(13), "Minimum overlapping bases (%) to trigger targeted assembly")
	    ("maxMismatch,m", po::value<int>(&maxMismatch)->default_value(4), "Maximum number of mismatches (%) allowed to trigger targeted assembly")
	    ("min-size,i", po::value<int>(&minSize)->default_value(50), "Minimum SV size (in bp)")
		("max-size,j", po::value<int>(&maxSize)->default_value(5000000), "Maximum SV size (in bp)")

	    ("find-cnv", po::value<std::string>(&cnv)->default_value("on"), "Find CNV through read depth analysis")
	    ("find-small", po::value<std::string>(&fs)->default_value("on"), "Find small SVs")
	    ("find-large", po::value<std::string>(&fl)->default_value("on"), "Find large SVs")
	    ("exclude,e",  po::value<std::string>(&regionsExclude)->required(), "Exclude regions,centromers, gaps (BED)")
	    ("threads,t", po::value<int>(&threads)->default_value(1), "Number of CPU cores")
		("debug,d", "Debug")
	    ("version,v", "Display version number")
	    ("help,h", "Display help message\n");

	po::variables_map vm;
	po::positional_options_description input_options;
	po::store(po::parse_command_line(argc, argv, description), vm);
	po::store(po::command_line_parser(argc, argv).options(description).positional(input_options).run(), vm);

	try {
		if(vm.count("help")){
			std::cout << description;
			return 0;
		}
		if(vm.count("version")){
			std::cout << "Program: GRAPES - SV detection on NGS data" << endl;
			std::cout << "version: v0.9.2" << endl;
			return 0;
		}
		if(vm.count("debug")) {
			debug = true;
		}
	    po::notify(vm);
		if(vm.count("bam")){
			now = time(0);
			// convert now to string form
			dt = ctime(&now);

			std::cout << endl << " INFO: GRAPES SV analysis - " << dt << endl;
			std::cout << " ---------- Input Parameters ----------" << endl;
			std::cout << " INFO: bamFile input: " << bamFile << endl;
		}
		if (wes == "on") {
			cnv == "off";
		}
		if(vm.count("genome")){
			std::cout << " INFO: Genome input: " << reference << endl;
		}
		if (vm.count("find-small") || vm.count("find-large") || vm.count("find-complex") || vm.count("find-cnv")) {
			if(cnv != "on" && cnv != "off") {
				cout << " ERROR: please provide 'on' or 'off' to find-cnv option "<<  endl;
				return 0;
			}
			std::cout << " INFO: Find CNVs: " << fl << endl;
			if(fs != "on" && fs != "off") {
				cout << " ERROR: please provide 'on' or 'off' to find-small option "<<  endl;
				return 0;
			}
			std::cout << " INFO: Find Small SVs: " << fs << endl;		
			if(fl != "on" && fl != "off") {
				cout << " ERROR: please provide 'on' or 'off' to find-large option" << endl;
				return 0;
			}
			std::cout << " INFO: Find Large SVs: " << fl << endl;
			if(fc != "on" && fc != "off") {
				cout << " ERROR: please provide 'on' or 'off' to find-complex option" << endl;
				return 0;
			}
			std::cout << " ---------------------------------------" << endl;		
		}

	} catch (po::error& e) {
		std::cout << description;
 		return 1;
	}
	BamReader br;
	if (br.Open(bamFile)) {
	}
	else {
		cout << " ERROR: please provide a valid input bamFile file" << endl;
		return 0;
	}

	RefGenome ref;
	ref.LoadIndex(reference);
	if (ref.LoadIndex(reference)) {
	}
	else {
		cout << " ERROR: please provide a valid Genome file in FASTA format" << endl;
		return 0;
	} 
	fs::path full_path( fs::initial_path<fs::path>() );
	full_path = fs::system_complete( fs::path( argv[0] ) );
	fs::path exe_path_boost =  full_path.parent_path();

	long genomeSize = 2966866909;

	std::string logF = outdir + "/" + "grapes.log";
	std::ofstream logFile;
	logFile.open(logF);
	//cout << outdir << "\t" << logF << endl;

	logFile << endl << " INFO: GRAPES SV analysis - " << dt << endl;
	logFile << " ---------- Input Parameters ----------" << endl;
	logFile << " INFO: bamFile input: " << bamFile << endl;

	std::cout << "\n INFO: Analyzing " << bamFile << endl;
	logFile << "\n INFO: Analyzing " << bamFile << endl;

	std::string getCountsByWindow = "on"; 

	if (wes == "on") {
		getCountsByWindow = "off"; 
	}

	omp_set_num_threads(threads);
	omp_set_nested(1);

	//std::string contig = "GCAGAAAAGCTGAAATTTCTAAAAATCAGAGCAACTCTTCTCCTCCAAAGGAACGCAGCTCCTCACCAGCAACGGAACAAAGCTGTAATAACAAACTTCTCTGAGCTAAAGGAGGATGTTCGAAC";
	//std::string refSequence = "aagagagtagtggttctcccagaatggagtttgagatctgagaacggacagactgcctcctcaagtgggtccctgacccctgagtagcctaactgcgagacacctcccagtaggggccgactgacacctcacacagccaggtgcccctctgagatgaagcttccagagaaaggatcaggcaggaacatttgccgttctgcaatatttgcggttctgcagcctctgctggtgatacccaggaaacagggtctggagtggacctccagcaaactccaacagacctgcagctgaggggcctgactgttagaaggaaaactaacaaacagaaaggacatccacaccaaaaccccatctgcacgtcaccatcatcaaagaccaaaggtagataaaaccacaatgatggggagaaaacagagcagaaaagctgaaaattctaaaaatcagagcaactcttctcctccaaaggaacgcagctcctcaccagcaacggaacaaagctggacggagaatgactttgacgagttgagagaagaaggcttcagacaatcagtaataacaaacttctctgagctaaaggaggatgttcgaacccatcgcaaagaagctaaaaaccttgaaaaaagattagacaaatggctaactagaataaacagcatagagaagatgttaaatgacctgatggagctgaaaaccatggcacgagaactacatgatgcatgcacaagcttcagtagccaattcgatcaagtggaagaaagggtatcagtgattgaagatcaaatgaatgaaatgaagcgagaagagaagtttagagaaaaaagattaaaaagaaacgaacaaagcctccaagaaatatgggactatgtgaaaagaccaaatctacatctgattggtgtacctgaaagtgacggggagaatggaaccaagttggaaaacactctgcaggatatcatccaggagaacttccccaacctagcaaggcaggccaacattcaaattca";

	//std::vector<string> contigVec;
	//contigVec.push_back(contig);

	/*int Pstart = 10000;
	int Pend   = 10000;
	string CHR = "chr1";
	alignContig A (contigVec, refSequence, CHR, Pstart, Pend);
	A.Align(vcf_v);

	cout << vcf_v.size() << endl;

	for (auto l : vcf_v) {

cout << l.chr << "\t" << l.start << "\t" << l.end << "\t" << l.precision << ";SVTYPE=" << l.svtype << ";MAPQ=" << l.mapq << ";KDIV=" << l.kdiv << ";BREAKREADS=" << l.breakReads << ";ASSEMBLED=" << l.assembled << ";PE=" <<
	l.discordants << ";CSDISC=" << l.cumulativeSize << ";NINS=" << 0 << ";RDratio=" << l.RDratio << ";RDmad=" << l.RDmad << ";RDsupp=" << l.hasRDsupport << ";LOHsupp=" << l.LOHsupport << "\n";

	}
	*/
	//return 0;
	
	//hashAligner H (14, 2, 1);
	//std::multimap<std::string, int> queryHash = H.indexContig(contig);
	

	//std::vector<seed_t> seedV = H.findSeeds(RefSequence, contig, queryHash);

	//std::vector<seed_t> mseedV = H.mergeIntervals(seedV);
    	//H.createDAG(mseedV);

	//return 0;

	// Exclude region
	excludeRegion blackList(regionsExclude);
	map<string, vector<GenomicRange>> mapOfRegions = blackList.readAndSave();

	now = time(0);
	dt = ctime(&now);
	logFile << " INFO: Extracting discordant/split reads and read-depth - " << dt << endl;

	clusterDiscordant clustObj (bamFile, outdir, output, nDiscordants, nSds, genomeSize, reference, getCountsByWindow, maxSize );
	clustObj.extractDiscordant(mapOfRegions);
	std::vector<int> readCounts_somatic  = clustObj.getSomaticCounts();
	std::vector<int> readCounts_germinal = clustObj.getGerminalCounts();
	int binSize = clustObj.getBinSize();
	float mean_coverage = clustObj.getMeanCoverage();
	
	if (wes == "off") {
	    std::tie(nDiscordants, nBreakReads) = adjustMinimumBreakReads(mean_coverage, nDiscordants, nBreakReads);
	}

	// Deprecated Jul19. Now using reference coverage profile https://doi.org/10.3389/fgene.2015.00045
	/*if (cnv == "on" && wes == "off") {

		std::string outCNV = outdir + "/" + output + ".CNV.bed";

		int medianSomaticCounts;
		if (readCounts_somatic.size() == 0) {
			medianSomaticCounts = 0;
		} 
		else {
			medianSomaticCounts = computeMedian(readCounts_somatic);
		}
		int medianGerminalCounts;
		if (readCounts_germinal.size() == 0) {
			medianGerminalCounts = 0;
		} 
		else {
			medianGerminalCounts = computeMedian(readCounts_germinal);
		}
		std::cout << " INFO: Median somatic chromosome counts per window: " << medianSomaticCounts << "\n";
		std::cout << " INFO: Median germinal chromosome counts per window: " << medianGerminalCounts << "\n";

		logFile << " INFO: Median somatic chromosome counts per window: - " << medianSomaticCounts << "\n";
		logFile << " INFO: Median germinal chromosome counts per window: - " << medianGerminalCounts << "\n";

		// Normalizing counts
		std::string count_file = outdir + "/" + "coverage_windows.bed";
		Cnv cnv1 (outdir, count_file, medianSomaticCounts, medianGerminalCounts, 3);

		if (!is_file_exist(outCNV.c_str())) {
			now = time(0);
			dt = ctime(&now);			
			std::cout << " INFO: Normalizing GC content " << endl;
			logFile << " INFO: Normalizing GC content " << dt << endl;
			cnv1.normalize_gc();

			now = time(0);
			dt = ctime(&now);
			std::cout << " INFO: Normalizing Mappability " << endl;
			logFile << " INFO: Normalizing Mappability " << dt << endl;

			cnv1.normalize_mappability();

			// Segmenting using external software (HMMseg)
			std::string cmd = "perl " +  exe_path_boost.string() + "/HMMseg.pl " + outdir + "/" + "copy_ratios.bed " + std::to_string(binSize) + " " + outCNV + " " + output + " " + outdir + " " + regionsExclude + " " + bamFile + " " + reference + " " + std::to_string(mean_coverage) ;

			system(cmd.c_str());
		}
		vcf_v = cnv1.dumpCNV(output, bamFile, reference);
	}*/
	//return 0;

	// Analyzing large SV	
	if (fl == "on") {
		now = time(0);
		dt = ctime(&now);
		logFile << " INFO: Clustering FR pairs " << dt << endl;
		std::string discBam = outdir + "/" + output + ".FR.bam";
		clustObj.cluster(discBam, FR_ranges);

		now = time(0);
		dt = ctime(&now);
		logFile << " INFO: Clustering RF pairs " << dt << endl;
		discBam = outdir + "/" + output + ".RF.bam";
		clustObj.cluster(discBam, RF_ranges);

		now = time(0);
		dt = ctime(&now);
		logFile << " INFO: Clustering FF pairs " << dt << endl;
		discBam = outdir + "/" + output + ".FF.bam";
		clustObj.cluster(discBam, FF_ranges);

		now = time(0);
		dt = ctime(&now);
		logFile << " INFO: Clustering RR pairs " << dt << endl;
		discBam = outdir + "/" + output + ".RR.bam";
		clustObj.cluster(discBam, RR_ranges);
	}
	std::string srBam = outdir + "/" + output + ".SR.bam";

	totalFR = clustObj.getTotalFR();
	totalRF = clustObj.getTotalRF();
	totalFF = clustObj.getTotalFF();
	totalRR = clustObj.getTotalRR();

	totalSR = clustObj.getTotalSR();

	std::map<std::string, std::vector<GenomicRange>> mapBreakRanges;

	// BreakRead clustering
	now = time(0);
	dt = ctime(&now);
	logFile << " INFO: Clustering Split-Reads " << dt << endl;
	clusterSC srClust (srBam, outdir, 10, softClipLenPerc, nBreakReads, totalSR);
	srClust.extractSC( breakReadClusters, mapBreakRanges, mapOfRegions, reference );

	readLength = srClust.getReadLength();
	minOlapBases  = readLength * minOlapPerc/100;
	maxMismatches = readLength * maxMismatch/100;

	std::string dInfo = outdir + "/" + output + ".discordantInfo.txt";
	std::ofstream drpInfo;

	drpInfo.open(dInfo);
	drpInfo << totalFR << "\t" << totalRF << "\t" << totalFF << "\t" << totalRR << endl;
	drpInfo.close();

	// Analyzing large SV	
	if (fl == "on") {

		largeSV lsv (bamFile, reference, minOlapBases, maxMismatches );
		
		now = time(0);
		dt = ctime(&now);
		logFile << " INFO: Searching large deletions " << dt << endl;
		std::cout << " INFO: Searching large deletions "<< endl;
		lsv.callStructVar( FR_ranges, mapBreakRanges, breakReadClusters, vcf_v, "DEL", "undef");

		now = time(0);
		dt = ctime(&now);
		logFile << " INFO: Searching large duplications " << dt << endl;
		std::cout << endl << " INFO: Searching large duplications "<< endl;
		lsv.callStructVar( RF_ranges, mapBreakRanges, breakReadClusters, vcf_v, "DUP", "undef");

		now = time(0);
		dt = ctime(&now);
		logFile << " INFO: Searching large inversions (FF) " << dt << endl;
		std::cout << endl << " INFO: Searching large inversions (FF) "<< endl;
		lsv.callStructVar( FF_ranges, mapBreakRanges, breakReadClusters, vcf_v, "INV", "RIGHT");

		now = time(0);
		dt = ctime(&now);
		logFile << " INFO: Searching large inversions (RR) " << dt << endl;
		std::cout << endl << " INFO: Searching large inversions (RR) "<< endl;
		lsv.callStructVar( RR_ranges, mapBreakRanges, breakReadClusters, vcf_v, "INV", "LEFT");
	}

	std::string fastq = outdir + "/" + output + ".fastq";
	std::string VCF   = outdir + "/" + output + ".vcf";

	// Analyzing small SV	
	if (fs == "on") {
		now = time(0);
		dt = ctime(&now);

		if (wes == "on") {
			logFile << " INFO: Assemblying break-reads " << dt << endl;
			std::cout << endl << " INFO: Assembling break-reads "<< endl;
			smallSV ssv (bamFile, reference, minOlapBases, maxMismatches, readLength );
			//ssv.callSmallSvExome (breakReadClusters, mapBreakRanges, fastq, vcf_v, wes, outdir);
			ssv.callSmallSV( breakReadClusters, mapBreakRanges, fastq, vcf_v, wes, outdir );
		}
		else {
			logFile << " INFO: Assemblying break-reads " << dt << endl;
			std::cout << endl << " INFO: Assembling break-reads "<< endl;
			smallSV ssv (bamFile, reference, minOlapBases, maxMismatches, readLength );
			ssv.callSmallSV( breakReadClusters, mapBreakRanges, fastq, vcf_v, wes, outdir );
		}
	}
	std::string vcfOut = outdir + "/" + output + ".tmp.rawcalls.bed";
	writeRawCalls(vcf_v, genomeSize, bamFile, reference, vcfOut, minSize);

	return 0;
}

//####### END #######
void writeRawCalls( std::vector<vcf_t>& vcf_v, long& genomeSize, std::string& bamFile, std::string& reference, std::string& vcfOut, const int& minSize) {

	// Coordinate sorting
	std::sort(vcf_v.begin(),vcf_v.end(), comp_vcf);

	std::ofstream vcf_out;
	vcf_out.open(vcfOut.c_str(), std::ofstream::app);

	int numInserts;
	string reciprocalSupport = "no";

	int N = vcf_v.size();
	int ntotal = 0;

	for (auto& l : vcf_v) {
		
		int size = l.end-l.start;
		if (size < minSize) {
			continue;
		}
		if (l.svtype == "DEL") {
			numInserts = totalFR;
		}
		if (l.svtype == "DUP") {
			numInserts = totalRF;
		}		
		if (l.svtype == "INV") {
			numInserts = totalFF + totalRR;
		}

vcf_out << l.chr << "\t" << l.start << "\t" << l.end << "\t" << l.precision << ";SVTYPE=" << l.svtype << ";MAPQ=" << l.mapq << ";KDIV=" << l.kdiv << ";BREAKREADS=" << l.breakReads << ";ASSEMBLED=" << l.assembled << ";PE=" <<
	l.discordants << ";CSDISC=" << l.cumulativeSize << ";NINS=" << numInserts << ";RDratio=" << l.RDratio << ";RDmad=" << l.RDmad << ";RDsupp=" << l.hasRDsupport << ";LOHsupp=" << l.LOHsupport << "\n";
	}
}

std::pair<int, int> adjustMinimumBreakReads(float& meanCoverage, int& nDiscordants, int& nBreakReads) {

	if (meanCoverage <= 10) {
		nDiscordants = 1;
		nBreakReads  = 1;
	}
	if (meanCoverage > 10 && meanCoverage<=20) {
		nDiscordants = 2;
		nBreakReads  = 2;
	}
	if (meanCoverage > 20 && meanCoverage<=30) {
		nDiscordants = 3;
		nBreakReads  = 3;
	}
	if (meanCoverage > 30 && meanCoverage<=40) {
		nDiscordants = 4;
		nBreakReads  = 3;
	}
	if (meanCoverage > 40 && meanCoverage<=50) {
		nDiscordants = 5;
		nBreakReads  = 3;
	}
	if (meanCoverage > 50 && meanCoverage<=60) {
		nDiscordants = 6;
		nBreakReads  = 3;
	}
	if (meanCoverage > 60 && meanCoverage<=70) {
		nDiscordants = 7;
		nBreakReads  = 3;
	}
	if (meanCoverage > 70 && meanCoverage<=80) {
		nDiscordants = 8;
		nBreakReads  = 3;
	}
	if (meanCoverage > 80 && meanCoverage<=90) {
		nDiscordants = 9;
		nBreakReads  = 3;
	}
	if (meanCoverage > 100) {
		nDiscordants = 10;
		nBreakReads  = 3;
	}
	return std::make_pair(nDiscordants, nBreakReads);										
}






