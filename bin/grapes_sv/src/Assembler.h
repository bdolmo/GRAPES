#include <iostream>
#include <string> // carreguem la llibreria d'strings
#include <fstream> // llibreria per obrir inputs
#include <vector>
#include <list>
#include <algorithm>
#include <map>
#include <unordered_map>
extern bool debug;

struct overlapData {

	std::string seq1;
	std::string seq2;
	std::string contigSeq;
	int matches;
	int mismatches;
	int gaps;
	bool isRevComp;
};

struct Nucleotide {
	int A = 0;
	int C = 0;
	int T = 0;
	int G = 0;
	int totalBases = 0;
};
struct Contig {
	std::string seq;
	std::vector<Nucleotide> Consensus;

	bool operator ==(const Contig& lhs)
	{
		return seq == lhs.seq;
	}
};
struct candidateRead {
	std::string Sequence;
	int offset;
	std::string Seed;
	int mismatches;
	std::vector<Nucleotide> Consensus;
	Contig ContigStruct;
};

 class Assembler {

	// Mètodes
	public:
		// Constructor
		Assembler(std::vector<std::string>&, int, int, int, int);

		std::multimap<std::string, std::string> createReadTable();
		std::unordered_multimap<std::string, Contig> createPrefixTable();
	    std::pair<int, Contig> Extend( std::unordered_multimap<std::string, Contig>&,std::multimap<std::string, std::string>&, std::map<std::string, int>&, Contig&, int, bool);

		std::vector<std::string> ungappedGreedy();

		std::vector<std::string> gappedGreedy();
		overlapData findBestOverlap();

		// For testing purposes only!, needs lots of work 
 		std::vector<std::string> deBruijnGraph();

		int getTotalReads();

		int getNumAssembled();

	private:
		std::vector<std::string> reads;
		int minOlapLength;
		int maxMismatches;
		int maxGaps;
		int kSize;

		int totalReads;
		int totalAssembled;
		std::vector<Contig> ContigList;

	
 };
