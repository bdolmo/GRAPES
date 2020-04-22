#include <iostream>
#include <vector>
#include <string> 
#include <fstream>
#include "vcf_t.h"

extern bool debug;

struct  alignmentStruct {

	// Positions are 0-based!
	int score;
	int refStart;
	int refEnd;
	int queryStart;
	int queryEnd;
	int queryLength;
	std::string querySeq;
	std::string strand;

	int mismatches;
	int length;
	int refLength;
	std::string cigar;
	std::string refSeq;

	int contigStart;
	int contigEnd;
	int numClip;
	int numMatched;
	int multimaps;
	int repeatsOnRead;
	bool operator<(const alignmentStruct& x) const {
		return (contigStart < x.contigStart);
   	}
};


class alignContig {

	public:
		// Constructor
		alignContig(std::vector<std::string>&, std::string&, std::string&, int&, int&, int&, std::string&);

		std::vector<vcf_t> Align(int&, int&, double&, int&);

		std::vector<vcf_t> detectSV(std::vector<alignmentStruct>&, int&, int&, double&, int&, int&, std::string&, std::string&);

	private:
		std::vector<std::string> contigs;
		std::string reference;
		std::string chr;
		int genomicStart;
		int genomicEnd;
		int readLength;
		std::string genome;
};

