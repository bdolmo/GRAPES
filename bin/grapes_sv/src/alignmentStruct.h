#include <iostream>
#include <string>


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

	bool operator<(const alignmentStruct& x) const {
		return (contigStart < x.contigStart);
   	}
};

