#ifndef READPAIRINFO_H
#define READPAIRINFO_H

#include <iostream>
#include <string>


// This structure contains information from both reads of a pair

struct ReadPairInfo {

	std::string readName;
	bool																																																																																																																																																																																																													 is_chimeric;	
	// 1st mate
	std::string first_mate_chr;
	int first_mate_start;
	int first_mate_mapq;
	std::string first_mate_cigar;
	int first_mate_alignPos;
	int first_mate_alignEnd;
	std::string first_mate_strand;
	

	// 2nd mate
	std::string second_mate_chr;
	int second_mate_start;
	int second_mate_mapq;
	std::string second_mate_cigar;
	int second_mate_alignPos;
	int second_mate_alignEnd;
	std::string second_mate_strand;
};

#endif
