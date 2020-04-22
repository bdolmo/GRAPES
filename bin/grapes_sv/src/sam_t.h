#ifndef SAM_T_H
#define SAM_T_H


#include <iostream>
#include <string>

//std::string IntToString ( int& );

 struct sam_t {

	std::string read_name;
	std::string chr;
	int pos;
	int align_pos;
	int align_end;
	std::string strand;
	std::string order;
	std::string cigar;
	int mapq;
	std::string seq;
	int seqLen;

	uint8_t* bitSeq;
	std::string qual;
	int clipped;
	std::string coordinate;
	int numMatched;
	std::string sctype;
	bool isRevComp;
		
	bool operator<(const sam_t& x) const {
		return (read_name < x.read_name);
	}
	inline bool comp(const sam_t& lhs, const sam_t& rhs){
	  return std::tie(lhs.align_pos) < std::tie(rhs.align_pos);
	}

 };


#endif
