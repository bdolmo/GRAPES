#ifndef GENOMICRANGE_T_H
#define GENOMICRANGE_T_H

#include <iostream>
#include <string>


struct GenomicRange {
	std::string chromosomeA;
	std::string chromosomeB;
	int start;
	int end;
	int supporting;
	int t_cov;
	std::string softclip_type;
	double discordant_ratio;
	float pvalue_discordant;
	int cumulativeSize;
	int copy_number;
	int mapq;

	int chromosome_start;
	int chromosome_end;
	int chr_centromere_start;
	int chr_centromere_end;

	bool operator<(const GenomicRange& x) const {
		return (start < x.start);
	}
	// lexicographical comparison provides strict weak ordering
	inline bool comp_gr(const GenomicRange& lhs, const GenomicRange& rhs){
	  return std::tie(lhs.chromosomeA, lhs.chromosomeB, lhs.start, lhs.end, lhs.supporting, lhs.t_cov, lhs.softclip_type, lhs.discordant_ratio, lhs.pvalue_discordant, lhs.cumulativeSize, lhs.copy_number, lhs.mapq, lhs.chromosome_start, lhs.chromosome_end, lhs.chr_centromere_start, lhs.chr_centromere_end) < std::tie(rhs.chromosomeA, rhs.chromosomeB, rhs.start, rhs.end, rhs.supporting, rhs.t_cov, rhs.softclip_type, rhs.discordant_ratio, rhs.pvalue_discordant, rhs.cumulativeSize, rhs.copy_number, rhs.mapq, rhs.chromosome_start, rhs.chromosome_end, rhs.chr_centromere_start, rhs.chr_centromere_end);
	}
};

#endif
