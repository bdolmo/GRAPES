#ifndef VCF_T_H
#define VCF_T_H

#include <iostream>
#include <string>


 struct vcf_t {

	std::string chr;
	int start;
	int end;
	std::string precision;
	std::string svtype;
	int breakReads;
	int assembled;
	double covPvalue;
	int discordants;
	double alleleBalance;
	
	// Parameters only for discordant pair clusters
	int cumulativeSize;
	float discPvalue;
	int mapq;
	double kdiv;
	bool hasRDsupport;
	std::string RDratio;
	std::string RDmad;
	std::string LOHsupport;
		
	bool operator<(const vcf_t& x) const {
		return (chr < x.chr);
	}
	inline bool comp_vcf(const vcf_t& lhs, const vcf_t& rhs){
	  return std::tie(lhs.chr, lhs.start, lhs.end, lhs.precision, lhs.svtype, lhs.breakReads, lhs.assembled, lhs.covPvalue, lhs.discordants, lhs.alleleBalance, lhs.cumulativeSize, lhs.discPvalue, lhs.mapq, lhs.kdiv, lhs.hasRDsupport, lhs.RDratio, lhs.RDmad, lhs.LOHsupport) < std::tie(rhs.chr, rhs.start, rhs.end, rhs.precision, rhs.svtype, rhs.breakReads, rhs.assembled, rhs.covPvalue, rhs.discordants, rhs.alleleBalance, rhs.cumulativeSize, rhs.discPvalue, rhs.mapq, rhs.kdiv, rhs.hasRDsupport, rhs.RDratio, rhs.RDmad, rhs.LOHsupport);
	}
 };


#endif
