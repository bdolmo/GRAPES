#ifndef SAMPLE_METRICS_T_H
#define SAMPLE_METRICS_T_H

#include <string>

 struct sample_metrics_t  {

	std::string name;
	int total_reads;
	int total_readsX;
	double mean_coverage;
	double median_coverage;		
	
	double mean_coverageX;
	double median_coverageX;

	double mean_countsX;
	double median_countsX;

	double mean_counts;
	double median_counts;

	double mean_insertSize;
	double sd_insertSize;		
	
	bool operator<(const sample_metrics_t& x) const {
		return (name < x.name);
	}
	inline bool comp(const sample_metrics_t& lhs, const sample_metrics_t& rhs){
	  return std::tie(lhs.name, lhs.total_reads,  lhs.total_readsX, lhs.mean_coverage, lhs.mean_coverageX, lhs.median_coverageX, lhs.mean_countsX, lhs.median_countsX, lhs.median_coverage, lhs.mean_counts, lhs.median_counts, lhs.mean_insertSize, lhs.sd_insertSize) < std::tie(rhs.name, rhs.total_reads, rhs.total_readsX, rhs.mean_coverage, rhs.mean_coverageX, rhs.median_coverageX, rhs.mean_countsX, rhs.median_countsX, rhs.median_coverage, rhs.mean_counts, rhs.median_counts, rhs.mean_insertSize, rhs.sd_insertSize);
	}
 };


#endif
