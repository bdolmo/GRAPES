#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <fstream>

#include "cnv.h"
#include "vcf_t.h"
#include "Utils.cpp"
#include "varCall.h"

int median (vector<int>&);
void printProgBar( int);

std::vector<vcf_t> Cnv::dumpCNV( std::string outName, std::string& bamFile, std::string& reference) {

   std::string outCNV = outDir + "/" + outName + ".CNV.bed";

   fstream cnv_vcf;
   cnv_vcf.open (outCNV);
   std::string line;

   std::vector<vcf_t> vcf_v;

   //std::cout << " INFO: Annotating LOH events across deletions" << std::endl; 	

   int totalCNV = 0;
   if (cnv_vcf.is_open()) {
 	while (std::getline(cnv_vcf, line)) {
		totalCNV++;
	}
   }
   cnv_vcf.close();

   int countLines = 0;
   cnv_vcf.open (outCNV);
   if (cnv_vcf.is_open()) {
 	while (std::getline(cnv_vcf, line)) {

		//countLines++;
		//float p = (countLines / (float) totalCNV) * (float) 100;
		//printProgBar(p);

		std::vector<string> tmp = split(line, '\t');

		varCall SNV(bamFile, reference, 50, 13, 2, 12, 0.8);

		int resLOH = 3;
		if (tmp[4] == "DEL") {
			//resLOH = SNV.isLOH(tmp[0], atoi(tmp[1].c_str()), atoi(tmp[2].c_str()), tmp[4]);
		}

		vcf_t call;
		call.chr       = tmp[0];
		call.start     = atoi(tmp[1].c_str());
		call.end       = atoi(tmp[2].c_str());
		call.precision = tmp[3] ;
		call.svtype    = tmp[4];
		call.mapq      = atoi(tmp[5].c_str());
		call.breakReads= 0;
		call.assembled = 0;
		call.RDratio   = tmp[6];
		call.RDmad     = tmp[7];
		call.hasRDsupport = true;

		if (resLOH == 3) {
			call.LOHsupport = '.';
		}
		else if (resLOH == 2) {
			call.LOHsupport = "yes";
		}
		else if (resLOH == 1) {
			call.LOHsupport = "no";
		}

		double meanPvalue    = 0.00;
		double pvalue_upstream;
		double pvalue_downstream;
		double pvalue_twosided;
		int counts_5prime, counts_inner, counts_3prime;
		call.covPvalue = meanPvalue;
		call.discordants    = 0;
		call.alleleBalance  = 0;
		call.discPvalue     = 0.00;
		call.cumulativeSize = 0.00;
		vcf_v.push_back(call);
	}
   }
   std::cout << "\n"; 
   return vcf_v;
}


void Cnv::normalize_gc() {

   std::fstream raw;
   raw.open (raw_counts);
   std::map<float, std::vector<int>> mapGC {
	{10.0, vector<int>() },
	{15.0, vector<int>() },
	{20.0, vector<int>() },
	{25.0, vector<int>() },
	{30.0, vector<int>() },
	{35.0, vector<int>() },
	{40.0, vector<int>() },
	{45.0, vector<int>() },
	{50.0, vector<int>() },
	{55.0, vector<int>() },
	{60.0, vector<int>() },
	{65.0, vector<int>() },
	{70.0, vector<int>() },
	{75.0, vector<int>() },
	{80.0, vector<int>() },
	{85.0, vector<int>() },
	{90.0, vector<int>() },
	{95.0, vector<int>() },
	{100.0,vector<int>() },
   };

   // Reading coverage file filling a map with GC content as primary key and an array of counts as a secondary key
   std::map<float, double> mapMedianGC;
   std::vector<int> count_v;
   std::string line;
   if (raw.is_open()) {
         while ( getline(raw, line) ) {
		vector<string> tmp = split (line, '\t');
		int counts = atoi(tmp[3].c_str());
		float gc   = stof(tmp[4].c_str());
		int map    = atoi(tmp[5].c_str());
		for (auto const& x : mapGC) {
			if ( gc <= x.first && gc > x.first-5 ) {
				mapGC[x.first].push_back(counts);
			}
		}
	 }
   }
   raw.close();

   // Calculating median counts per bin
   for (auto & x: mapGC) {
           std::sort(x.second.begin(), x.second.end());
	   int size = x.second.size();

	   if (size == 0) { 
		continue;
	   }
	   int median;
	   if (size  % 2 == 0)
	   {
	      median = (x.second[size / 2 - 1] + x.second[size / 2]) / 2;
	   }
	   else 
	   {
	      median = x.second[size / 2];
	   }
	   mapMedianGC.insert(std::make_pair(x.first , median));
   }

   std::ofstream normalized_gc_file;
   raw.open (raw_counts);
   normalized_gc_file.open(outDir + "/normalized_gc.bed");

   if ( raw.is_open()) {
         while ( getline(raw, line) ) {

		vector<string> tmp = split (line, '\t');
		int counts = atoi(tmp[3].c_str());
		float gc = stof(tmp[4].c_str());
		int map  = atoi(tmp[5].c_str());

		double medianGC;
		for (auto const& x : mapMedianGC) {
			if ( gc <= x.first && gc > x.first -5 ) {
				medianGC = x.second;
			}
		}
		double normalized_gc;
		double count_ratio;
		if ( !isChrSomatic ( tmp[0] )) {
			normalized_gc =  medianGC ? (double) counts * (medianGerminalCounts/medianGC) : counts;
			//count_ratio = normalized_gc/medianGerminalCounts;
		}
		else {
			normalized_gc = medianGC ? (double) counts * (medianSomaticCounts/medianGC) : counts;
			//count_ratio = normalized_gc/medianSomaticCounts;
		}

		normalized_gc_file << tmp[0] << "\t" << tmp[1] << "\t" << tmp[2] << "\t" << gc << "\t" << map << "\t" << counts << "\t" << normalized_gc << "\n";
	 }
   }
   raw.close();
   normalized_gc_file.close();
}


void Cnv::normalize_mappability() {


   std::fstream normalized_gc_file;
   normalized_gc_file.open(outDir + "/normalized_gc.bed");

   std::map<float, std::vector<int>> mapMap {
	{0.0, vector<int>() },
	{5.0, vector<int>() },
	{10.0, vector<int>() },
	{15.0, vector<int>() },
	{20.0, vector<int>() },
	{25.0, vector<int>() },
	{30.0, vector<int>() },
	{35.0, vector<int>() },
	{40.0, vector<int>() },
	{45.0, vector<int>() },
	{50.0, vector<int>() },
	{55.0, vector<int>() },
	{60.0, vector<int>() },
   };

   // Reading coverage file filling a map with GC content as primary key and an array of counts as a secondary key
   std::map<float, double> mapMedianMap;
   std::vector<int> count_v;
   std::string line;
   if (normalized_gc_file.is_open()) {
         while ( getline(normalized_gc_file, line) ) {

		//chr12	1309677	1309777	42.5743	60	38	32.8182

		vector<string> tmp = split (line, '\t');
		float gc   = stof(tmp[3].c_str());
		int map    = atoi(tmp[4].c_str());
		int counts = atoi(tmp[6].c_str());

		for (auto const& x : mapMap) {
			if ( map <= x.first && map > x.first-5 ) {
				mapMap[x.first].push_back(counts);
			}
		}
	 }
   }
   normalized_gc_file.close();

   // Calculating median counts per bin
   for (auto & x: mapMap) {

           std::sort(x.second.begin(), x.second.end());
	   int size = x.second.size();

	   if (size == 0) { 
		continue;
	   }
	   int median;
	   if (size  % 2 == 0)
	   {
	      median = (x.second[size / 2 - 1] + x.second[size / 2]) / 2;
	   }
	   else 
	   {
	      median = x.second[size / 2];
	   }
	   mapMedianMap.insert(std::make_pair(x.first , median));

   }

   std::ofstream copy_ratios;
   copy_ratios.open(outDir + "/copy_ratios.bed");
   normalized_gc_file.open(outDir + "/normalized_gc.bed");

   if (normalized_gc_file.is_open()) {
         while ( getline(normalized_gc_file, line) ) {

		vector<string> tmp = split (line, '\t');
		float gc = stof(tmp[3].c_str());
		int map  = atoi(tmp[4].c_str());
		int counts = atoi(tmp[6].c_str());

		double medianMap;
		for (auto const& x : mapMedianMap) {
			if ( map <= x.first && map > x.first -5 ) {
				medianMap = x.second;
			}
		}
		double normalized_map;
		double count_ratio;
		if ( !isChrSomatic ( tmp[0] )) {
			normalized_map = medianMap ? (double) counts * (medianGerminalCounts/medianMap) : counts;
			count_ratio = normalized_map/medianGerminalCounts;
		}
		else {
			normalized_map = medianMap ?  (double) counts * (medianSomaticCounts/medianMap) : counts;
			count_ratio = normalized_map/medianSomaticCounts;
		}
		copy_ratios << tmp[0] << "\t" << tmp[1] << "\t" << tmp[2] << "\t" << gc << "\t" << map << "\t" << normalized_map << "\t" << count_ratio << "\n";
	 }
   }
   normalized_gc_file.close();


}










