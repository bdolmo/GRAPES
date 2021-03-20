
#ifndef UTILS_CPP_
#define UTILS_CPP_

#include <iostream>
#include <vector>
#include <string> 
#include <map>
#include "math.h"
#include <iomanip>
#include <algorithm>
#include <cmath>
#include <numeric>
#include <stdint.h>
#include <stdbool.h>
#include <assert.h>

#include "kfunc.h"
#include <boost/algorithm/string/replace.hpp>
#include <boost/algorithm/string.hpp>

#include "SeqLib/RefGenome.h"
#include "SeqLib/BamReader.h"
#include "SeqLib/BWAWrapper.h"
#include "SeqLib/BamRecord.h"
#include "SeqLib/BFC.h"



#define TWO_BIT_MASK (3)
#define BITS_PER_BYTE (8)
#define BIG_ENOUGH (1024)

using namespace std;
string IntToString ( int& );
using namespace SeqLib;

//#########################################
int inline mostFrequentPosition ( vector<int> positions ) {

	map<int, int> m;
	int maxCount = 0;
	int currentMax;
	int mostCommon = 0;
	for(int i=0;i < positions.size(); i++)
	{
	    int updated = m[positions[i]]++;  //Increment the value of key for counting occurances        
	    updated++; // due to post increment 
	    if (maxCount < updated) {
		 maxCount = updated;
		 currentMax = i;
		 mostCommon = positions[i];
	    }
	}
	return mostCommon;
}

int inline whichReadLength ( std::string & bamFile ) {

	// opening bam file	
	BamReader reader;
	reader.Open(bamFile);
	BamRecord r;
	int counter = 0;
	vector<int> array;
	while (reader.GetNextRecord(r)) {
		counter++;
		array.push_back(r.Length());
		if (counter > 100) {
			break;
		}
	}
	int readLength = mostFrequentPosition(array);
	return readLength;
}


inline float computeGC ( std::string chr, int start, int end, std::string& reference ) {
	
	RefGenome ref;
	ref.LoadIndex(reference);

	std::string sequence = ref.QueryRegion(chr, start, end);
	int count = 0;
	for (auto& ntd : sequence) {
		if (ntd == 'C' || ntd == 'c' || ntd == 'G' || ntd == 'g' ) {
			count++;
		}
	}
	float gc_content = (count/ (float) sequence.length()) * (float) 100;
	return gc_content;
}

//######################################### Calculate median
int inline calculateMedian ( vector<int>& Array ) {

   int size = Array.size();
   int median;
   sort(Array.begin(), Array.end());

   if (size  % 2 == 0)
   {
      median = (Array[size / 2 - 1] + Array[size / 2]) / 2;
   }
   else 
   {
      median = Array[size / 2];
   }

   return median;
}

//######################################### return 1 is chr is somatic; 0 if no
int inline isChrSomatic ( string chrom ) {

	if (chrom == "23") { 
		return 0;
	}
	if (chrom == "chr23") {
		return 0;
	}
	if (chrom == "chrX") {
		return 0;
	}
	if (chrom == "X") {
		return 0;
	}
	return 1;
}

//######################################### return 1 is chr is somatic; 0 if no
int inline isChrX ( string chrom ) {

	if (chrom == "23") { 
		return 1;
	}
	if (chrom == "chr23") {
		return 1;
	}
	if (chrom == "chrX") {
		return 1;
	}
	if (chrom == "X") {
		return 1;
	}
	return 0;
}

//######################################### return chromosme lexicographical format
string inline returnChromLexicoGraphical ( string chrom ) {

	if (chrom == "23") { 
		chrom = "chrX";
	}
	if (chrom == "24") { 
		chrom = "chrY";
	}
	if (chrom == "chr23") {
		chrom = "chrX";
	}
	if (chrom == "chr24") {
		chrom = "chrY";
	}
	if (chrom == "chrM") {
		chrom = "chrM";
	}
	if (chrom == "M") {
		chrom = "chrM";
	}
	if (chrom == "MT") {
		chrom = "chrM";
	}
	return chrom;
}

//#########################################
double inline getCoverage( std::string& bamFile, string chr, int start, int end) {

	int st_A = start-1;
	int st_B = start;
	int ed_A = end-1;
	int ed_B = end;
	std::string upstream_area, downstream_area;

	upstream_area   = chr + ":" + IntToString(st_A) + "-" + IntToString(st_B);
	downstream_area = chr + ":" + IntToString(ed_A) + "-" + IntToString(ed_B);

	BamReader br;
	br.Open(bamFile);
	BamRecord r;

	std::vector<std::string> regions;
	regions.push_back(upstream_area);
	regions.push_back(downstream_area);
	int count = 0;
	for (auto& i : regions ) {
		count = 0;
		GenomicRegion gr(i, br.Header());
		br.SetRegion(gr);
		while (br.GetNextRecord(r)) {
	      		count++;
			//if (count > 200) { break; }
		}
	}
	double cov;
	if (count == 0) {
		cov = 0.00;
	}
	else {	
		cov = count/2;
	}
	return cov;		
}
//#########################################
// Split function (works faster than boost implementation)
vector<string> inline split( std::string const& original, char separator )
{
    std::vector<std::string> results;
    std::string::const_iterator start = original.begin();
    std::string::const_iterator end = original.end();
    string::const_iterator next = std::find( start, end, separator );
    while ( next != end ) {
        results.push_back( string( start, next ) );
        start = next + 1;
        next = std::find( start, end, separator );
    }
    results.push_back( std::string( start, next ) );
    return results;
}

//#########################################
vector<string> inline splitSoftClip( string read ) {

 vector <string> tmp;
 tmp = split (read, '\t');
 string cigar    = tmp[7];
 string sequence = tmp[11];

 int softclipped, matched;
 string clipped_seq, matched_seq;
 vector<string> s;
 vector<string> m;
 if (tmp[0] == "LEFT") {
	s = split(cigar, 'S');
	m = split(s[1], 'M');

	softclipped = atoi(s[0].c_str());
	matched     = atoi(m[0].c_str());

	clipped_seq = sequence.substr(0,softclipped);
	matched_seq = sequence.substr(softclipped);
 }
 if (tmp[0] == "RIGHT") {
	m = split(cigar, 'M');
	s = split(m[1], 'S');

	matched     = atoi(m[0].c_str());
	softclipped = atoi(s[0].c_str());

	clipped_seq = sequence.substr(matched);
	matched_seq = sequence.substr(0, matched);
 }

 vector<string> seqs;//this vector will be returned 

 seqs.push_back(matched_seq);
 seqs.push_back(clipped_seq);

 return seqs;
}


//#########################################
string inline revComp(string seq) {

	string rev_seq(seq);
	reverse(rev_seq.begin(),rev_seq.end());
	for (int i = 0; i < rev_seq.length(); i++) {
		switch (rev_seq[i]) 
		{
			case 'A': rev_seq[i] = 'T'; break;
			case 'T': rev_seq[i] = 'A'; break;
			case 'C': rev_seq[i] = 'G'; break;
			case 'G': rev_seq[i] = 'C'; break;
			case 'a': rev_seq[i] = 't'; break;
			case 't': rev_seq[i] = 'a'; break;
			case 'c': rev_seq[i] = 'g'; break;
			case 'g': rev_seq[i] = 'c'; break;
		}
	}
	return rev_seq; 
}



//#########################################
 float inline poisson_pmf(int k, double lambda) {
    // https://en.wikipedia.org/wiki/Poisson_distribution#Definition
     return pow(M_E, k * log(lambda) - lambda - lgamma(k + 1.0));
 }

//#########################################

double inline shannonEntropy ( std::string teststring ) {

   std::map<char , int> frequencies ;
   for ( char c : teststring )
     frequencies[ c ] ++ ;
   int numlen = teststring.length( ) ;
   double infocontent = 0 ;
   for ( std::pair<char , int> p : frequencies ) {
      double freq = static_cast<double>( p.second ) / numlen ;
      infocontent += freq * log2( freq ) ;
   }
   infocontent *= -1 ;

   return infocontent;
}
//#########################################

double inline computeDiscordantClusterSignificance ( long int& cumulativeSizes, long& genomeSize, int& numInserts, int& nDiscordants  ) {
	// Based on BreakDancer's confidence interval

	long int csni = cumulativeSizes*numInserts;

	if (csni < 0) {
		return 10e-20;
	}
	double lambda = static_cast<double>(cumulativeSizes*numInserts)/genomeSize;

	double pvalue   = poisson_pmf(nDiscordants, lambda);

	return pvalue;
}

//#########################################
 std::tuple<double, double, double, int, int, int, int, int, int> inline calculatePvalueCoverage( std::string chr, int start , int end, std::string svtype, std::string precision, std::string bamFile)  {

	// This function calculates the significance (poisson dist) between inner-breakpoint counts and outter-breakpoint counts
	//         <- outer       inner ->   <- inner      outer ->
	//eg.    #############|_________________________|##############
	//       ######### break1 ################### break2 ##########

	BamReader br;
	br.Open(bamFile);
	BamRecord r;

	if (svtype == "IMPRECISE") {
		if (start+100 < end-100) {
			start = start+100;
			end = end+100;
		}
	}
	int size = end-start;

	if (size > 1000) {
		size = 1000;
	}

	int start_upstream = start - 2000;
	int end_upstream   = start;
	std::string upstream_area = chr + ":" + std::to_string(start_upstream) + "-" + std::to_string(end_upstream);

	int inner_start = start;
	int inner_end   = start + size;
	std::string inner_area = chr + ":" + std::to_string(inner_start) + "-" + std::to_string(inner_end);

	int start_downstream = end;
	int end_downstream   = end + 2000;
	std::string downstream_area = chr + ":" + std::to_string(start_downstream) + "-" + std::to_string(end_downstream);

	//cout << "Upstream " << upstream_area << "\t" << inner_area << "\t" << downstream_area << endl;

	std::vector<std::string> regions;
	regions.push_back(upstream_area);
	regions.push_back(inner_area);
	regions.push_back(downstream_area);

	std::vector<int> v_counts;
	std::vector<double> v_coverage;

	std::vector<int> v_cov5prime;
	std::vector<int> v_covInner;
	std::vector<int> v_cov3prime;

	int readLength = whichReadLength(bamFile);
	double coverage = 0.00;

	int c = 0;
	for (auto& i : regions ) {
	
		GenomicRegion gr(i, br.Header());
		br.SetRegion(gr);
		int count = 0;
		int depth = 0;
		std::vector<int> Positions;

		while (br.GetNextRecord(r)) {
			int pos    = r.Position();
			int ori = r.PairOrientation();
			int proper = r.ProperPair();
			if (r.Interchromosomal()) {
				continue;
			}
			if (r.CountBWAChimericAlignments() > 0) { 
				//continue; 
			}
			if (r.CountBWASecondaryAlignments() > 0) { 
				//continue; 
			}
			if (c==1) {
				if (r.NumClip() >= 10) {
					continue;
				}
				if (ori != 0) {
					continue;
				}
				if (!proper) {
					continue;
				}
			}
			
			if (r.MapQuality() < 5) {
				continue;
			}			
			Positions.push_back(pos);
		    count++;
		}

		int st;
		int ed;
		if (c == 0) {
			st = start_upstream;
			ed = end_upstream;
		}
		if (c == 1) {
			st = inner_start;
			ed = inner_end;
		}
		if (c == 2) {
			st = start_downstream;
			ed = end_downstream;
		}

		//cout << c << "\t" <<  st << "  => " << ed << endl;
		for (int j = st; j<=ed; j++) {
			for (auto& p : Positions) {
					if (c == 0) {
					//	cout << p << "\t" << p+readLength << "\t" << j << endl;
					}
				if (p <= j && p+readLength >= j) {	
					depth++;
		
					if (c == 0) {
					//	cout << p << "\t" << depth << endl;
					}
				}
				if (p > j) {
				//	cout << p << " > " << j << "\t" << st << " " << ed << endl;
					break;
				}
			}
			//cout << j <<  "\t" << depth << endl;
			if (c == 0) {
				//cout << depth << endl;
				v_cov5prime.push_back(depth);
			}
			if (c == 1) {

				v_covInner.push_back(depth);
			}
			if (c == 2) {
				v_cov3prime.push_back(depth);
			}
			depth = 0;
		}
		v_counts.push_back(count);
		c++;
	}

				
	int medianCov5prime = calculateMedian(v_cov5prime);
	int medianCovInner  = calculateMedian(v_covInner);
	int medianCov3prime = calculateMedian(v_cov3prime);
	//cout << medianCov5prime << "\t" << medianCovInner << "\t" << medianCov3prime << endl;


   	double fisher_left_p, fisher_right_p, fisher_twosided_p;
   	kt_fisher_exact(medianCov5prime, medianCovInner, medianCovInner, medianCov3prime, &fisher_left_p, &fisher_right_p, &fisher_twosided_p);
	return std::make_tuple(fisher_left_p, fisher_right_p, fisher_twosided_p, medianCov5prime, medianCovInner, medianCov3prime, medianCov5prime, medianCovInner, medianCov3prime);
 }




#endif /* UTILS_CPP_ */
