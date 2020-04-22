#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <algorithm>
#include <tuple>
#include <sstream>
#include <regex>
#include <fstream>

#include "SV.h"
#include "sam_t.h"
#include "vcf_t.h"
//#include "reAlignSoftClip.h"
#include "SeqLib/RefGenome.h"
#include "SeqLib/BamReader.h"
#include "SeqLib/BWAWrapper.h"
#include "SeqLib/BamRecord.h"
#include "SeqLib/GenomicRegion.h"
#include "Utils.cpp"

using namespace SeqLib;
inline bool comp_sam(const sam_t& lhs, const sam_t& rhs){
  return std::tie(lhs.align_pos) < std::tie(rhs.align_pos);
}

 float poisson_pmf(int k, double lambda) {
    // https://en.wikipedia.org/wiki/Poisson_distribution#Definition
     return pow(M_E, k * log(lambda) - lambda - lgamma(k + 1.0));
 }

 bool checkIndels (std::string chr, int pos, int start, std::string cigar, int mapq, int n_reads, int n_assembled, double kmer_diversity, std::vector<vcf_t>& vcf_v, std::string& bamFile) {

	// checking if cigar contains MIM or MDM
	std::regex insertion ("^[0-9]+M[0-9]+I[0-9]+M$");
	std::regex deletion ("^[0-9]+M[0-9]+D[0-9]+M$");

        std::smatch m;
	int found_ins = 0;
	int found_del = 0;
	found_ins = std::regex_search (cigar,m, insertion);
	found_del = std::regex_search (cigar,m, deletion);

	if (found_del != 0 || found_ins != 0) {
 		std::string svtype;
		std::vector<std::string> tmp = split( cigar, 'M' );
		std::vector<std::string> tmp2;
		if (found_ins != 0 ) {
			svtype = "INS";
			tmp2 = split( tmp[1], 'I' );
		}
		if (found_del != 0) {
			svtype = "DEL";
			tmp2 = split( tmp[1], 'D' );
		}
		int match1 = atoi(tmp[0].c_str());
		int indel  = atoi(tmp2[0].c_str());
		int match2 = atoi(tmp[2].c_str());
		int st     = start + pos -500 + match1;
		int ed     = start + pos -500 + match1 +indel;
		int size = ed - st;

		vcf_t call;
		call.chr =  chr;
		call.start = st;
		call.end = ed;
		call.precision = "PRECISE";
		call.svtype = svtype;
		call.breakReads = n_reads;
		call.assembled = n_assembled;
		call.discordants = 0;
		call.covPvalue = 0.00;
		call.alleleBalance = 0;

		call.cumulativeSize = 0;
		call.discPvalue = 0.00;
		call.mapq = mapq;
		call.kdiv = kmer_diversity;
		call.LOHsupport = ".";
		call.RDratio = ".";
		call.RDmad   = ".";
		call.hasRDsupport = false;
		vcf_v.push_back(call);

		return 1;
	}
	else {
		return 0;
	}
 }

 bool SV::classify(std::vector<vcf_t>& vcf_v, std::string& suspectedSVtype) {

    int addNtdPositions;

	if (svlength == "small") {
		addNtdPositions = 500;
	}
	else if (svlength == "medium") {
		addNtdPositions = 0;
	}
	else if (svlength == "large") {
		addNtdPositions = 1000;
	}
	std::regex soft_regex_SMS ("^[0-9]+S[0-9]+M[0-9]+S$");
	std::regex soft_regex_SIMS("^[0-9]+S[0-9]+M[0-9]+I[0-9]+S$");
	std::regex soft_regex_L ("^[0-9]+S[0-9]+M$");
	std::regex soft_regex_R ("^[0-9]+M[0-9]+S$");

	int mapq_average;

	//ofstream printOut;

	// Skipping un-assembled reads 
	if ( total_assembled == 0 && total_reads > 1) {
		return 0;
	}

	// if we have a contig that has not created an split read throught BWA alignment
	if (reads.size() == 1) {

		mapq_average = reads[0].mapq;
		// Check if CIGAR contains MIM or MDM indels
		if ( checkIndels(chrA, reads[0].pos, posA, reads[0].cigar, mapq_average, total_reads, total_assembled, kmer_diversity, vcf_v, bamFile) == 1 ) {

		}
		// Perform a separated split-read aligment. This can recover true SV when BWA fails split-read alignment
        else {
			std::smatch m;
			if (std::regex_search (reads[0].cigar, m, soft_regex_L)) {
				softclip_type = "LEFT";
			}
			else if (std::regex_search (reads[0].cigar, m, soft_regex_R)) {
				softclip_type = "RIGHT";
			}
		}
	}

 	// map by read name (map Of Name)
	std::map<std::string, std::string> mapName;
	std::multimap<std::string, sam_t> multimapName;

	std::vector<int> mapq_vector;

	for (int i = 0; i < reads.size(); i++) {
		mapName.insert(std::pair<std::string, std::string>(reads[i].read_name, reads[i].read_name));
		multimapName.insert(std::pair<std::string, sam_t>(reads[i].read_name, reads[i]));
	}
	std::vector<sam_t> ordered_hits;
	std::map<std::string, std::string> strands;
	for (std::map<std::string,std::string>::iterator it=mapName.begin(); it!=mapName.end(); ++it) {

		std::pair <std::multimap<std::string, sam_t>::iterator, std::multimap<std::string, sam_t>::iterator> ret;

  		ret = multimapName.equal_range(it->first);
		std::vector<sam_t> tmp_hits;
		for (std::multimap<std::string, sam_t>::iterator ret=multimapName.begin(); ret!=multimapName.end(); ++ret) {
			mapq_vector.push_back(ret->second.mapq);
			if (reads.size() == 1) {
				strands.insert(std::pair<std::string, std::string>(ret->second.strand, ret->second.strand));
				tmp_hits.insert(tmp_hits.begin(),ret->second);
				break;
			}			
			if (ret->second.read_name != it->first) {
				continue;
			}

			// skipping aligments with mapq score less than 10
			if (ret->second.mapq < 10) {
				continue;
			}

			tmp_hits.push_back(ret->second);
			strands.insert(std::pair<std::string, std::string>(ret->second.strand, ret->second.strand));
		}
		std::sort(tmp_hits.begin(), tmp_hits.end(), comp_sam);
		ordered_hits.insert( ordered_hits.end(), tmp_hits.begin(), tmp_hits.end() );
	}

		mapq_average = std::accumulate(mapq_vector.begin(), mapq_vector.end(), 0.0) / mapq_vector.size();
		std::string first_name;
		int first_align, first_length;
		std::string first_strand;
		std::string first_cigar;
		int first_clipped;
		
		std::string second_name;
		int second_align, second_length;
		std::string second_strand;
		std::string second_cigar;
		int second_clipped;

		int total_hits = 0;
		int succeed = 0;
		int count = 1;		

		if (ordered_hits.size() == 0) {
			return 0;
		}

		for (int i= 0; i < ordered_hits.size()-1; i++) {
			if (ordered_hits.size() == 1) {
				//Indels present in the alignment
				if ( checkIndels(chrA, ordered_hits[0].pos, posA, ordered_hits[0].cigar, mapq_average, total_reads, total_assembled, kmer_diversity, vcf_v, bamFile) == 1 ) {
				}
			}

			first_name    = ordered_hits[i].read_name;
			first_align   = ordered_hits[i].pos;
			first_length  = ordered_hits[i].align_end - ordered_hits[i].align_pos;
			first_strand  = ordered_hits[i].strand;
			first_cigar   = ordered_hits[i].cigar;
			first_clipped = ordered_hits[i].clipped;

			second_name    = ordered_hits[i+1].read_name;
			second_align   = ordered_hits[i+1].pos;
			second_length  = ordered_hits[i+1].align_end - ordered_hits[i+1].align_pos;
			second_strand  = ordered_hits[i+1].strand;
			second_cigar   = ordered_hits[i+1].cigar;
			second_clipped = ordered_hits[i+1].clipped;		

			float pvalue_upstream;
			float pvalue_downstream;

			if (multimapName.count(first_name) < 2) {
				continue;
			}
			if (first_name == second_name) {

				// Deletions and duplications
				if (first_strand == second_strand ) {
					if (second_align > first_align + first_length) {
						int start, end; 
						int tmpStart, tmpEnd;

						if (svlength == "large") {
							tmpStart = first_align + first_length + posA;
							tmpEnd   = second_align + posB - addNtdPositions;
						}
						if (svlength == "small"|| svlength == "medium") {
							tmpStart = first_align + first_length + posA - addNtdPositions +1;
							tmpEnd   = second_align + posA - addNtdPositions;
						}
						if (tmpStart > tmpEnd) {
							end = tmpStart;
							start = tmpEnd;
						}
						else {
							start = tmpStart;
							end   = tmpEnd;
						}
						int length =  end-start;
						double meanPvalue = 0.00;
						int size = end-start;
						vcf_t call;
						call.chr =  chrA;
						call.start = start;
						call.end = end;
						call.precision = "PRECISE";
						call.svtype = "DEL";
						call.breakReads = total_reads;
						call.assembled = total_assembled;
						call.covPvalue = meanPvalue;
						call.discordants = nDiscordants;
						call.alleleBalance = 0;
						call.mapq = mapq_average;
						call.kdiv = kmer_diversity;
						call.RDratio = ".";
						call.RDmad   = ".";
						call.hasRDsupport = false;
						call.LOHsupport = ".";
						vcf_v.push_back(call);
//vcf_out << l.chr << "\t" << l.start << "\t" << l.end << "\t" << l.precision << ";SVTYPE=" << l.svtype << ";MAPQ=" << l.mapq << ";KDIV=" << l.kdiv << ";BREAKREADS=" << l.breakReads << ";ASSEMBLED=" << l.assembled << ";PE=" <<
//	l.discordants << ";CSDISC=" << l.cumulativeSize << ";NINS=" << numInserts << ";RDratio=" << l.RDratio << ";RDmad=" << l.RDmad << ";RDsupp=" << l.hasRDsupport << ";LOHsupp=" << l.LOHsupport << "\n";

						//printOut << call.chr << "\t" << call.start << "\t" << call.end << "\t" << call.precision << ";SVTYPE=" << call.svtype << ";MAPQ=" << call.mapq << ";KDIV=" << call.kdiv;
						//printOut << ";BREAKREADS=" << call.breakReads << ";ASSEMBLED=" << call.assembled << ";PE=" << call.discordants << 
						if ( size >=50 && (suspectedSVtype == "DEL" || suspectedSVtype == "undetermined")) {
							succeed = 1;
						}
					}
					// Duplications
					if (second_align < first_align + first_length) {
						int start, end;
						int tmpStart, tmpEnd;
						if (svlength == "large" ) {
							end = first_align + first_length + posB -addNtdPositions;
							start = second_align + posA;
						}
						if (svlength == "small"|| svlength == "medium") {
							tmpEnd = second_align + posA - addNtdPositions;
							tmpStart = first_align + first_length + posA-addNtdPositions+1;
							if (tmpStart > tmpEnd) {
								end = tmpStart;
								start = tmpEnd;		
							}
							else {
								start = tmpStart;
								end   = tmpEnd;
							}
						}
						int length =  end-start;
						double meanPvalue = 0.00;
						int size = end-start;
						vcf_t call;
						call.chr =  chrA;
						call.start = start;
						call.end = end;
						call.precision = "PRECISE";
						call.svtype = "DUP";
						call.breakReads = total_reads;
						call.assembled = total_assembled;
						call.covPvalue = meanPvalue;
						call.discordants = nDiscordants;
						call.alleleBalance = 0;
						call.discPvalue = pvalue_discordant;
						call.mapq = mapq_average;
						call.kdiv = kmer_diversity;
						call.RDratio = ".";
						call.RDmad = ".";
						call.hasRDsupport = false;
						call.LOHsupport = ".";
						vcf_v.push_back(call);
				
						if (size >=50 && (suspectedSVtype == "DUP" || suspectedSVtype == "undetermined")) {
							succeed = 1;
						}
					}
				}
				// Inversions
				else if ( first_strand != second_strand )  {

					int tmpStart = first_align + first_length;
					int tmpEnd   = second_align + second_length;

					int start = tmpStart + posA - 500 +1;
					int end   = tmpEnd + posA - 500;

					if (svlength == "small" ||  svlength == "medium") {

						if (second_align == first_align + first_length) {

							if (first_strand == "-" && second_strand == "+") {
								std::smatch m;
								//SMS && SMS
								if (std::regex_search (first_cigar, m, soft_regex_SMS) && std::regex_search (second_cigar, m, soft_regex_SMS)) {
									//std::cout << "SMS && SMS" << "\n";
									tmpStart =  second_align + posA + second_length - addNtdPositions;
									tmpEnd   =  first_align + posA + first_length - addNtdPositions;
								}
								// MS
								if (std::regex_search (first_cigar, m, soft_regex_R)) {
									tmpStart =  second_align + posA + second_length - addNtdPositions;
									tmpEnd   =  first_align + posA + first_length - addNtdPositions;
								}
								else {
									tmpStart =  second_align + posA - addNtdPositions;
									tmpEnd   =  first_align + posA - addNtdPositions;
								}
							}
							if (first_strand == "+" && second_strand == "-") {
								std::smatch m;
								// SMS
								if (std::regex_search( first_cigar, m, soft_regex_SMS) && !std::regex_search( second_cigar, m, soft_regex_SMS) ) {
									tmpStart = second_align  + posA - addNtdPositions;
									tmpEnd   = first_align + posA - addNtdPositions;
								}
								//SMS && SMS
								else if (std::regex_search (first_cigar, m, soft_regex_SMS) && std::regex_search (second_cigar, m, soft_regex_SMS)) {
									tmpStart =  second_align + second_length + posA - addNtdPositions;
									tmpEnd   =  first_align + first_length + posA - addNtdPositions;
								}
								// MS with SMS
								else if (std::regex_search (first_cigar, m, soft_regex_R) && std::regex_search( second_cigar, m, soft_regex_SMS)) {
									tmpStart =  second_align + posA + second_length - addNtdPositions;
									tmpEnd   =  first_align + posA + first_length - addNtdPositions;
								}
								// MS no SMS
								else if (std::regex_search (first_cigar, m, soft_regex_R) && !std::regex_search( second_cigar, m, soft_regex_SMS)) {	
									tmpStart =  second_align + posA - addNtdPositions;
									tmpEnd   =  first_align + posA - addNtdPositions;
								}
								else {
									tmpStart =  second_align + posA - addNtdPositions;
									tmpEnd   =  first_align + posA - addNtdPositions;
								}
							}
						}

						if (second_align < first_align + first_length) {

							if (first_strand == "-" && second_strand == "+") {
								std::smatch m;
								// SMS
								if (std::regex_search( first_cigar, m, soft_regex_SMS)  && std::regex_search( second_cigar, m, soft_regex_L) ) {
									tmpStart = second_align +  posA - addNtdPositions;
									tmpEnd   = first_align + posA - addNtdPositions;
								}
								if (std::regex_search( first_cigar, m, soft_regex_SMS)  && !std::regex_search( second_cigar, m, soft_regex_SMS) && !std::regex_search( second_cigar, m, soft_regex_L) ) {
									tmpStart = second_align +  posA - addNtdPositions;
									tmpEnd   = first_align +  posA - addNtdPositions;
								}
								//SMS && SMS
								else if ( (std::regex_search( first_cigar, m, soft_regex_SMS) || std::regex_search( first_cigar, m, soft_regex_SIMS)) && (std::regex_search( second_cigar, m, 									soft_regex_SMS) || std::regex_search( second_cigar, m, soft_regex_SIMS))) {
									tmpStart = second_align + second_length +  posA - addNtdPositions;
									tmpEnd   = first_align + first_length + posA - addNtdPositions;
								}
								// SM
								else if (std::regex_search (first_cigar, m, soft_regex_L) &&  !std::regex_search( second_cigar, m, soft_regex_SMS) ) {
									tmpStart = second_align  + posA - addNtdPositions;
									tmpEnd   = first_align  - first_length + posA - addNtdPositions;
								}
								// MS SMS
								else if (std::regex_search (first_cigar, m, soft_regex_R) && std::regex_search( second_cigar, m, soft_regex_SMS)) {
									tmpStart =  second_align + posA + second_length - addNtdPositions;
									tmpEnd   =  first_align + posA + first_length - addNtdPositions;
								}
								// MS	
								else if (std::regex_search (first_cigar, m, soft_regex_R) && !std::regex_search( second_cigar, m, soft_regex_SMS)) {
									tmpStart = second_align + second_length + posA - addNtdPositions;
									tmpEnd   = first_align + first_length + posA - addNtdPositions;
								}
								else if (!std::regex_search(first_cigar, m, soft_regex_R) && !std::regex_search(first_cigar, m, soft_regex_L) && !std::regex_search(first_cigar, m, soft_regex_SMS) 									&& !std::regex_search(first_cigar, m, soft_regex_SIMS)) {
									continue;
								}
								else {
									tmpStart = second_align +  posA- addNtdPositions;
									tmpEnd   = first_align  +  posA - addNtdPositions;
								}
							}
							if ( first_strand == "+" && second_strand == "-") {
								std::smatch m;
								// SMS
								if (std::regex_search( first_cigar, m, soft_regex_SMS) && !std::regex_search( second_cigar, m, soft_regex_SMS) ) {
									tmpStart = second_align +  posA - addNtdPositions;
									tmpEnd   = first_align +  posA - addNtdPositions;
								}
								// SMS && SMS
								else if (std::regex_search( first_cigar, m, soft_regex_SMS) && std::regex_search( second_cigar, m, soft_regex_SMS)) {

									tmpStart = second_align + second_length +  posA - addNtdPositions;
									tmpEnd   = first_align + first_length + posA - addNtdPositions;
								}
								// SM
								else if (std::regex_search (first_cigar, m, soft_regex_L) && !std::regex_search( second_cigar, m, soft_regex_SMS) ) {
									tmpStart = second_align + posA - addNtdPositions;
									tmpEnd   = first_align  + posA - addNtdPositions;
								}
								// SM
								else if (std::regex_search (first_cigar, m, soft_regex_L) && std::regex_search( second_cigar, m, soft_regex_SMS) ) {
									tmpStart = second_align + posA - addNtdPositions;
									tmpEnd   = first_align  + posA - addNtdPositions;
								}
								// MS SMS
								else if (std::regex_search (first_cigar, m, soft_regex_R) && std::regex_search( second_cigar, m, soft_regex_SMS)) {
									tmpStart =  second_align + posA + second_length - addNtdPositions;
									tmpEnd   =  first_align  + posA + first_length  - addNtdPositions;
								}
								else if (!std::regex_search(first_cigar, m, soft_regex_R) && !std::regex_search(first_cigar, m, soft_regex_L) && !std::regex_search(first_cigar, m, soft_regex_SMS)	&& !std::regex_search(first_cigar, m, soft_regex_SIMS)) {
									continue;
								}
								// MS
								else {
									tmpStart = second_align + second_length + posA- addNtdPositions;
									tmpEnd   = first_align  + first_length + posA - addNtdPositions;
								}
							}
						}
					        if (second_align > first_align + first_length) {

							//cout << "tipus 3" << endl;

							if (first_strand == "-" && second_strand == "+") {
								std::smatch m;
								if (std::regex_search( first_cigar, m, soft_regex_SMS)) {
									//std::cout << "SMS" << "\n";
									tmpStart = first_align  +  posA - addNtdPositions;
									tmpEnd   = second_align +  posA - addNtdPositions;
								}
								else if (std::regex_search (first_cigar, m, soft_regex_L)) {
									//std::cout << "SM" << "\n";
									tmpStart = first_align  + posA - addNtdPositions;
									tmpEnd   = second_align + posA - addNtdPositions;
								}
								else {
									//std::cout << "MS" << "\n";
									tmpStart = first_align  + first_length  + posA - addNtdPositions;
									tmpEnd   = second_align + second_length + posA - addNtdPositions;
								}
							}
							// first_strand (+)      second_strand (-)
							if (first_strand == "+" && second_strand == "-") {
								std::smatch m;
								if (std::regex_search( first_cigar, m, soft_regex_SMS) ) {
									//std::cout << "SMS" << "\n";
									std::vector<std::string> M = split( first_cigar, 'M' );
									std::vector<std::string> S = split( M[0], 'S' );
									tmpStart = first_align  + first_length + posA - addNtdPositions;
									tmpEnd   = second_align + first_length + posA - addNtdPositions;
								}
								else if (std::regex_search (first_cigar, m, soft_regex_L)) {
									//std::cout << "SM" << "\n";
									tmpStart = first_align + posA - addNtdPositions;
									tmpEnd   = second_align  + posA - addNtdPositions;
								}
								else  {
									//std::cout << "MS" << "\n";
									tmpStart = first_align  + first_length  + posA - addNtdPositions;
									tmpEnd   = second_align + second_length + posA - addNtdPositions;
								}
							}
						}
					}
					if (svlength == "large") {
						if (second_align == first_align + first_length) {

							if (first_strand == "-" && second_strand == "+") {
								tmpStart =  second_align + posA + second_length;
								tmpEnd   = first_align + posA + first_length;	
							}
							else  {
								tmpStart =  first_align + posA + first_length;
								tmpEnd   = second_align + posA + second_length;	
							}
						}
						if (second_align < first_align + first_length) {

							if (first_strand == "-" && second_strand == "+") {
								tmpStart = second_align + posA;
								tmpEnd   = first_align + posB - addNtdPositions;
							}
							else {
								std::smatch m;
								if (std::regex_search (first_cigar, m, soft_regex_L)) {
									tmpStart = second_align + posA;
									tmpEnd   = first_align  + posB - addNtdPositions;
								}
								else {
									tmpStart = second_align + second_length + posA;
									tmpEnd   = first_align + first_length + posB - addNtdPositions;
								}
							}
							//continue;
						}
					        if (second_align > first_align + first_length ) {

							if (first_strand == "-" && second_strand == "+") {
								tmpStart = first_align + first_length + posA;
								tmpEnd   = second_align + second_length + posB - addNtdPositions;
							}
							// first_strand (+)      second_strand (-)
							else {
								std::smatch m;
								if (std::regex_search (first_cigar, m, soft_regex_L)) {
									tmpStart = first_align + posA;
									tmpEnd   = second_align  + posB - addNtdPositions;
								}
								else  {
									tmpStart = first_align + posA + first_length;
									tmpEnd   = second_align + second_length + posB - addNtdPositions;
								}

							}
						}
					}
					if (tmpStart > tmpEnd) {
						end   = tmpStart;
						start = tmpEnd;
					} else {
						start = tmpStart;
						end = tmpEnd;
					}	
					int length = end - start;
					if (length == 0) {
						continue;
					}
					int size = end-start;

					vcf_t call;
					call.chr =  chrA;
					call.start = start;
					call.end = end;
					call.precision = "PRECISE";
					call.svtype = "INV";
					call.breakReads = total_reads;
					call.assembled = total_assembled;
					call.covPvalue = 0;
					call.discordants = nDiscordants;
					call.alleleBalance = 0;
					call.discPvalue = pvalue_discordant;
					call.mapq = mapq_average;
					call.kdiv = kmer_diversity;
					call.RDratio = ".";
					call.RDmad = ".";
					call.hasRDsupport = false;
					call.LOHsupport = ".";
					vcf_v.push_back(call);
					if (size >=50 && (suspectedSVtype == "INV" || suspectedSVtype == "undetermined")) {
						succeed = 1;
					}
				}
		}
	}
	return succeed;
 }


