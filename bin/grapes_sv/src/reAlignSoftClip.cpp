#include <iostream>
#include <vector>
#include <string> 
#include "Utils.cpp"
#include <sstream>
#include "ssw_cpp.h"
#include "reAlignSoftClip.h"
#include <regex>
#include <fstream>
#include "vcf_t.h"

std::vector<int> returnAlignment (const StripedSmithWaterman::Alignment& alignment);

 std::vector<std::string> splitSoftClip( std::string read, std::string cigar, std::string softclip_type ) {
	std::vector<std::string> s;
	std::vector<std::string> m;
	int soft_n;
	int match_n;

	std::string softclipped_seq;
	std::string matched_seq;
	
	if(softclip_type == "LEFT") {
		s = split(cigar, 'S');
		m = split(s[1], 'M');

		soft_n = atoi(s[0].c_str());
		match_n= atoi(m[0].c_str());
		
		softclipped_seq = read.substr(0, soft_n);
		matched_seq = read.substr(soft_n);

	}
	if(softclip_type == "RIGHT") {
		m = split(cigar, 'M');
		s = split(m[1], 'S');

		soft_n = atoi(s[0].c_str());
		match_n= atoi(m[0].c_str());
		
		matched_seq = read.substr(0, match_n);
		softclipped_seq = read.substr(match_n);
	}
	std::vector<std::string> output;
	output.push_back(matched_seq);
	output.push_back(softclipped_seq);

	return output;
 }

 void reAlignSoftClip::Align(std::vector<vcf_t>& vcf_v) {

 		//vcf_out.open(VCF, std::fstream::app);

		//std::ofstream vcf_out;
		//vcf_out.open(VCF.c_str(), std::ofstream::app);

        	std::smatch m;
		std::regex soft_regex_L ("^[0-9]+M[0-9]+S$");
		std::regex soft_regex_R ("^[0-9]+M[0-9]+S$");

		if ( !std::regex_search (cigar, m, soft_regex_L) && !std::regex_search (cigar, m, soft_regex_R) ) {
			return;
		}

		int addNtdPositions;

		if (svlength == "small") {
			addNtdPositions = 500;
		}
		else if (svlength == "medium") {
			addNtdPositions = posB - posA;
		}
		else if (svlength == "large") {
			addNtdPositions = 1000;
		}

		std::vector<std::string> split_v = splitSoftClip(read, cigar, softclip_type );

		int nreads     = 0;
		int sameStrand = 0;
		int diffStrand = 0;
		int badAlign   = 0;
		int minBadAlign= 0;

		int delScore = 0; 
		int dupScore = 0;
		int invScore = 0;
		int st, ed;
		std::vector <int> STARTS;
		std::vector <int> ENDS;
		std::string svtype;			
		int score_read1, score_read2, mismatches_read1, mismatches_read2;
		double min_score_read1, min_score_read2, max_score_read2;
		char strand_read1, strand_read2;

		int ref_start_read1, ref_end_read1, query_start_read1, query_end_read1;
		int ref_start_read2, ref_end_read2, query_start_read2, query_end_read2;

		int length_read1 = split_v[0].length();
		int length_read2 = split_v[1].length();

		std::string rev_seq = revComp(ref_sequence);

		//read1 is always the matched part
		//read2 is always the clipped part

		for ( int i = 0; i < 2; i++ ) {

			nreads++;
			string seq = split_v[i];	

			// Declares a default Aligner
			StripedSmithWaterman::Aligner ssw_aligner;

			// Declares a default filter
			StripedSmithWaterman::Filter filter;

			// Declares an alignment that stores the result
			StripedSmithWaterman::Alignment alignment;

			// Aligns the query to the ref
			ssw_aligner.Align(seq.c_str(), ref_sequence.c_str(), ref_sequence.length(), filter, &alignment);

			std::vector<int> forward_results;
			forward_results = returnAlignment(alignment);
				
			ssw_aligner.Align(seq.c_str(), rev_seq.c_str(), ref_sequence.length(), filter, &alignment);

			std::vector<int> reverse_results;
			reverse_results = returnAlignment(alignment);

			if (i == 0) {
				min_score_read1 = length_read1*1.7;
			}
			else {
				min_score_read2 = length_read2*1.7;
				max_score_read2 = length_read2*2;
			}
	
			if ( reverse_results[0]-reverse_results[5]*3 > forward_results[0]) {
					if (i == 0) {
						score_read1   = reverse_results[0];
						strand_read1  = '-'; 
						ref_start_read1  = reverse_results[1];
						ref_end_read1    = reverse_results[2];
						query_start_read1= reverse_results[3];
						query_end_read1  = reverse_results[4];
						mismatches_read1 = reverse_results[5];	
					}
					if (i == 1) { 	
						score_read2   = reverse_results[0];
						strand_read2  = '-';
						ref_start_read2  = reverse_results[1];
						ref_end_read2    = reverse_results[2];
						query_start_read2= reverse_results[3];
						query_end_read2  = reverse_results[4]; 
						mismatches_read2 = reverse_results[5];
					}	
			}
			if (reverse_results[0] < forward_results[0]-forward_results[5]*3) {
					if (i == 0) { 	
						score_read1   = forward_results[0];
						strand_read1  = '+';
						ref_start_read1  = forward_results[1];
						ref_end_read1    = forward_results[2];
						query_start_read1= forward_results[3];
						query_end_read1  = forward_results[4]; 
						mismatches_read1 = forward_results[5];
					}
					if (i == 1) { 	
						score_read2   = forward_results[0];
						strand_read2  = '+'; 
						ref_start_read2  = forward_results[1];
						ref_end_read2    = forward_results[2];
						query_start_read2= forward_results[3];
						query_end_read2  = forward_results[4]; 
						mismatches_read2 = forward_results[5];
					}
			}
		}

		// count for inversions
			if ( strand_read1 != strand_read2 ) {
				diffStrand+=2;
				if (softclip_type == "RIGHT") {
					st = ref_end_read1;
					ed   = ref_start_read2;
					STARTS.push_back( st );
					ENDS.push_back( ed );
				}
				if (softclip_type == "LEFT") {
					st = ref_start_read1;
					ed   = ref_end_read2;
					STARTS.push_back( st );
					ENDS.push_back( ed );
				}
			}

			// count if bad quality alignment

			//std::cout << "min_score_r1:" << min_score_read1 << " score_read1:" << score_read1 << " min_score_r2:" << min_score_read2 <<  " score_read2:" << score_read2 << "\n";

			if (score_read1 < min_score_read1 or score_read2 < min_score_read2 ) { badAlign++; }

			if (score_read2 > max_score_read2 ) {

				badAlign++;
			}
			if (score_read2 <  max_score_read2*-1 ) {
				badAlign++;
			}

			// count for others
			if ( strand_read1 == strand_read2 ) {
				sameStrand+=2;
				if (softclip_type == "RIGHT") {
					if (ref_end_read1 < ref_start_read2) {
						delScore++;
						st = ref_end_read1;
						ed = ref_start_read2;
						STARTS.push_back( st );
						ENDS.push_back( ed );
					}	
					if (ref_end_read1 > ref_start_read2) {
						dupScore++;
						st = ref_start_read2;
						ed = ref_end_read1;
						STARTS.push_back( st );
						ENDS.push_back( ed );
					}
				}	
				if (softclip_type == "LEFT") {
					if (ref_start_read1 > ref_end_read2 ) {
						delScore++;
						st   = ref_end_read2;
						ed   = ref_start_read1;
						STARTS.push_back( st );
						ENDS.push_back( ed );
					}	
					if ( ref_end_read2 > ref_start_read1 ) {
						dupScore++;
						st = ref_start_read1;
						ed = ref_end_read2;
						STARTS.push_back( st );
						ENDS.push_back( ed );
					}	
				}
			}

		//std::cout << "delscore: " << delScore << "\t" << "dupScore: " << dupScore << "\n";

		double percSameStrand = 100*(sameStrand/nreads);
		double percDiffStrand = 100*(diffStrand/nreads);
		minBadAlign = (nreads/3)+0.5;

		//std::cout << "badalign: " << badAlign << "\t" << "minbadalign: " << minBadAlign << "\n";
		//std::cout << "perSameStrand: " << percSameStrand << "\t" << "perDiffStrand: " << percDiffStrand << "\n";


		if (badAlign > minBadAlign && badAlign > 0 && minBadAlign >= 0) {
			return;
		}

		int tmpStart = mostFrequentPosition(STARTS);
		int tmpEnd = mostFrequentPosition(ENDS);
		int start_coord, end_coord;

		
		if (svlength == "small") {
		 	start_coord = pos - addNtdPositions + tmpStart +1;
			end_coord = pos - addNtdPositions + tmpEnd +1 ;
		}
		else if (svlength == "medium" || svlength == "large") {
			start_coord = tmpStart + posA +1;
			end_coord   = tmpEnd + posB + 1 - addNtdPositions;

		}
		int len = end_coord - start_coord;

		if (percDiffStrand > 70) {
			svtype = "INV";
			if (tmpEnd > tmpStart) {
				int tstart = tmpStart-tmpEnd+pos;
				end_coord   = start_coord;
				start_coord = tstart;
			}
			if (tmpEnd < tmpStart) {
				int tend = tmpStart - tmpEnd + pos;
				end_coord = tend;
			}

			//vcf_out << chr << "\t" << start_coord << "\t"<< end_coord << "\t" << "PRECISE;SVTYPE=INV;SIZE=" << len << ";READS="<< nreads/2 << endl;
			delScore = 0;
			dupScore = 0;
			invScore = 0;
			//continue;
		}
		if (delScore > dupScore) {
			svtype = "DEL";
			//vcf_out << chr << "\t" << start_coord << "\t"<< end_coord << "\t" << "PRECISE;SVTYPE=DEL;SIZE=" << len << ";READS="<< nreads/2 << endl;
			//std::cout << chr << "\t" << start_coord << "\t"<< end_coord << "\t" << "PRECISE;SVTYPE=DEL;SIZE=" << len << ";READS="<< nreads/2 << endl;
			delScore = 0;
			dupScore = 0;
			invScore = 0;
			//continue;
		}
		if (dupScore > delScore) {
			svtype = "DUP";
			delScore = 0;
			dupScore = 0;
			invScore = 0;
			//continue;
		}
		tmpStart = start_coord;
		tmpEnd = end_coord;

		if (tmpEnd < tmpStart ) {
			end_coord = tmpStart;
			start_coord = tmpEnd;
		}

	

		vcf_t call;
		call.chr =  chr;
		call.start = start_coord;
		call.end = end_coord;
		call.precision = "PRECISE";
		call.svtype = svtype;
		call.breakReads = nreads/2;
		call.assembled = 0;
		call.covPvalue = 0.00;
		call.discordants = 0;
		call.alleleBalance = 0;
		std::cout << "MAQPE\t" << mapq << "\n";
		call.mapq = mapq;
		vcf_v.push_back(call);

}




