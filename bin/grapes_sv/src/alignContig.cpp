#include <iostream>
#include <vector>
#include <string> 
#include "Utils.cpp"
#include <sstream>
#include "ssw_cpp.h"
//#include "reAlignSoftClip.h"
#include <regex>
#include <fstream>
#include "vcf_t.h"
#include "alignContig.h"
#include <cmath>
#include <unordered_map>

// Seqlib headers
#include "SeqLib/RefGenome.h"
#include "SeqLib/BamReader.h"

#include <boost/algorithm/string.hpp>

struct unalignedBlock  {
	int start;
	int end;
	string seq;
	int length = seq.length();

	// Relative positions to initial contig sequence
	int contigStart;
	int contigEnd;
};

struct summaryCigar {
	int numS;
	int numM;
	int numX;
	int numI;
	int numD;
	string collapsed;
};

summaryCigar cigar2Summary ( std::string& cigar) {

	std::string concat;
	//10S30M5S
	int numS = 0;
	int numM = 0;
	int numX = 0;
	int numI = 0;
	int numD = 0;

	//cout << globalVar::

	summaryCigar sCigar;
	for (auto& i : cigar)  {
		
		if (i == 'D' or i == 'I' or i == 'M' or i == 'X' or i == 'S' or i == '=') {
			if (i == 'D') {
				numD++;
			}
			if (i == 'I') {
				numI++;
			}
			if (i == 'M') {
				numM++;
			}
			if (i == 'X') {
				numX++;
			}
			if (i == 'S') {
				numS++;
			}
			concat += i;

		}
		else  {
		}
	}
	sCigar.numS = numS;
	sCigar.numD = numD;
	sCigar.numI = numI;
	sCigar.numM = numM;
	sCigar.numX = numX;
	sCigar.collapsed = concat;
	return sCigar;
}

int checkSeqOccurrences (std::string& sequence, std::string& reference) {

	std::unordered_multimap<string,string> hashTable;
	boost::to_upper(reference);
	std::vector<char> Alphabet {'A', 'C', 'T', 'G'};
	std::vector<std::string> seqFull;
	int seqLength = sequence.length();

	// Here we can generate all possible 1-change kmers
	for (int i = 0; i < sequence.length(); i++ ) {
		for (auto & base : Alphabet) {
			std::string newSeq = sequence;
			newSeq[i] = base;
			seqFull.push_back(newSeq);
		}
	}
	//	cout <<" ok8 " << endl;

	for (int i = 0; i < reference.length()-seqLength; i++) {
		std::string kmer = reference.substr(i, seqLength);
 		hashTable.insert(std::make_pair(kmer,kmer));
	}
	//	cout <<" ok8 out" << endl;

	int numOccurrences = 0;
	numOccurrences += hashTable.count(sequence);
	for (auto& seq : seqFull) {
		if (seq != sequence) {
			numOccurrences += hashTable.count(seq);
		}
	}

	//int numOccurrences = hashTable.count(sequence);
	return numOccurrences;
}

string getAlignedFromCigar ( std::string& cigar, std::string sequence) {

	int flag = 0;
	int ini;
	int end;

	std::string concat;
	//10S30M5S
	for (auto& i : cigar)  {
		if (i == 'D' or i == 'I' or i == 'M' or i == 'X' or i == 'S' or i == '=') {
			int n = std::stoi(concat);
			if (i == 'S') {
				if (flag ==1) {
					end = n;
				}
				if (flag == 0) {
					ini = n;
					flag =1;
				}
			}
			concat = "";				
		}
		else  {
			concat += i;
		}
	}
}

alignmentStruct returnAlignment (const StripedSmithWaterman::Alignment& alignment, int& queryLength, int& refLength, string& refSeq, string& querySeq, string& strand){

	alignmentStruct a;

	a.score      = alignment.sw_score;
	a.refStart   = alignment.ref_begin;
	a.refEnd     = alignment.ref_end;
	
	a.queryStart = alignment.query_begin;
	a.queryEnd   = alignment.query_end;
	a.querySeq   = querySeq;
	int len      = querySeq.length();
	a.queryLength= len;
	a.strand     = strand;
	a.mismatches = alignment.mismatches;
	a.cigar      = alignment.cigar_string;

	a.length     = a.queryEnd-a.queryStart;
	a.refLength  = refLength;
	a.refSeq     = refSeq;
	
	a.numClip;
	a.numMatched;
	a.multimaps;
	a.repeatsOnRead;

 	return a;
}


int treatAligment(alignmentStruct& fwd, alignmentStruct& rev, std::vector<alignmentStruct>& splitReadVec, std::vector<unalignedBlock>& unalignedBlockVec, std::string& contigSeq, std::string& reference, int& readLength) {

	alignmentStruct best;
	if (fwd.score >= rev.score) {
		best = fwd;
	}
	else {
		best = rev;
	}
	string alignedSeq;
	if (best.queryStart < best.querySeq.length()) {
		alignedSeq = best.querySeq.substr(best.queryStart, best.queryEnd-best.queryStart+1);
	}
	int alignedSize = alignedSeq.length();
	float maxMismatches = alignedSeq.length()*0.02;
	int roundMaxMismatches = round(maxMismatches);

	std::vector<std::string> kdivVec;
	kdivVec.push_back(alignedSeq);
	double kdiv = getKmerDiversity(kdivVec);

	//float minAcceptedLength = readLength*0.13;
	float minAcceptedLength = 30.00;
	if (debug) {
		cout << "Kdiv:" << kdiv << "\t" << "Mida:" << alignedSeq.length()<< "\t" << minAcceptedLength << "->" <<  "max_mismatches:" << roundMaxMismatches << "\t" << "mismatches:" << best.mismatches << endl;
	}
	if (best.mismatches <= roundMaxMismatches && alignedSeq.length() > minAcceptedLength && kdiv > 0.05) {
		
		//cout << "Forma adjancencia!\n";

		StripedSmithWaterman::Aligner ssw_aligner;
		StripedSmithWaterman::Filter filter;
		StripedSmithWaterman::Alignment localAlignment;
		ssw_aligner.Align(alignedSeq.c_str(), contigSeq.c_str(), contigSeq.length(), filter, &localAlignment);

		int aSeqLen   = alignedSeq.length();
		int contigLen = contigSeq.length();
		string strand = "+";
		alignmentStruct res = returnAlignment(localAlignment, aSeqLen, contigLen, contigSeq, alignedSeq, strand);
		best.contigStart   = res.refStart;
		best.contigEnd     = res.refEnd;
		best.multimaps     = checkSeqOccurrences(alignedSeq, reference);
		best.repeatsOnRead = checkSeqOccurrences(contigSeq, reference);
		splitReadVec.push_back(best);
	}
	else {
		alignedSize = 0;
	}
	// Aquí ja tenim el millor alineament entre la porció del contig i la referència
	// Ara extreurem les parts del contig que han quedat sense alinear i les introduirem al nou vector
	if (debug) {
		cout << "Segment:" << alignedSeq <<  "mismatches:" << best.mismatches << "\t" << " Query_start: " << best.queryStart << " Query_end: " << best.queryEnd << " Query_length: " << best.queryLength << "\t" << best.cigar << "\t" << "Ref_start: " << best.refStart << " Ref_end: " << best.refEnd << "\tRefLen: " << best.refLength  << endl;
	}
	// If there is an unaligned prefix
	if (best.queryStart > 0) {
		unalignedBlock pblock;
		pblock.start = 0;
		pblock.end   = best.queryStart;


		//cout << best.querySeq.length() << "\t" << best.queryStart << endl;
		pblock.seq   = best.querySeq.substr (0, best.queryStart);
		//cout << "ok2 out" << endl;

		if (pblock.seq.length() > minAcceptedLength) {
			unalignedBlockVec.push_back(pblock);
		}
	}
	// If there is an unaligned suffix
	if (best.queryEnd < best.queryLength-1) {

		unalignedBlock sblock;
		sblock.start = best.queryEnd;
		sblock.end   = best.queryLength-1;

		//cout << best.queryEnd << "\t" << best.queryLength << endl;

		sblock.seq   = best.querySeq.substr (best.queryEnd);
		//cout << "ok3 out" << endl;

		if (sblock.seq.length() > minAcceptedLength) {
			unalignedBlockVec.push_back(sblock);
		}
	}
	return alignedSize;
 }

 //######################################
 struct indel {
	 string type;
	 int start;
	 int end;
	 int size;
 };

//######################################
std::vector<indel> cigar2indel(std::string& cigar, int& pos) {

	std::string concat;
	int gaps = 0;
	int matches = 0;
	int globalPosition = 0;

	std::string globalConcat;
	std::vector<indel> indelVec;

	for (auto& i : cigar)  {
		if (i == 'D' or i == 'I' or i == 'M' or i == 'X' or i == 'S' or i == '=') {
			if (i == 'D' or i == 'I') {
				int n = std::stoi(concat);
				indel INDEL;
				if (i == 'D') {
					INDEL.type = "DEL";
					INDEL.start = pos+globalPosition;
					INDEL.end   = pos+globalPosition+n;
					INDEL.size  = n;
					indelVec.push_back(INDEL);
				}
				if (i == 'I') {
					INDEL.type = "INS";
					INDEL.start = pos+globalPosition;
					INDEL.end   = pos+globalPosition+1;
					INDEL.size  = n;
					indelVec.push_back(INDEL);
				}
			}
			if (i == '=') {
				int n = std::stoi(concat);
				//globalPosition+=std::stoi(concat);
				matches +=n;	
			}
			globalPosition+=std::stoi(concat);
			concat = "";				
		}
		else  {
			concat += i;
		}
	}
	return indelVec;
}

//######################################

 int passesCigarRules(std::string& chr, std::string& cigar, int& pos, int& startCoordinate, std::string& contigSeq, std::string& refSeq, std::string& genome) {

	boost::to_upper(contigSeq);
	boost::to_upper(refSeq);

	std::string concat;
	int gaps = 0;
	int matches = 0;
	int globalPosition = 0;

	int refPosition = 0;
	std::string globalConcat;
	std::vector<indel> indelVec;
	bool passes = false;

	//expanded refseq to check multimapping segments
	RefGenome ref;
	ref.LoadIndex(genome);

	int expandedStart   = startCoordinate-1000;
	int expandedEnd     = startCoordinate+1000;

	// chr3   +    38202755  38202793	 chr3   +    38203395  38203442

	std::string expandedRef = ref.QueryRegion(chr, expandedStart, expandedEnd);

	for (auto& i : cigar)  {
		if (i == 'D' or i == 'I' or i == 'M' or i == 'X' or i == 'S' or i == '=') {

			int n = std::stoi(concat);
			int m = std::stoi(concat);
			if (i =='I') {
				m = 0;
			}
			if (i =='D') {
				n = 0;
			}

			if (i == 'M') {
				int start = pos+refPosition;
				if (n >= 30) {
					std::string overlapRef    = refSeq.substr(start, n);
					std::string overlapContig = contigSeq.substr(globalPosition, n);

					std::vector<std::string> kdivVec;

					kdivVec.push_back(overlapContig);
					double kdiv = getKmerDiversity(kdivVec);

					int multimaps = checkSeqOccurrences(overlapContig, expandedRef);
					int mismatches = 0;

					float maxMismatches = overlapRef.length()*0.02;
					int roundMaxMismatches = round(maxMismatches);

					for (int i = 0; i < overlapRef.length(); i++) {
						if (overlapRef[i] != overlapContig[i]) {
							mismatches++;
						}
					}
					if (debug) {
						cout << "Kdiv" << "\t" << kdiv << endl;
						cout << "Multis" << "\t" << multimaps << endl;
						cout << overlapRef << endl;
						cout << overlapContig << endl;
						cout << roundMaxMismatches << endl;
						cout << mismatches << endl;
					}
					if (mismatches > roundMaxMismatches) {
						return 0;
					}
					else if (multimaps > 1) {
						return 1;
					}
					else if (kdiv < 0.08) {
						return 2;
					}
				}
			}
			if (i == '=') {
				int n = std::stoi(concat);
				//globalPosition+=std::stoi(concat);
				matches +=n;	
			}
			globalPosition+=n;
			refPosition+=m;
			concat = "";				
		}
		else  {
			concat += i;
		}
	}
	return 3;
 } 
 

 //######################################

 std::vector<vcf_t> alignContig::detectSV(std::vector<alignmentStruct>& splitReadVec, int& totalReads, int&totalReadsAssembled, double&kmerDiv, int& sizeV, int& meanMapQ, std::string& contigSeq, std::string& refSeq) {


	std::string debugF = "/home/bernat/Escritorio/HG002_CALLERS/GRAPES/GRAPES_DEBUG.txt";
	std::ofstream debugFile;
	debugFile.open(debugF, std::fstream::app);
	std::vector<vcf_t> vcf_v;

	int sizeSplitReadVec = splitReadVec.size();

	std::regex insertion ("(^[3-9]\\d|\\d{3,})+M[0-9]+I([3-9]\\d|\\d{3,})M([0-9]S)?$");
	std::regex deletion ("(^[3-9]\\d|\\d{3,})+M.*([3-9]\\d|\\d{3,})M([0-9]S)?$");

	std::regex soft_regex_L ("^[0-9]+S[0-9]+M$");
	std::regex soft_regex_R ("^[0-9]+M[0-9]+S$");
	std::regex soft_regex_SMS ("^[0-9]+S[0-9]+M[0-9]+S$");
	std::regex match_regex("^[0-9]+M$");
    std::smatch m;
	if (sizeSplitReadVec == 1) {
		
		alignmentStruct segmentA = splitReadVec[0];
		int segPos = segmentA.refStart+genomicStart;
		int valid_ins = std::regex_search (segmentA.cigar,m, insertion);
		int valid_del = std::regex_search (segmentA.cigar,m, deletion);
		if (debug) {
			cout << valid_del << endl;
		}

		if (!valid_ins && !valid_del) {
			debugFile << chr << "\t" << segPos << "\t" << segPos+1 << "\t" << segmentA.cigar << "\tFilter_Invalid_cigar" << endl;
		}

		std::vector<indel> indelsA = cigar2indel(segmentA.cigar, segPos);
  		int passes = passesCigarRules(chr, segmentA.cigar, segmentA.refStart, segPos, contigSeq, refSeq, genome);

		if (!valid_ins && !valid_del || passes != 3) {

			if (debug) {
				if (passes != 3) {
					cout << "Filtered by many mismatches!" << endl;
				}
			}
			if (passes == 0) {
				debugFile << chr << "\t" << segPos << "\t" << segPos+1 << "\t" << segmentA.cigar << "\tFilter_High_mismatches" << endl;
			}
			if (passes == 1) {
				debugFile << chr << "\t" << segPos << "\t" << segPos+1 << "\t" << segmentA.cigar << "\tFilter_Multimaps" << endl;
			}
			if (passes == 2) {
				debugFile << chr << "\t" << segPos << "\t" << segPos+1 << "\t" << segmentA.cigar << "\tFilter_Kdiv" << endl;
			}
			return vcf_v;
		}

		if (debug) {
			cout << "cigar is okay\n";
		}


		for (int i = 0; i< indelsA.size(); i++) {

			indel ind    = indelsA[i];
			int startPos = ind.start;
			int endPos   = ind.end;
			int size     = ind.size;
			string type  = ind.type;

			if (size >= 20) {
				vcf_t call;
				call.chr =  chr;
				call.start = startPos;
				call.end = endPos;
				call.precision = "PRECISE";
				call.svtype = type;
				call.breakReads = totalReads;
				call.assembled  = totalReadsAssembled;
				call.mapq =  meanMapQ;
				call.hasRDsupport = false;
				call.RDratio = ".";
				call.RDmad   = ".";
				call.LOHsupport = ".";
				call.kdiv = kmerDiv;

				double pvalue_upstream;
				double pvalue_downstream;
				double pvalue_twosided;
				int counts_5prime, counts_inner, counts_3prime;
				double meanPvalue = 0.00;
				call.covPvalue      = meanPvalue;
				call.discordants    = 0;
				call.alleleBalance  = 0;
				call.discPvalue     = 0;
				call.cumulativeSize = 0;
				//call.reciprocalSupport = "no";
				if (debug) {
					cout << "MIDA1 VARCALL: "<<  call.start << "\t" << call.end << "\t" << call.svtype << endl;
				}
				vcf_v.push_back(call);
			}
		}
	}
	
	if (sizeSplitReadVec > 1) {

		for (int i = 0; i < splitReadVec.size()-1; i++) {

			alignmentStruct segmentA = splitReadVec[i];
			alignmentStruct segmentB = splitReadVec[i+1];

			summaryCigar csA = cigar2Summary(segmentA.cigar);
			summaryCigar csB = cigar2Summary(segmentB.cigar);

			//cout << csA.collapsed << "\t" << csB.collapsed << endl;
			int segPos = segmentA.refStart+genomicStart;
			int startPos, endPos, size;
			string type;

			std::vector<indel> indelsA = cigar2indel(segmentA.cigar, segPos);
			int valid_ins = std::regex_search (segmentA.cigar,m, insertion);
			int valid_del = std::regex_search (segmentA.cigar,m, deletion);

			if (valid_ins or valid_del) {
				//for (auto& ind : indelsA) {
				for (int i = 0; i< indelsA.size(); i++) {
					indel ind = indelsA[i];
					
					startPos = ind.start;
					endPos   = ind.end;
					size     = endPos - startPos;
					type     = ind.type;

					if (size >= 20) {
						vcf_t call;
						call.chr =  chr;
						call.start = startPos;
						call.end = endPos;
						call.precision = "PRECISE";
						call.svtype = type;
						call.breakReads = totalReads;
						call.assembled  = totalReadsAssembled;
						call.mapq =  meanMapQ;
						call.hasRDsupport = false;
						call.RDratio = ".";
						call.RDmad   = ".";
						call.LOHsupport = ".";
						call.kdiv = kmerDiv;

						double pvalue_upstream;
						double pvalue_downstream;
						double pvalue_twosided;
						int counts_5prime, counts_inner, counts_3prime;
						double meanPvalue = 0.00;
						call.covPvalue      = meanPvalue;
						call.discordants    = 0;
						call.alleleBalance  = 0;
						call.discPvalue     = 0;
						call.cumulativeSize = 0;
						vcf_v.push_back(call);
					}
				}
			}
			if ( csA.numS <= 1 && csB.numS <= 1) {

				if (segmentA.strand == segmentB.strand) {

					if (segmentB.refStart > segmentA.refEnd) {
						if (segmentA.multimaps > 1 || segmentB.multimaps > 1) {
							debugFile << chr << "\t" << genomicStart+segmentA.refEnd+1 << "\t" << genomicStart+segmentB.refStart << "\t" << segmentA.cigar << "\tFilter_Multimaps" << endl;
						}
						if (segmentA.repeatsOnRead > 1 || segmentB.repeatsOnRead > 1) {
							debugFile << chr << "\t" << genomicStart+segmentA.refEnd+1 << "\t" << genomicStart+segmentB.refStart << "\t" << segmentA.cigar << "\tFilter_Repeats_on_read" << endl;
						}
					}
					// Defining deletions
					if (segmentB.refStart > segmentA.refEnd && segmentA.multimaps <= 1 && segmentA.repeatsOnRead <= 1 && segmentB.multimaps <= 1 && segmentB.repeatsOnRead <=1) {
						//cout << "DINS genomicStart:" << genomicStart << "\t" << segmentA.refStart << "\t" << genomicStart+segmentB.refStart<< endl;

						startPos = genomicStart+segmentA.refEnd+1;
						endPos   = genomicStart+segmentB.refStart;
						size     = endPos - startPos;
						type     = "DEL";
						//cout << size << endl;
					}
					// Defining duplications
					////////////////////////////+++++++++++++++++++++
					if (segmentB.refStart < segmentA.refEnd) {
						//cout << "Pos:" << pos << endl;
						//cout << "segmentB.refStart:" << segmentB.refStart << endl;
						//cout << " segmentA.refEnd" << segmentA.refEnd << endl;
						startPos = genomicStart+segmentB.refStart+1;
						endPos   = genomicStart+segmentA.refEnd+2;
						size     = endPos - startPos+1;
						type     = "DUP";
					}
				}
				// Defining inversions
				if (segmentA.strand == "+" && segmentB.strand == "-") {
					startPos = genomicStart+segmentA.contigEnd;
					endPos   = genomicStart+segmentB.contigEnd;
					size     = endPos - startPos;
					type     = "INV";
				}
				if (segmentA.strand == "-" && segmentB.strand == "+") {
					startPos = genomicStart+segmentA.contigStart-1;
					endPos   = genomicStart+segmentB.contigStart-1;
					size     = endPos - startPos;
					type     = "INV";
				}
				if ( size >= 100) {
					vcf_t call;
					call.chr =  chr;
					call.start = startPos;
					call.end = endPos;
					call.precision = "PRECISE";
					call.svtype = type;
					call.breakReads = totalReads;
					call.assembled  = totalReadsAssembled;
					call.mapq =  meanMapQ;
					call.hasRDsupport = false;
					call.RDratio = ".";
					call.RDmad   = ".";
					call.LOHsupport = ".";
					call.kdiv = kmerDiv;

					double pvalue_upstream;
					double pvalue_downstream;
					double pvalue_twosided;
					int counts_5prime, counts_inner, counts_3prime;
					double meanPvalue = 0.00;
					call.covPvalue      = meanPvalue;
					call.discordants    = 0;
					call.alleleBalance  = 0;
					call.discPvalue     = 0;
					call.cumulativeSize = 0;
					//call.reciprocalSupport = "no";
					vcf_v.push_back(call);
					if (debug) {
						cout << "MIDA2 VARCALL: "<<  call.start << "\t" << call.end << "\t" << call.svtype << endl;
					}
				}
			}
			// Defining inversions
			//cout << splitReadVec[i].contigStart << "\t" << splitReadVec[i].contigEnd << "\t" << splitReadVec[i].cigar << "\t" << splitReadVec[i].refStart << "\t" << splitReadVec[i].refEnd << endl;
			//cout << splitReadVec[i+1].contigStart << "\t" << splitReadVec[i+1].contigEnd << "\t" << splitReadVec[i+1].cigar << "\t" << splitReadVec[i+1].refStart << "\t" << splitReadVec[i+1].refEnd << endl;
		}
	}
	return vcf_v;
 }
 
 //######################################
 // Retornmen a una idea inicial que teniem:
 // Realitzar aliniaments locals a tots els segments possibles del contig tinguent en compte que:
 // 1. Realitzarem un alineament local per cada segment que falti a cada read, acceptant un mínim score/matches
 // 2. Per acceptar el contig necessitem que es cobreixi > 90% de la longitud del contig
 // 3- S'anirà annotant cada segment mapejat ordenat segons posició relativa en el contig, i es farà el variant call per cada parella contigua de segments.

std::vector<vcf_t> alignContig::Align(int& totalReads, int& totalReadsAssembled, double& kmerDiv, int& meanMapQual) {

	// Get reverse complement of the reference sequence
	std::string revSequence = revComp(reference);

	std::vector<vcf_t> vcf_v;

	for (int i = 0; i < contigs.size();i++) {

		string contig = contigs[i];

		if (debug) {
			cout << " Resolent contig " << contig << "\t" << genomicStart<<  endl;
		}

		bool isFullyChecked = false;

		// First unaligned block is the entire contig sequence
		unalignedBlock wblock;
		wblock.start = 0;
		wblock.end   = contig.length()-1;

		wblock.seq   = contig.substr (wblock.start, contig.length());

		wblock.length= contig.length();
		std::vector<unalignedBlock> unalignedBlockVec;
		unalignedBlockVec.push_back(wblock);
		int refLen = reference.length();
		std::vector<alignmentStruct> splitReadVec;
		int alignedSize = 0;		

		while (isFullyChecked == false) {
		
			unalignedBlock seqBlock = unalignedBlockVec[0];
			
			// Declares a default Aligner
			StripedSmithWaterman::Aligner ssw_aligner(2, 20, 18, 18);

			// Declares a default filter
			StripedSmithWaterman::Filter filter;

			// Declares an alignment that stores the result
			StripedSmithWaterman::Alignment localAlignment;

			// Aligns the query to the ref
			ssw_aligner.Align(seqBlock.seq.c_str(), reference.c_str(), reference.length(), filter, &localAlignment);

			string strand = "+";
			alignmentStruct fwd = returnAlignment(localAlignment, seqBlock.length, refLen, reference, seqBlock.seq, strand);
		
			// El mateix pero pel reverse-complement
			// Aligns the query to the Ref-complement ref
			strand = "-";
			ssw_aligner.Align(seqBlock.seq.c_str(), revSequence.c_str(), revSequence.length(), filter, &localAlignment);
	
			alignmentStruct rev = returnAlignment(localAlignment, seqBlock.length, refLen, reference, seqBlock.seq, strand);

			alignedSize+= treatAligment(fwd, rev, splitReadVec, unalignedBlockVec, contig, reference, readLength);

			unalignedBlockVec.erase(unalignedBlockVec.begin());
			if (unalignedBlockVec.size() == 0) {
				isFullyChecked = true;
			}
		}

		int sizeRV = splitReadVec.size();
		std::sort(splitReadVec.begin(), splitReadVec.end());
			
		float csize = (float)alignedSize/contig.length();
		//cout << "Mida per si acceptem:" << csize << endl;
		//cout << "Mida vector segments:" << sizeRV << endl;
		//cout << "Contig length:" << contig.length() << "\tAligned_bases:" << alignedSize << endl;
		//cout << "Proporció del contig mapejada:" << csize << endl;
		std::vector<vcf_t> vars;
		if (csize > 0.9 && sizeRV > 1) {
			vars = detectSV(splitReadVec, totalReads, totalReadsAssembled, kmerDiv, sizeRV, meanMapQual, contig, reference);
			vcf_v.insert(vcf_v.end(), vars.begin(), vars.end());
		}

		if (sizeRV <= 1 || vars.size() == 0) {

			unalignedBlockVec.clear();
			unalignedBlockVec.push_back(wblock);
			splitReadVec.clear();
			alignedSize = 0;	
			isFullyChecked = false;
			while (isFullyChecked == false) {
			
				unalignedBlock seqBlock = unalignedBlockVec[0];
				
				// Declares a default Aligner
				StripedSmithWaterman::Aligner ssw_aligner(2, 8, 12, 0.5);

				// Declares a default filter
				StripedSmithWaterman::Filter filter;

				// Declares an alignment that stores the result
				StripedSmithWaterman::Alignment localAlignment;

				// Aligns the query to the ref
				ssw_aligner.Align(seqBlock.seq.c_str(), reference.c_str(), reference.length(), filter, &localAlignment);

				string strand = "+";
				alignmentStruct fwd = returnAlignment(localAlignment, seqBlock.length, refLen, reference, seqBlock.seq, strand);
			
				// El mateix pero pel reverse-complement
				// Aligns the query to the Ref-complement ref
				strand = "-";
				ssw_aligner.Align(seqBlock.seq.c_str(), revSequence.c_str(), revSequence.length(), filter, &localAlignment);
		
				alignmentStruct rev = returnAlignment(localAlignment, seqBlock.length, refLen, reference, seqBlock.seq, strand);

				alignedSize+= treatAligment(fwd, rev, splitReadVec, unalignedBlockVec, contig, reference, readLength);

				unalignedBlockVec.erase(unalignedBlockVec.begin());
				if (unalignedBlockVec.size() == 0) {
					isFullyChecked = true;
				}
			}
			sizeRV = splitReadVec.size();
			std::sort(splitReadVec.begin(), splitReadVec.end());
			csize = (float)alignedSize/contig.length();
		}
		if (debug) {
			cout << "Proporció del contig mapejada:" << csize << endl;
		}
		//if (csize > 0.9 && sizeRV > 1) {
		if (csize > 0.9) {
			std::vector<vcf_t> vars = detectSV(splitReadVec, totalReads, totalReadsAssembled, kmerDiv, sizeRV, meanMapQual, contig, reference);
			vcf_v.insert(vcf_v.end(), vars.begin(), vars.end());
		}
		//cout << endl;
	}
	return vcf_v;
	
 }


///Per dups
//0-153	0-(final->153)	score:308
//154-305	(inici->29) -180	score:304



