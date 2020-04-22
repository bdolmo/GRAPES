#include <string>
#include <vector>
#include <iostream>
#include <map>
#include <regex>
#include <string.h>
#include "SeqLib/RefGenome.h"
#include "SeqLib/BamReader.h"
#include "SeqLib/BamRecord.h"
#include "ssw_cpp.h"
#include "clusterSC.h"
#include "sam_t.h"
#include "Utils.cpp"
#include "GenomicRange.h"
#include "excludeRegion.h"
#include "dna2bit.h"
#include "alignContig.h"
#include <boost/progress.hpp>

using namespace SeqLib;

std::string IntToString ( int& );

int clusterSC::getReadLength(){
	return readLength;
}

inline alignmentStruct returnAlignment (const StripedSmithWaterman::Alignment&,  int&, int&, string&, string&, string&);

int whichReadLength ( std::string & bamFile ) {

	// opening bam file	
	BamReader reader;
	reader.Open(bamFile);
	BamRecord r;
	int counter = 0;
	vector<int> array;
	while (reader.GetNextRecord(r)) {
		counter++;
		array.push_back(r.Length());
		if (counter > 500) {
			break;
		}
	}
	int readLength = mostFrequentPosition(array);
	return readLength;
}


float GetSoftClipQual ( std::string& cigar, std::string& sctype, std::string& seq, std::string& qual ) {

	int softclipped, matched;
	string clippedSeq, clippedQual;
	int sum = 0;
 	vector<string> s;
 	vector<string> m;
	if (sctype == "LEFT") {
		s = split(cigar, 'S');
		m = split(s[1], 'M');
		softclipped = atoi(s[0].c_str());
		matched     = atoi(m[0].c_str());
		clippedSeq  = seq.substr(0, softclipped);
		clippedQual = qual.substr(0, softclipped);
	}
	if (sctype == "RIGHT") {
		m = split(cigar, 'M');
		s = split(m[1], 'S');
		matched     = atoi(m[0].c_str());
		softclipped = atoi(s[0].c_str());
		clippedSeq  = seq.substr(matched);
		clippedQual = qual.substr(matched);
	}
	for (auto base : clippedQual) {
		// transforming to ASCII 33-offset
		int asci_decimal = base-33 >= 0 ? base-33 : 0;
		sum += asci_decimal;	
	}
	float avg = ((float)sum)/clippedQual.length();
	return avg;	
}

alignmentStruct rescueSplitReads ( BamRecord& r, string& reference, string& bamFile, int& lengthChromo ) {

	int chr = r.ChrID();
	// opening bam file	
	BamReader reader;
	reader.Open(bamFile);
	BamHeader myHeader;
	myHeader = reader.Header();
	string chrom =  myHeader.IDtoName(r.ChrID());
	chrom = returnChromLexicoGraphical( chrom );

	int Start = r.Position()-250 >= 0 ? r.Position()-250 : r.Position();
	int End   = r.Position()+250 <= lengthChromo ? r.Position()+250 : lengthChromo;

	RefGenome ref;
 	ref.LoadIndex(reference);

	int minEdge;	

	string localReference = ref.QueryRegion(chrom, Start, End);

	// Declares a default Aligner
	StripedSmithWaterman::Aligner ssw_aligner(2,15, 25, 25);

	// Declares a default filter
	StripedSmithWaterman::Filter filter;

	// Declares an alignment that stores the result
	StripedSmithWaterman::Alignment localAlignment;
	string readSequence = r.Sequence();
	int refLen = localReference.length();

	// Aligns the query to the ref
	ssw_aligner.Align(readSequence.c_str(), localReference.c_str(), refLen, filter, &localAlignment);

	string strand = "+";
	alignmentStruct fwd = returnAlignment(localAlignment, refLen, refLen, localReference, readSequence, strand);

	fwd.refStart = Start+fwd.refStart;	
	fwd.refEnd = Start+fwd.refEnd;	

	int matches = 0;
	int clipped = 0;
	int globalPosition = 0;
	std::string concat;
	for (auto& i : fwd.cigar)  {
		if (i == 'S' or i == 'M') {
			if (i == 'S') {
				int n = std::stoi(concat);
				clipped+=n;
			}	
			if (i == 'M') {
				int n = std::stoi(concat);
				matches +=n;	
			}
			globalPosition+=std::stoi(concat);	
			concat = "";	
		}
		else  {
			concat += i;
		}
	}
	fwd.numMatched = matches;	
	fwd.numClip = clipped;
	if (!fwd.refStart) {
		fwd.refStart = r.Position();
		fwd.refEnd   = r.Position() + r.Length();
		fwd.cigar    = r.CigarString();
	}
	return fwd;
}

int getMismatches( std::vector<std::string>& MDtag) {

	int mismatches = 0;
	bool isIndel = false;

	for (auto&element : MDtag) {

		std::string tmp = "";
		//int readPos = 0;

		for (int i=0; i<element.length(); i++) {

			if (element[i] == 'A' || element[i] == 'C' || element[i] == 'T' || element[i] == 'G' || element[i] == '^') {

				if (tmp == "") {
					tmp = '0';
				}
				int tmpInt = std::stoi(tmp);
				//readPos+= std::stoi(tmp);

				if (MDtag.size() > 1) {
					isIndel = true;
				}

				if (isIndel == false) { 
					mismatches++;
				}
				tmp = "";
			}
			else {
				isIndel = false;
				tmp = tmp + element[i];
			}
		}
	}	
	return mismatches;
}

void clusterSC::extractSC( std::map<std::string, std::vector<sam_t>>& breakReadClusters, std::map<std::string, std::vector<GenomicRange>>&  map_break_ranges, std::map<std::string, std::vector<GenomicRange>>& mapOfRegions, std::string& reference ) {

	// opening bam file	
	BamReader reader;
	reader.Open(bamFile);

	BamHeader myHeader;
	myHeader = reader.Header();
	BamRecord r;

	std::string softclip_type;
    std::smatch m;
	//std::regex soft_regex_L ("^[0-9]+S[0-9]+M$");
	//std::regex soft_regex_R ("^[0-9]+M[0-9]+S$");

	std::regex soft_regex_L ("^[0-9]+S.*");
	std::regex soft_regex_R (".*[0-9]+S$");

	std::regex soft_regex_DEL ("^[0-9]+M[0-9]+D[0-9]+M$");
	std::regex soft_regex_INS ("^[0-9]+M[0-9]+I[0-9]+M$");

	int edge_pos, pos_A = 0;
	std::string chrA; 
	std::string chrB;

	std::string coordinate; 
	std::string coordLeft; 
	std::string coordRight; 
	std::string coordIndel;
	
	int flag = 0;
	int breakReads = 0;

	std::vector<sam_t> localBreakIndel;
	std::vector<sam_t> localBreakLeft;
	std::vector<sam_t> localBreakRight;

	std::cout << " INFO: Clustering BreakReads " << "\n";

   	std::string clustSR = outDir + "/" + "clusters_SR.bed";
	ofstream sr_clust;
	sr_clust.open (clustSR);

	readLength = whichReadLength(bamFile);
	int minSoftLength = (softClipLenPerc/100) * (float) readLength;
	int counter = 0;

	//std::map<std::string, std::string> progressChr;

    std::cout << " INFO: Minimum soft-clipped reads required: " << minSoftClusterSize << "\n";

  	boost::progress_display show_progress( totalSR );

	while (reader.GetNextRecord(r)) {
		++show_progress;
		counter++;
		//float p = (counter / (float) totalSR) * (float) 100;
		//printProgBar(p);


		//cout << r.Position() << "\t" << r.CigarString() << endl;


		if (r.ChrID() < 0 || r.ChrID() >= 25 || r.MateChrID()  >= 25 ) { 
			continue; 
		}
		if (r.ChrName() == "chrM" || r.ChrName() == "MT") { 
			continue; 
		}

		if (r.CountBWASecondaryAlignments() > 3) { 
			//cout << "Salta" << endl;
			continue; 
		}

		std::string tag = "MD";
		std::vector<std::string> MDtag = r.GetSmartStringTag(tag);
		
		// Skipping reads with more than 5 mismatches
		//int numMismatches = r.Length() - r.NumClip() -r.NumMatchBases();
		int numMismatches  = 0;

		if (r.MaxDeletionBases() == 0 && r.MaxInsertionBases()== 0) {
			numMismatches = getMismatches(MDtag);
		}

		if (numMismatches > 6) {
			//cout << numMismatches << "\t"<< "Salta2" << endl;

			//std::cout << "Mismatches: " << numMismatches << "\t" << r.Qname() << "\t"<< r.Position() << "\t" << r.CigarString() << std::endl;
			continue;
		}

		// skipping reads with mapq less than 20
		if (r.MapQuality()  < 10 ) { 
			continue; 
		} 

		// skipping PCR duplicates
		if (r.DuplicateFlag() == true ) { 
			//cout << "Salta3" << endl;

			continue; 
		}

		// BreakReads
		//if (r.NumSoftClip() > 0 && r.NumSoftClip() < minSoftLength ) { 
		//	cout << "Salta10\n";

		//	continue; 
		//}
		// Min Deleted bases required
		//if (r.MaxDeletionBases() > 0 && r.MaxDeletionBases() < 20) { 
			//continue; 
		//}
		// Min Inserted bases required
		//if (r.MaxInsertionBases() > 0 && r.MaxInsertionBases() < 20) { 
			//continue; 
		//}

		int chromosome = r.ChrID();
		string chrom =  myHeader.IDtoName(r.ChrID());
		chrom = returnChromLexicoGraphical( chrom );
        int lengthChrom = myHeader.GetSequenceLength(chrom) ;
	
		int readPos = r.Position();
		int readLen = r.Length();
		int readPosEnd = readPos+readLen;

		bool isOnRegion = liesOnBlackRegion( mapOfRegions[chrom], chrom, readPos, readPosEnd);

		if (isOnRegion) {
			continue;
		}

		if (r.CountNBases() > 0 || mapOfRegions.count( chrom ) == 0) {
			totalSR++;
			continue;
		} 

		sam_t break_read;
		break_read.read_name = r.Qname();
		break_read.chr       = chrom;
		break_read.pos       = r.Position();
		break_read.align_pos = r.AlignmentPosition();
		break_read.align_end = r.AlignmentEndPosition();
		break_read.cigar     = r.CigarString();

		break_read.strand    = r.ReverseFlag() == true ? '-' : '+';
		break_read.order     = r.FirstFlag()   == true ? '1' : '2';
		int mapq_ascii       = r.MapQuality();
		break_read.mapq      = mapq_ascii;
		break_read.isRevComp = false;

		int orientation = r.PairOrientation();
		string Sequence = r.Sequence();

		break_read.seqLen     = Sequence.length();
		dna2bit DNA(Sequence);
		uint8_t* DNA_C = DNA.compressDNA(Sequence);

		break_read.bitSeq    = DNA_C;
		break_read.qual      = r.Qualities();
		break_read.clipped   = r.NumClip();
		break_read.numMatched= r.NumMatchBases();
		std::string tmpCigar = r.CigarString();

		//std::cout << r.Qname() << "\t"<< r.Position() << "\t" << r.CigarString() << "\t" << numMismatches << "\t" << totalSR << "\t" << counter << std::endl;
		// Attempt to rescue informative split-reads at complex regions
		/*if (numMismatches > 0 || r.MaxDeletionBases() > 0 || r.MaxInsertionBases() > 0) {

			alignmentStruct alignment =  rescueSplitReads ( r, reference, bamFile, lengthChromo); 
			break_read.pos       = alignment.refStart;
			//break_read.align_pos = alignment();
			//break_read.align_end = r.AlignmentEndPosition();
			break_read.clipped   = alignment.numClip;
			break_read.numMatched= alignment.numMatched;
			break_read.cigar     = alignment.cigar;
			tmpCigar= alignment.cigar;
		}*/	

		if (r.NumClip() >= 10 || r.MaxDeletionBases() >= 10 || r.MaxInsertionBases() >= 10) {

			std::string tmps = tmpCigar;

			// gapped alignments
			if (r.MaxDeletionBases() >= 10 || r.MaxInsertionBases() >= 10) {

				if (std::regex_search (tmps, m, soft_regex_DEL) || std::regex_search (tmpCigar, m, soft_regex_INS))  {

	    			if (std::regex_search (tmps, m, soft_regex_DEL)) { 
						softclip_type = "DEL"; 
					}
					if (std::regex_search (tmps, m, soft_regex_INS)) { 
						softclip_type = "INS"; 
					}

					break_read.sctype  = softclip_type;
					std::vector<std::string> tmp = split( r.CigarString(), 'M' );
					int st = atoi(tmp[0].c_str());
					pos_A = r.Position() + st;

					if (flag == 0 || coordIndel == "") {
						edge_pos = r.Position() + st;
						flag = 1;
						coordIndel =  break_read.chr + "\t" + std::to_string(pos_A) + "\t" + softclip_type;
					}
					if (pos_A >= edge_pos-10  && pos_A < edge_pos +10) {
						breakReads++;
						localBreakIndel.push_back(break_read);

						GenomicRange grange;
						grange.chromosomeA = chrA;
						grange.start = edge_pos;
						grange.end   = edge_pos+1;
						grange.softclip_type = softclip_type;
					    map_break_ranges[chrA].push_back(grange);
					}
				}
			}

			// clipped alignments
			if (std::regex_search (tmpCigar, m, soft_regex_L)) {
				
				softclip_type = "LEFT";

				break_read.sctype  = softclip_type;
				pos_A = break_read.pos;
				break_read.coordinate = break_read.chr + "\t" + std::to_string(pos_A) + "\t" + softclip_type;

				if (flag == 0 ) {
					coordLeft = break_read.chr + "\t" + std::to_string(pos_A) + "\t" + softclip_type;
					chrA  =  break_read.chr;
					edge_pos  = pos_A;
					flag = 1;
				}
				float avgQual = GetSoftClipQual ( break_read.cigar , softclip_type, Sequence, break_read.qual );
				if (avgQual > 15 && pos_A >= edge_pos -10  && pos_A < edge_pos +10 && chrA == break_read.chr) {

					breakReads++;
					localBreakLeft.push_back(break_read);
					chrB =  break_read.chr;
					GenomicRange grange;
					grange.chromosomeA = chrA;
					grange.start = edge_pos;
					grange.end   = edge_pos+1;
					grange.softclip_type = "LEFT";
					coordLeft = break_read.chr + "\t" + std::to_string(edge_pos) + "\t" + softclip_type;

					if (orientation == 1 && r.NumHardClip() == 0 || orientation == 3 && r.NumHardClip() == 0) {
						break_read.isRevComp = true;
					}
					if (chrA == chrB) {

						if (debug) {
							//cout << chrA << "\t" << pos_A << "\t" << pos_A+1 << "\t" << "LEFT" << "\t" << Sequence << endl;
						}

						map_break_ranges[chrA].push_back(grange);
					}
				}
			} 
			else if (std::regex_search (tmpCigar, m, soft_regex_R)) {

				softclip_type = "RIGHT";
				break_read.sctype = softclip_type;
				pos_A = break_read.pos + break_read.numMatched;
				break_read.coordinate =  break_read.chr + "\t" + std::to_string(pos_A) + "\t" + softclip_type;

				if (flag == 0 ) {
					coordRight = break_read.chr + "\t" + std::to_string(pos_A) + "\t" + softclip_type;
					chrA  =   break_read.chr;
					edge_pos = pos_A;
					flag = 1;
				}
				float avgQual = GetSoftClipQual ( break_read.cigar , softclip_type, Sequence, break_read.qual );
				if (avgQual > 15 && pos_A >= edge_pos-10 && pos_A < edge_pos +10 && chrA == break_read.chr) {

					breakReads++;
					localBreakRight.push_back(break_read);
					chrB =  break_read.chr;
					GenomicRange grange;
					grange.chromosomeA = chrA;
					grange.start = edge_pos;
					grange.end   = edge_pos+1;
					grange.softclip_type = "RIGHT";
					coordRight = break_read.chr + "\t" + std::to_string(edge_pos) + "\t" + softclip_type;
					if (chrA == chrB) {
						if (debug) {
						//	cout << chrA << "\t" << pos_A << "\t" << pos_A+1 << "\t" << "RIGHT" << "\t" << Sequence << endl;
						}
						
						map_break_ranges[chrA].push_back(grange);
					}
				}
				if (orientation == 1 && r.NumHardClip() == 0 || orientation == 3 && r.NumHardClip() == 0) {
					break_read.isRevComp = true;
				}
			}

			// Clustering soft-clips within a 20bp window 
			if (pos_A < edge_pos-10 || pos_A > edge_pos +10 || chrA != chrB || counter == totalSR) {

				// This avoids incorporating soft-clips at the very beggining of the chromosomes
				int minEdge = pos_A-500;

				if (breakReads >= minSoftClusterSize && breakReads < 200 && minEdge > 0) {
					if (localBreakIndel.size() >= minSoftClusterSize) {
						if (debug) {
						//	cout << "ENTRANT: " <<  coordIndel << "\t" << breakReads << "\t" << localBreakIndel.size() << endl;
						}
						sr_clust << coordIndel << "\t" << breakReads << "\t" << localBreakIndel.size() << endl;
						breakReadClusters.insert(std::pair<std::string, std::vector<sam_t>>(coordIndel, localBreakIndel));
					}
					if (localBreakLeft.size() >= minSoftClusterSize) {

						if (debug) {
						//	cout  << "ENTRANT: " <<  coordLeft << "\t" <<  breakReads << "\t" << localBreakLeft.size() << "\t" << localBreakRight.size() << endl;
						}						
						sr_clust << coordLeft << "\t" <<  breakReads << "\t" << localBreakLeft.size() << "\t" << localBreakRight.size() << endl;
						breakReadClusters.insert(std::pair<std::string, std::vector<sam_t>>(coordLeft,  localBreakLeft));
					}
					if (localBreakRight.size() >= minSoftClusterSize) {

						if (debug) {
						//	cout << "ENTRANT: " << coordRight << "\t" << breakReads << "\t" << localBreakLeft.size() << "\t" << localBreakRight.size()<< endl;
						}	
						sr_clust << coordRight << "\t" << breakReads << "\t" << localBreakLeft.size() << "\t" << localBreakRight.size()<< endl;
						breakReadClusters.insert(std::pair<std::string, std::vector<sam_t>>(coordRight, localBreakRight));
					}
				}

				localBreakLeft.clear();
				localBreakIndel.clear();
				localBreakRight.clear();

				if (softclip_type == "LEFT") {
					coordLeft  = break_read.chr + "\t" + std::to_string(pos_A) + "\t" + softclip_type;
					localBreakLeft.push_back(break_read);
				}
				else if (softclip_type == "RIGHT") {
					coordRight = break_read.chr + "\t" + std::to_string(pos_A) + "\t" + softclip_type;
					localBreakRight.push_back(break_read);
				}
				else if (softclip_type == "DEL") {
					coordIndel = break_read.chr + "\t" + std::to_string(pos_A) + "\t" + softclip_type;
					localBreakIndel.push_back(break_read);
				}
				else if (softclip_type == "INS") {
					coordIndel = break_read.chr + "\t" + std::to_string(pos_A) + "\t" + softclip_type;
					localBreakIndel.push_back(break_read);
				}

				flag = 0;
				breakReads = 1;				
				edge_pos = pos_A;
				chrA = chrB;
			}
	    }
	}
	breakReads = 0;
}


