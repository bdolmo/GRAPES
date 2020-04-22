#include <iostream>
#include <vector>
#include <string>
#include <stdlib.h>
#include <string.h>
#include <math.h>

// grapes headers
#include "Assembler.h"
#include "Aligner.h"
#include "dna2bit.h"
#include "GenomicRange.h"
#include "sam_t.h"
#include "SV.h"
#include "Trimmer.h"
#include "Utils.cpp"
#include "vcf_t.h"
#include "varCall.h"
#include "largeSV.h"
#include "callSV.h"

// Seqlib headers
#include "SeqLib/RefGenome.h"
#include "SeqLib/BamReader.h"
#include "SeqLib/BWAWrapper.h"
#include "SeqLib/BamRecord.h"
#include "SeqLib/BFC.h"
#include "SeqLib/FermiAssembler.h"

#include <boost/progress.hpp>


using namespace std;
using namespace SeqLib;

std::string returnSoftEdge(std::string& svtype, std::string& location, std::string& softEdge) {
	
	if (svtype == "INV") {
		return softEdge;
	}
	if (svtype == "DEL" && location == "UP") {
		softEdge = "RIGHT";
	}
	if (svtype == "DUP" && location == "UP") {
		softEdge = "LEFT";
	}
	if (svtype == "DEL"  && location == "DOWN") {
		softEdge = "LEFT";
	}
	if (svtype == "DUP" && location == "DOWN") {
		softEdge = "RIGHT";
	}
	return softEdge;
}

// This function return a cluster if it contains at least (-r) minimum soft-clips.
// This is because at clusterSC.cpp we clustered SC reads with different subtypes
//std::vector<std::string> checkClusterMinSoft

void largeSV::callStructVar( std::map<std::string, std::vector<GenomicRange>>& discordantList, std::map<std::string, std::vector<GenomicRange>>& breakReadLocations, std::map<std::string, std::vector<sam_t>>& breakReadClusters, std::vector<vcf_t>& vcf_v, std::string svtype, std::string softEdge) {

	int N = 0;
	int count = 0;

	string dummy = "INIT";
	dna2bit DNA(dummy);

	RefGenome ref;
	ref.LoadIndex(reference);

	string genomeLoc;
   	std::map<std::string, std::string> visitedBreak;

	for (auto& x : discordantList) {
		for (auto& y : x.second) {
			N++;
		}
	}
  	boost::progress_display show_progress( N );

 	for (auto& chrom : discordantList) {

	   	for (auto& coord : chrom.second) {

			count++;
			++show_progress;

			vector<GenomicRange> upstream_v;
			vector<GenomicRange> downstream_v;
			std::vector<std::string> cluster;

			// Defining upstream, downstream and midstream searching coordinates for soft-clip cluster catching
			std::string upstream_ref;
			std::string downstream_ref;
			std::string pseudo_ref;

			string chr = returnChromLexicoGraphical( chrom.first );	

			int start_pseudo, end_pseudo, addNtdPositions; 
	
			map<std::string, std::vector<sam_t>>::iterator it1;
			map<std::string, std::vector<sam_t>>::iterator it2;
			map< map<std::string, std::vector<sam_t>>::iterator,  map<std::string, std::vector<sam_t>>::iterator> iterators_map;

			std::string svlength;
				
			// The idea here is to add soft-clips to every discordant cluster and perform local assembly
			// For small/medium SVs
			if (coord.end - coord.start <= 2000 ) {
				svlength = "medium";
				int n = coord.end-coord.start;
				int m = coord.end-coord.start;
				int diff = (n)/2;
				int loffset = (m)/2;

				start_pseudo    = coord.start- diff;
				end_pseudo      = coord.end  + diff;
				addNtdPositions = end_pseudo - start_pseudo;

				// Finding soft-clip clusters using binary search
				upstream_v   = binarySearch( breakReadLocations[chrom.first], breakReadLocations[chrom.first].size(), chrom.first, coord.start-loffset, coord.start+loffset );
				downstream_v = binarySearch( breakReadLocations[chrom.first], breakReadLocations[chrom.first].size(), chrom.first, coord.end-loffset+1, coord.end+loffset );
				pseudo_ref   = ref.QueryRegion(chr, start_pseudo, end_pseudo);
			}
			// For larger SVs
			else {
				svlength = "large";
				// Looking for breakpoints at upstream region
				upstream_v   = binarySearch( breakReadLocations[chrom.first], breakReadLocations[chrom.first].size(), chrom.first, coord.start-500, coord.start+500);

				// Looking for breakpoints at downstream region
				downstream_v   = binarySearch( breakReadLocations[chrom.first], breakReadLocations[chrom.first].size(), chrom.first, coord.end-500, coord.end+500);
				start_pseudo   = coord.start-500;
				end_pseudo     = coord.end-500;
				upstream_ref   = ref.QueryRegion(chrom.first, coord.start-500, coord.start+500);
				downstream_ref = ref.QueryRegion(chrom.first, coord.end-500, coord.end+500);
				pseudo_ref     = upstream_ref + downstream_ref;
			}
			// Now selecting most appropiate upstream soft-clusters and assemble them
			genomeLoc = "UP";
			softEdge = returnSoftEdge(svtype, genomeLoc, softEdge);
			for (auto& c: upstream_v) {
				std::string mkey = chr + "\t" + std::to_string(c.start) + "\t" + c.softclip_type;

				if (visitedBreak.count(mkey)>0) {
					continue;
				}
				if (softEdge != "undef" && softEdge != c.softclip_type) {
					continue;
				}
		 		it1 = breakReadClusters.find(mkey);
				if(it1 != breakReadClusters.end()) {

					for (auto& read: it1->second) {
						string seq  = DNA.decompressDNA(read.bitSeq, read.seqLen);
						Trimmer trim1(seq, read.qual);
						seq = trim1.trim();
						cluster.push_back(seq);
					}
				}
			}
			// Now selecting most appropiate downstream soft-clusters and assemble them 
			genomeLoc = "DOWN";
			softEdge = returnSoftEdge(svtype, genomeLoc, softEdge);
			for (auto& c: downstream_v) {

				std::string mkey = chr + "\t" + std::to_string(c.start) + "\t" + c.softclip_type;

				if (visitedBreak.count(mkey)>0) {
					continue;
				}
				if (softEdge != "undef" && softEdge != c.softclip_type) {
					continue;
				}
		 		it2 = breakReadClusters.find(mkey);
				if(it2 != breakReadClusters.end()) {

					for (auto& read: it2->second) {

						string seq  = DNA.decompressDNA(read.bitSeq, read.seqLen);
						Trimmer trim1(seq, read.qual);
						seq = trim1.trim();
						int nTrim = trim1.getTotalTrimmed();

						if (svtype == "INV") {
							cluster.push_back(revComp(seq));
						}
						else {
							cluster.push_back(seq);
						}
					}
				}
			}

			if (cluster.size() >  0) {

				Assembler A(cluster, minOlapBases, 2, 2, 20);
				std::vector<std::string> contigs = A.ungappedGreedy();
				int numAssembledA = A.getNumAssembled();
	
				int totalAssembled = numAssembledA;
				double kmerDiversity = 0.0;
				kmerDiversity = getKmerDiversity(contigs);
			    int totalReads= cluster.size();

				int i = 0;
				vector<sam_t> v_sam;
				for (auto& x: contigs) {

					BamRecord r;
					BamRecordVector results; // alignment results (can have multiple alignments)
					bool hardclip = false;
					float secondary_cutoff = 0.95; // secondary alignments must have score >= 0.9*top_score
					int secondary_cap = 5; // max number of secondary alignments to return

					// Make an in-memory BWA-MEM index of region
					BWAWrapper bwamem;
					UnalignedSequenceVector usv = {{chrom.first, pseudo_ref}};
					bwamem.ConstructIndex(usv);
					string Qname = chrom.first + "_" + std::to_string(i);

					bwamem.SetAScore(2);
					bwamem.SetGapOpen(30);
					bwamem.SetGapExtension(0.5);			
					bwamem.SetBandwidth(1000);	
					bwamem.SetMismatchPenalty(8);		
					bwamem.SetReseedTrigger(1.5);
					bwamem.Set3primeClippingPenalty(5);
					bwamem.Set5primeClippingPenalty(5);
					bwamem.AlignSequence(contigs[i], Qname, results, hardclip, secondary_cutoff, secondary_cap);

					for (auto& i : results) {
		
						sam_t sam_read;
						sam_read.read_name = i.Qname();
						sam_read.chr       = i.ChrName();
						sam_read.pos       = i.Position();
						sam_read.align_pos = i.AlignmentPosition();
						sam_read.align_end = i.AlignmentEndPosition();
						sam_read.cigar     = i.CigarString();
						sam_read.strand    = i.ReverseFlag() == true ? '-' : '+';
						sam_read.order     = i.FirstFlag()   == true ? '1' : '2';
						int mapq_ascii     = i.MapQuality();
						sam_read.mapq      = mapq_ascii;
						sam_read.seq       = i.Sequence();
						v_sam.push_back(sam_read);
					}
					i++;
				}
				SV call(v_sam, chr, chr, start_pseudo, end_pseudo, "undef",  totalReads, totalAssembled, bamFile, pseudo_ref, svlength,  coord.supporting, coord.pvalue_discordant, kmerDiversity);
				int foundSV = call.classify(vcf_v, svtype);
				if (foundSV) {

					// Erasing used clusters
					genomeLoc = "UP";
					softEdge = returnSoftEdge(svtype, genomeLoc, softEdge);
					for (auto& c: upstream_v) {

						if (softEdge != "undef" && softEdge != c.softclip_type) {
							continue;
						}
						std::string mkey =  chrom.first + "\t" + std::to_string(c.start) + "\t" + c.softclip_type;	
						visitedBreak.insert(std::pair<std::string, std::string>(mkey, mkey));
						breakReadClusters.erase(mkey);
					}
					genomeLoc = "DOWN";
					softEdge = returnSoftEdge(svtype, genomeLoc, softEdge);
					for (auto& c: downstream_v) {

						if (softEdge != "undef" && softEdge != c.softclip_type) {
							continue;
						}
						std::string mkey =  chrom.first + "\t" + std::to_string(c.start) + "\t" + c.softclip_type;	
						visitedBreak.insert(std::pair<std::string, std::string>(mkey, mkey));
						breakReadClusters.erase(mkey);
			 		}
				 }
				 else {

					int len = coord.end - coord.start;
					vcf_t call;
					call.chr =  chrom.first;
					call.start = coord.start;
					call.end = coord.end;
					call.precision = "IMPRECISE";
					call.svtype = svtype;
					call.breakReads = 0;
					call.assembled  = 0;
					call.mapq =  coord.mapq;
					call.hasRDsupport = false;
					call.RDratio = ".";
					call.RDmad   = ".";
					call.LOHsupport = ".";

					double pvalue_upstream;
					double pvalue_downstream;
					double pvalue_twosided;
					int counts_5prime, counts_inner, counts_3prime;
					double meanPvalue = 0.00;
					call.covPvalue      = meanPvalue;
					call.discordants    = coord.supporting;
					call.alleleBalance  = coord.discordant_ratio;
					call.discPvalue     = coord.pvalue_discordant;
					call.cumulativeSize = coord.cumulativeSize;
					//call.reciprocalSupport = "no";
					vcf_v.push_back(call);
				 }
			    }
			    else {
					int len = coord.end - coord.start;
					vcf_t call;
					call.chr =  chrom.first;
					call.start = coord.start;
					call.end = coord.end;
					call.precision = "IMPRECISE";
					call.svtype = svtype;
					call.breakReads = 0;
					call.assembled = 0;
					call.mapq =  coord.mapq;
					call.hasRDsupport = false;
					call.RDratio = ".";
					call.RDmad   = ".";
					call.LOHsupport = ".";

					double pvalue_upstream;
					double pvalue_downstream;
					double pvalue_twosided;
					int counts_5prime, counts_inner, counts_3prime;
					double meanPvalue  = 0.00;
					call.covPvalue     = meanPvalue;
					call.discordants   = coord.supporting;
					call.alleleBalance = coord.discordant_ratio;
					call.discPvalue    = coord.pvalue_discordant;
					call.cumulativeSize= coord.cumulativeSize;
					//call.reciprocalSupport = "no";
					vcf_v.push_back(call);
			    }
			}
		}
}
