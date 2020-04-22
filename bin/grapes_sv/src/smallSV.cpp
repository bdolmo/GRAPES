#include <iostream>
#include <vector>
#include <string>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <omp.h>

// grapes headers
#include "Assembler.h"
#include "Aligner.h"
#include "alignContig.h"
#include "dna2bit.h"
#include "GenomicRange.h"
#include "sam_t.h"
#include "SV.h"
#include "Utils.cpp"
#include "vcf_t.h"
#include "varCall.h"
#include "smallSV.h"
#include "Trimmer.h"

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

/*int getNumClusters ( map<string, vector<sam_t>>& breakReadClusters ) {

	vector<GenomicRange> breakClusters = binarySearch( breakReadLocations[chr], breakReadLocations[chr].size(), chr, start, end );

}*/

void Analyze ( string& bamFile, vector<string>& contigs, string& chr, int& start, string& sctype, int& totalReads, int& totalAssembled, string& refSequence, std::vector<vcf_t>& vcf_v, std::string& wes, double& kmerDiversity ) {

	std::vector<sam_t> v_sam;

	for (int i = 0; i < contigs.size(); i++) {

		BamRecord r;
		BamRecordVector brv;
		bool hardclip = false;
		float secondary_cutoff = 0.90; // secondary alignments must have score >= 0.9*top_score
		int secondary_cap = 10; // max number of secondary alignments to return

		// Make an in-memory BWA-MEM index of region
		BWAWrapper bwa;
		UnalignedSequenceVector usv = {{chr, refSequence}};
		bwa.ConstructIndex(usv);
		string num = std::to_string(i);
		string Qname = chr + "_" + num;
		BamRecordVector results; // alignment results (can have multiple alignments)
		
		bwa.SetAScore(2);
		bwa.SetGapOpen(30);
		bwa.SetGapExtension(1);			
		bwa.SetBandwidth(1000);	
		bwa.SetMismatchPenalty(4);		
		bwa.SetReseedTrigger(1.5);
		bwa.Set3primeClippingPenalty(5);
		bwa.Set5primeClippingPenalty(5);
			
		BamHeader myHeader;
		myHeader = bwa.HeaderFromIndex();

		bwa.AlignSequence(contigs[i], Qname, results, hardclip, secondary_cutoff, secondary_cap);
		vector<string> bwa_results;

		for (auto& i : results) {

			sam_t sam_read;
			sam_read.read_name = i.Qname();
			sam_read.chr       = myHeader.IDtoName(i.ChrID());
			sam_read.pos       = i.Position();
			sam_read.align_pos = i.AlignmentPosition();
			sam_read.align_end = i.AlignmentEndPosition();
			sam_read.cigar     = i.CigarString();
			sam_read.strand    = i.ReverseFlag() == true ? '-' : '+';
			sam_read.order     = i.FirstFlag()   == true ? '1' : '2';
			int mapq_ascii     = i.MapQuality();
			sam_read.mapq      = mapq_ascii;
			sam_read.seq       = i.Sequence();
			sam_read.clipped   = i.NumClip();
	
			v_sam.push_back(sam_read);
		 }
	 }

	 SV call(v_sam, chr, chr, start, start, sctype, totalReads, totalAssembled, bamFile, refSequence, "small", 0, 0.00, kmerDiversity);
	 std::string undetermined = "undetermined";
	 if (call.classify(vcf_v, undetermined) == 1) {

	 }
	 else   {
		//cout << "uncalled" << chr << "\t" << start << "\t"  << "\t" << totalReads<< "\tTOTAL_ASSEMBLED:" <<totalAssembled << "\n" << "\n";
	 }
}

void smallSV::callSmallSV (map<string, vector<sam_t>>& breakReadClusters, std::map<std::string, std::vector<GenomicRange>>& breakReadLocations, std::string& fastq, std::vector<vcf_t>& vcf_v, std::string& wes, std::string& outdir) {

   // Goal here is to jointly consider all brek-reads that fall inside +/-500 bp from each soft-clipping position
   // as assembly candidates.

   int N = breakReadClusters.size();
   int count = 0;
	
   string dummy = "INIT";
   dna2bit DNA(dummy);

   RefGenome ref;
   ref.LoadIndex(reference);

   // opening bam file	
   BamReader reader;
   reader.Open(bamFile);

   BamHeader myHeader;
   myHeader = reader.Header();

   std::ofstream Fastq;
   Fastq.open(fastq.c_str(), std::ofstream::app);

   std::multimap<std::string, std::string> visitedBreak;

   std::string logF = outdir + "/" + "grapes.log";
   std::ofstream logFile;
   logFile.open(logF, std::ofstream::app);
   int kmerSize = readLength*0.15;
   cout << " INFO: kmer-size for assemblying: " << kmerSize << endl;

   std::multimap<std::string, std::string> breaksInWindows;

   for (auto& breakC : breakReadClusters) {
	    vector<string> tmpCoord = split ( breakC.first, '\t');
	    string chr = tmpCoord[0];

	    int pos   = std::stoi(tmpCoord[1]);		
	    int start = pos - 250;
	    int end   = pos + 250;

	    vector<GenomicRange> breakClusters = binarySearch( breakReadLocations[chr], breakReadLocations[chr].size(), chr, start, end );
   	    for (auto& c: breakClusters) {

		std::string mkey = chr + "\t" + std::to_string(c.start) + "\t" + c.softclip_type;
		visitedBreak.insert(std::pair<std::string, std::string>(mkey, mkey));
		if (visitedBreak.count(mkey) <= 1) {
			breaksInWindows.insert(std::pair<std::string, std::string>(breakC.first, mkey));
		}
	   }
    }

  	boost::progress_display show_progress( breaksInWindows.size() );

	// Using dynamic schedule for assigning chunks on runtime. This works better when we have unbalanced tasks between threads.
	// http://jakascorner.com/blog/2016/06/omp-for-scheduling.html
	#pragma omp parallel 
	{
	    std::vector<vcf_t> vcf_v_private;
	
	    int tid = omp_get_thread_num();
	
  	    //fill vec_private in parallel
		#pragma omp for schedule(dynamic, 1)
 		//for (std::multimap<string,string>::iterator breakC=breaksInWindows.begin(); breakC!=breaksInWindows.end(); ++breakC) {	
		//for (auto& breakC : breaksInWindows) {
   		 for(int i = 0; i < breaksInWindows.size(); i++) {

			auto breakC = breaksInWindows.begin();
       		advance(breakC, i);

			vector<string> tmpCoord = split ( breakC->first, '\t');
			string chr = tmpCoord[0];

			int pos   = std::stoi(tmpCoord[1]);		
			int start = pos - 500;
			int end   = pos + 500;

			int totalReads;
			int totalReadsAssembled;
			std::vector<string> contigsJoint;
			double kmerDiversity;

			std::vector<std::string> breakReadVec;
			vector<string> tmpClust = split ( breakC->second, '\t');
			std::string mkey = tmpClust[0] + "\t" + tmpClust[1] + "\t" + tmpClust[2];
			std::map<std::string, std::vector<sam_t>>::iterator it1 = breakReadClusters.find(mkey);
			int meanMapQual;	

			if (debug) {
				cout << mkey << endl;
			}

			if(it1 != breakReadClusters.end()) {

				std::vector<int> mapqQualities;

				for (auto& read: it1->second) {
					string seq  = DNA.decompressDNA(read.bitSeq, read.seqLen);

					mapqQualities.push_back(read.mapq);
					Trimmer trim1(seq, read.qual);
					string trimmedseq = trim1.trim();

					if ( trimmedseq.length() >= kmerSize ) {
						breakReadVec.push_back(trimmedseq);
					}							
					if (debug) {
						cout << "Read:" << trimmedseq << endl;
					}	
				}

				meanMapQual  = accumulate( mapqQualities.begin(), mapqQualities.end(), 0.0)/mapqQualities.size(); 
				Assembler A(breakReadVec, minOlapBases, 2, 2, kmerSize);
				std::vector<std::string> contigs = A.ungappedGreedy();
				int numAssembledA = A.getNumAssembled();

				for (auto& c : contigs) {
					//cout << "COntig assemblat: " << c << endl;
				}
				//cout << numAssembledA << endl;
				totalReads = breakReadVec.size();
				totalReadsAssembled = numAssembledA;
				kmerDiversity = getKmerDiversity(contigs);
				if (totalReadsAssembled > 0)  {
					for (auto& seq : contigs ) {
						if (debug) {
							cout << "CONTIG:" << seq << "\t" << kmerDiversity << "\t" << totalReadsAssembled << endl;
						}
						contigsJoint.push_back(seq);
					}
				}
			}
			if (debug) {
				cout << "Fora surt" << endl;
			}

			std::vector<vcf_t>smallSV;

			if (contigsJoint.size() > 0) {
				RefGenome ref;
 				ref.LoadIndex(reference);

        		int lengthChrom = myHeader.GetSequenceLength(chr) ;

				start = start >= 0 ? start : 0;
				end   = end <= lengthChrom ? end : lengthChrom;

				std::string localReference = ref.QueryRegion(chr, start, end);

				alignContig A (contigsJoint, localReference, chr, start, end, readLength, reference);
				smallSV = A.Align(totalReads, totalReadsAssembled, kmerDiversity, meanMapQual);
				for (auto& b : smallSV) {
					vcf_v_private.push_back(b);
				}
			}
			#pragma omp critical 
			{
				++show_progress;
			}
			//}*/
				// Recovering additional SVs using bwamem
				if (smallSV.size() == 0) {
					int num = 0;
					for (auto& seq : contigsJoint ) {
						#pragma omp critical
						{
							Fastq << "@r" << num << "_" << chr << "_" << pos << "_" << tmpClust[2] << "_" << totalReads << "_" << totalReadsAssembled << "_" << kmerDiversity << endl;
							Fastq << seq << endl;
							Fastq << "+" << endl;
							Fastq << fillQualString(seq) << endl;
						}
						++num;
					}
				}
				contigsJoint.clear();
				totalReadsAssembled = 0;
				totalReads = 0;
			 //}*/
		//}

	}
	#pragma omp critical
		vcf_v.insert(vcf_v.end(), vcf_v_private.begin(), vcf_v_private.end());
    }
   	//cout << vcf_v.size() << endl;
}

void smallSV::callSmallSvExome (map<string, vector<sam_t>>& breakReadClusters, std::map<std::string, std::vector<GenomicRange>>& breakReadLocations, std::string& fastq, std::vector<vcf_t>& vcf_v, std::string& wes, std::string& outdir) {

   // Goal here is to jointly consider all brek-reads that fall inside +/-500 bp from each soft-clipping position
   // as assembly candidates.

   int N = breakReadClusters.size();
   int count = 0;
	
   string dummy = "INIT";
   dna2bit DNA(dummy);

   RefGenome ref;
   ref.LoadIndex(reference);

   std::ofstream Fastq;
   Fastq.open(fastq.c_str(), std::ofstream::app);

   std::multimap<std::string, std::string> visitedBreak;

   std::string logF = outdir + "/" + "grapes.log";
   std::ofstream logFile;
   logFile.open(logF, std::ofstream::app);

   for (auto& breakC : breakReadClusters) {

	    vector<string> tmpCoord = split ( breakC.first, '\t');
	    string chr = tmpCoord[0];

	    int pos   = std::stoi(tmpCoord[1]);		
	    int start = pos - 500;
	    int end   = pos + 500;

	    vector<GenomicRange> breakClusters = binarySearch( breakReadLocations[chr], breakReadLocations[chr].size(), chr, start, end );

	    int totalReads;
	    int totalReadsAssembled;
        std::vector<string> contigsJoint;
	    double kmerDiversity;

	    std::vector<std::string> leftCluster;
	    std::vector<std::string> rightCluster;

	    for (auto& read : breakC.second) {

		  string Seq = DNA.decompressDNA(read.bitSeq, read.seqLen);

		  // Trimming low quality bases
		  Trimmer trim1(Seq, read.qual);
		  Seq = trim1.trim();
		  //int total_trimmed = trim1.getTotalTrimmed();

		  if (read.sctype == "LEFT" || read.sctype == "DEL" || read.sctype == "INS") {
			leftCluster.push_back(Seq);
		  }
		  if (read.sctype == "RIGHT") {
			rightCluster.push_back(Seq);
		  }
	    }

	    int totalAssembledLeft;
	    int totalAssembledRight;
	    vector<string> contigsLeft;
	    vector<string> contigsRight;

	    chr = returnChromLexicoGraphical ( chr );
	    string refSequence = ref.QueryRegion(chr, start, end);

    	if (leftCluster.size() > 0) {

			Assembler assemblyLeft(leftCluster, minOlapBases, 3, 2, 20);
			vector<string> contigsLeft = assemblyLeft.ungappedGreedy();

			totalAssembledLeft   = assemblyLeft.getNumAssembled();

			int totalReadsLeft   = leftCluster.size();
			double kmerDiversity = getKmerDiversity(contigsLeft);
			int num = 0;
			string softClass = "LEFT";
			Analyze (bamFile, contigsLeft, chr, pos, softClass , totalReadsLeft, totalAssembledLeft, refSequence, vcf_v, wes, kmerDiversity);
	    }
		if (rightCluster.size() > 0) {

			Assembler assemblyRight(rightCluster, minOlapBases, 3, 2, 20);
			vector<string> contigsRight = assemblyRight.ungappedGreedy();

			totalAssembledRight = assemblyRight.getNumAssembled();

			int totalReadsRight  = rightCluster.size();
			double kmerDiversity = getKmerDiversity(contigsRight);
			int num = 0;
			string softClass = "RIGHT";
			Analyze (bamFile, contigsRight, chr, pos,softClass , totalReadsRight, totalAssembledRight, refSequence, vcf_v, wes, kmerDiversity);
		}
	    leftCluster.clear();
	    rightCluster.clear();
     }
}


