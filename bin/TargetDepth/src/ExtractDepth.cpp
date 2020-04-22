#include <iostream>
#include <string.h>
#include <fstream>
#include <vector>
#include <map>

#include "SeqLib/BamHeader.h"
#include "SeqLib/RefGenome.h"
#include "SeqLib/BamReader.h"
#include "SeqLib/BamRecord.h"
#include "SeqLib/BamWriter.h"
#include "SeqLib/BWAWrapper.h"
#include "SeqLib/BamWalker.h"

#include "ExtractDepth.h"
#include "bed_t.h"
#include "sample_metrics_region_t.h"

using namespace SeqLib;

 int ExtractDepth::getTotalReads() {
	return totalReads;
 }
 int ExtractDepth::getTotalReadsX() {
	return totalReadsX;
 }
 double ExtractDepth::getMeanCov() {
	return meanCov;
 }
 double ExtractDepth::getMeanCovX() {
	return meanCovX;
 }
 double ExtractDepth::getMeanCount() {
	return meanCount;
 }
 double ExtractDepth::getMeanCountX() {
	return meanCountX;
 }
 double ExtractDepth::getMeanInsertSize() {
	return meanIsize;
 }
 double ExtractDepth::getStdDevInsertSize() {
	return stdDevIsize;
 }

 void ExtractDepth::extract ( std::map<std::string, std::vector<sample_metrics_region_t>>& mapSummaryCount, std::map<std::string, std::vector<sample_metrics_region_t>>& mapSummaryCoverage) {

	BamReader br;
	br.Open(sampleNamePath);
	BamRecord r;

	std::ofstream coverageStream;
	std::string covFile = outdir + "/" + sampleName + "_coverage.bed";
	if (extractCoverage == true) {
		coverageStream.open (covFile, std::ofstream::out);
	}

	std::ofstream countStream;
	std::string countFile = outdir + "/" + sampleName + "_counts.bed";
	std::ofstream isizeStream;
	std::string isizeFile = outdir + "/" + sampleName + "_isizes.bed";
	if (extractCounts == true) {
		countStream.open (countFile, std::ofstream::out);
		isizeStream.open (isizeFile, std::ofstream::out);
	}

	totalReads  = 0;
	totalReadsX = 0;
	std::vector<int> cov_vec;
	std::vector<int> count_vec;

	std::vector<int> cov_vecX;
	std::vector<int> count_vecX;

	std::vector<double> isize_vec;
	int sumIsizes= 0;
	int sumCov   = 0;
	int sumCount = 0;

	int sumCovX  = 0;
	int sumCountX= 0;

	meanIsize  = 0.00;
	stdDevIsize= 0.00;

	for (auto& l : RegionsVector ) {

		std::string coordinate = l.chr + ":" + std::to_string(l.start) + "-" + std::to_string(l.end);
		std::string key = l.chr + "\t" + std::to_string(l.start) + "\t" + std::to_string(l.end) + "\t" + l.info + "\t" + std::to_string(l.gc);

		GenomicRegion gr(coordinate, br.Header());
		br.SetRegion(gr);
		int count  = 0;
		int countX = 0;
		int depth  = 0;
		int depthX = 0;
		int multimaps = 0;
		std::vector<int> insert_sizes;
		double mean_isize;
		int sum = 0;
	        int maxISize = 1000;
		double accum = 0.0;

		std::vector<int> Positions;
		int readSize = 0;
		while (br.GetNextRecord(r)) {
			std::string chrom =  r.ChrName();
			totalReads++;

			std::string seq = r.Sequence();
			readSize = seq.length();

			if (l.chr == "23" || l.chr == "chr23" || l.chr == "X" || l.chr == "chrX") { 
				totalReadsX++;
			}
			if (extractCoverage == true) {
				int pos = r.Position();
				Positions.push_back(pos);
			}
			if (extractCounts == true) {
				std::string chrom =  r.ChrName();
			
				if (l.chr == "23" || l.chr == "chr23" || l.chr == "X" || l.chr == "chrX") { 
 			       		countX++;
			        }
				else{
			       		count++;
				}      

			    int Isize = r.FullInsertSize();
				
			    if (Isize < maxISize) { 
			        insert_sizes.push_back(Isize); 
			        sum+=Isize;
			    }
			    if (r.MapQuality() == 0) {
			        multimaps++;
			    }
			}
		}
		if ( extractCoverage == true) {
			for (int i = l.start; i<=l.end; i++) {
				for (auto& p : Positions) {
					if (p <= i && p+readSize >= i) {
						if (l.chr == "chrX" || l.chr == "X" || l.chr == "chr23" || l.chr == "23") {
							depthX++;
						}
						else {
							depth++;
						}
					}
					if (p > i) {
						break;
					}
				}

				if (l.chr == "chrX" || l.chr == "X" || l.chr == "chr23" || l.chr == "23") {
					sumCovX+=depthX;
					cov_vecX.push_back(depthX);
					std::string coord = l.chr + "\t" +  std::to_string(i) + "\t" +  std::to_string(i+1);
					coverageStream << l.chr << "\t" << i << "\t" << i+1 << "\t"  << l.info << "\t" <<  std::to_string(l.gc) << "\t" << depthX << "\n";
				}
				else {
					sumCov+=depth;
					cov_vec.push_back(depth);
					std::string coord = l.chr + "\t" +  std::to_string(i) + "\t" +  std::to_string(i+1);
					coverageStream << l.chr << "\t" << i << "\t" << i+1 << "\t"  << l.info << "\t" <<  std::to_string(l.gc) << "\t" << depth  << "\n";
				}

				depth = 0;
				depthX= 0;
			}
		}
		if (extractCounts == true) {

			mean_isize =  sum / (double)insert_sizes.size();
			double percMultimaps = 0.00;

			if (l.chr == "chrX" || l.chr == "X" || l.chr == "chr23" || l.chr == "23") {
				sumCountX+=countX;
				count_vecX.push_back(countX);
				percMultimaps = (multimaps/ (double) countX) * (double) 100; 
				countStream << l.chr << "\t" << l.start << "\t" << l.end << "\t" << l.info << "\t" <<  std::to_string(l.gc) << "\t" << countX  << "\n";
			}
			else {
				sumCount+=count;
				count_vec.push_back(count);
				percMultimaps = (multimaps/ (double) count) * (double) 100; 
				countStream << l.chr << "\t" << l.start << "\t" << l.end << "\t" << l.info << "\t" <<  std::to_string(l.gc) << "\t" << count  << "\n";
			}
	
			for (auto& i : insert_sizes) {
				sumIsizes+=i;
				isize_vec.push_back(i);
		    		accum += (i - mean_isize) * (i - mean_isize);
			}
			double stdev = sqrt(accum / (double)(insert_sizes.size()-1));
			isizeStream << l.chr << "\t" << l.start << "\t" << l.end << "\t" << l.info << "\t" <<  std::to_string(l.gc) << "\t" << mean_isize << "\t" << stdev  << "\n";
		}
	}
	if (extractCounts == true) {
		meanIsize  =  (double)sumIsizes / isize_vec.size();
		meanCount  = sumCount / (double)count_vec.size();
		meanCountX = sumCountX / (double)count_vecX.size();

		double accum = 0.0;
		for (auto& i : isize_vec) {
		 	accum += (i - meanIsize) * (i - meanIsize);
		}
		stdDevIsize = sqrt(accum / (double)isize_vec.size()-1);

	}
	if (extractCoverage == true) {
		meanCov  = sumCov/(double)cov_vec.size();
		meanCovX = sumCovX/(double)cov_vecX.size();
	}
 }


