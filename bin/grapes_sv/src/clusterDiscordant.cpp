#include <string>
#include <vector>
#include <iostream>
#include <map>
#include <regex>
#include <string.h>
#include <numeric>

#include "SeqLib/BamHeader.h"
#include "SeqLib/RefGenome.h"
#include "SeqLib/BamReader.h"
#include "SeqLib/BamRecord.h"
#include "SeqLib/BamWriter.h"
#include "SeqLib/BWAWrapper.h"
#include "SeqLib/BamWalker.h"
#include "ssw_cpp.h"

#include "alignContig.h"
#include "clusterDiscordant.h"
#include "sam_t.h"
#include "pe_adjacency.h"
#include "GenomicRange.h"
#include "ReadPairInfo.h"
#include "excludeRegion.h"
#include "Utils.cpp"

using namespace SeqLib;

float poisson_pmf(int , double);
void buildPairInfo ( const std::string& , const std::string& , const int&, const int&, const int& , const std::string&, const std::string&, const int& , const int&, const bool&, const int&, std::map<std::string, ReadPairInfo>& );
double getCov ( std::string& bamFile, string chr, int start, int end);

inline alignmentStruct returnAlignment (const StripedSmithWaterman::Alignment&,  int&, int&, string&, string&, string&);
inline int getMismatches( std::vector<std::string>& );
inline alignmentStruct rescueSplitReads ( BamRecord&, string&, string&, int&);


 float clusterDiscordant::getMeanCoverage() {
	return mean_coverage;
 }
 int clusterDiscordant::getTotalFR() {
	return totalFR;
 }
 int clusterDiscordant::getTotalRF() {
	return totalRF;
 }
 int clusterDiscordant::getTotalFF() {
	return totalFF;
 }
 int clusterDiscordant::getTotalRR() {
	return totalRR;
 }
 int clusterDiscordant::getBinSize() {
	return binSize;
 }
 std::vector<int> clusterDiscordant::getGerminalCounts() {
	return readCounts_germinal;
 }
 std::vector<int> clusterDiscordant::getSomaticCounts() {
	return readCounts_somatic;
 }
 int clusterDiscordant::getTotalSR() {
	return totalSR;
 }
 int clusterDiscordant::getReadLength() {
	return readLength;
 }

double  getCov( std::string& bamFile, string chr, int start, int end) {

	int st_A = start-50;
	int st_B = start-49;
	int ed_A = end-70;
	int ed_B = end-50;
	std::string upstream_area, downstream_area;

	upstream_area   = chr + ":" + std::to_string(st_A) + "-" + std::to_string(st_B);
	double cov = -1;

	BamReader br;
	br.Open(bamFile);
	BamRecord r;

	int count = 0;
	GenomicRegion gr(upstream_area, br.Header());
	br.SetRegion(gr);
	while (br.GetNextRecord(r)) {
      		count++;
		if (count > 200) { break; }
	}
	if (count == 0) {
		cov = 0.00;
	}
	else {
		cov = count;
	}
	return cov;		
}

std::tuple<double, double> clusterDiscordant::calculateInsertSizeLimits () {

	BamReader reader;
	reader.Open(bamFile);
	BamRecord r;
	int count = 0;

	BamHeader myHeader;
	myHeader = reader.Header();

	std::vector<int> isizes;
	int sum = 0;
	const int maxISize = 1000;
	double accum = 0.0;
	
	int read_counts = 0;
	std::vector<int> positions;
	while (reader.GetNextRecord(r) ) {
	
		if (r.ChrID() < 1 || r.ChrID() >= 25 || r.MateChrID()  >= 25 ) { 
			continue; 
		}

		std::string chrom = myHeader.IDtoName(r.ChrID());
		if (chrom == "chrM" || chrom == "MT" || chrom == "M" ) { continue; }

		count++;
		if (count > 5000000) {
			break;
		}

		readLength = r.Length();
		read_counts++;
		positions.push_back(r.Position());
		if (r.ProperPair() && r.ProperOrientation())  {  
			count++;
			int insert_size = r.InsertSize() < 0 ? r.InsertSize() *-1 : r.InsertSize();
			if (insert_size >= maxISize) { 
				continue; 
			}
			sum+=insert_size;
			isizes.push_back(insert_size);
		}
	}

	mean =  sum / isizes.size();

	int distance = positions[positions.size()-1] - positions[0];

	mean_coverage = ( (double)read_counts/(double)distance )*readLength;

	std::cout << " INFO: mean coverage: " << mean_coverage << "\n";

	for (auto& i : isizes) {
    	accum += (i - mean) * (i - mean);
	}
	stdev = sqrt(accum / (isizes.size()-1));
	return std::make_tuple(mean, stdev );
}


void clusterDiscordant::cluster( std::string& bam, std::map<std::string, std::vector<GenomicRange>>& map_of_ranges ) {

	BamReader discordantBam;
	discordantBam.Open(bam);
	BamRecord r;

	BamHeader myHeader;
	myHeader = discordantBam.Header();

	std::string out_cluster;
    std::ofstream outfile;

	std::string line;

	std::regex soft_regex_L ("^[0-9]+S[0-9]+M$.*");
	std::regex soft_regex_R ("^[0-9]+M[0-9]+S$.*");
	int countDisc;
	if (bam ==  outdir + "/" +  sampleName + ".FR.bam") {
		std::cout << " INFO: Clustering FR pairs" << "\n";
		out_cluster =  outdir + "/" +  sampleName + "." + "FR.clusters.bed";
		countDisc = totalFR;
	}
	if (bam ==  outdir + "/" +  sampleName + ".RF.bam") {
		std::cout << " INFO: Clustering RF pairs" << "\n";
		out_cluster =  outdir + "/" +  sampleName + "." +  "RF.clusters.bed";
		countDisc = totalRF;
	}
	if (bam ==  outdir + "/" +  sampleName + ".FF.bam") {
		std::cout << " INFO: Clustering FF pairs" << "\n";
		out_cluster =  outdir + "/" +  sampleName + "." +  "FF.clusters.bed";
		countDisc = totalFF;
	}
	if (bam ==  outdir + "/" +  sampleName + ".RR.bam") {
		std::cout << " INFO: Clustering RR pairs" << "\n";
		out_cluster =  outdir + "/" +  sampleName + "." +  "RR.clusters.bed";
		countDisc = totalRR;
	}

    outfile.open (out_cluster);

	int flag    = 0;
	int trigger = 0;
	string first_chrom, second_chrom, first_chrom_mate, second_chrom_mate;
	int first_pos, second_pos, first_end, second_end;
	int start, end;

	std::vector<int> mapq_v;

	std::vector<PE_adjacency> adjacencies;
	PE_adjacency adj;

	std::multimap<std::string, std::string> mapReadName;

	int totalDisc = 0;

	std::smatch m;
	while (discordantBam.GetNextRecord(r)) {

	    BamRecordVector results;
		mapReadName.insert(std::pair<std::string, std::string>(	r.Qname(), r.Qname()));
		std::string readName = r.Qname();
		totalDisc++;
		

		if ( readPair_map[r.Qname()].first_mate_chr == "" || !readPair_map[r.Qname()].first_mate_start || !readPair_map[r.Qname()].second_mate_start  || readPair_map[r.Qname()].second_mate_chr == "") { 
			countDisc = countDisc-1;
			totalDisc = totalDisc-1;
	
			if (totalDisc < countDisc) {
				continue; 
			}
		}		
		if (flag == 0) {

			first_chrom      = myHeader.IDtoName(r.ChrID());
			first_chrom_mate = myHeader.IDtoName(r.MateChrID());

			int readPos = r.Position();
			int matePos = r.MatePosition();
			std::string tmpCigar = r.CigarString();

			if (std::regex_search (tmpCigar, m, soft_regex_R)) {
				readPos = r.Position() + r.NumMatchBases();
			}

			first_pos = readPos < matePos ? readPos : matePos;
			first_end = readPos > matePos ? readPos : matePos;

			flag = 1;
			adj.chromosome = first_chrom;
			adj.start      = first_pos;
			adj.end        = first_end;
			adj.supporting = 1;
			//continue;
		}
		if (flag == 1) {

			std::string readName        = r.Qname();
			std::string chromo_first    = readPair_map[r.Qname()].first_mate_chr;
			std::string chromo_second   = readPair_map[r.Qname()].second_mate_chr;
			int start_first             = readPair_map[r.Qname()].first_mate_start;
			int start_second            = readPair_map[r.Qname()].second_mate_start;
			int mapQ_first              = readPair_map[r.Qname()].first_mate_mapq;
			int mapQ_second             = readPair_map[r.Qname()].second_mate_mapq;
			int alignEnd_first          = readPair_map[r.Qname()].first_mate_alignEnd;
			int alignEnd_second         = readPair_map[r.Qname()].second_mate_alignEnd;
			std::string cigar_first     = readPair_map[r.Qname()].first_mate_cigar;
			std::string cigar_second    = readPair_map[r.Qname()].second_mate_cigar;
			std::string first_strand    = readPair_map[r.Qname()].first_mate_strand;
			std:: string second_strand  = readPair_map[r.Qname()].second_mate_strand;

			std::map<std::string, std::string> mapFirstStrand;
			std::map<std::string, std::string> mapSecondStrand;

			second_chrom      = myHeader.IDtoName(r.ChrID());
			second_chrom_mate = myHeader.IDtoName(r.MateChrID());

			int readPos = r.Position();
			int matePos = r.MatePosition();
			std::string tmpCigar = r.CigarString();

			if (std::regex_search (tmpCigar, m, soft_regex_R)) {
				readPos = r.Position() + r.NumMatchBases();
			}

			second_pos = readPos < matePos ? readPos : matePos;
			second_end = readPos > matePos ? readPos : matePos;
			
			std::string chromA =  second_chrom;
			std::string chromB =  second_chrom_mate;

			chromA = returnChromLexicoGraphical( chromA );
			chromB = returnChromLexicoGraphical( chromB );

			int st_first, st_second;
			int align_ed_first, align_ed_second;
			std::string cigar_st, cigar_ed;
			std::string st_strand, ed_strand;

			if ( start_first < start_second) {

				st_first = start_first;
				st_second = start_second;
				align_ed_first = alignEnd_first;
				align_ed_second = alignEnd_second;
				cigar_st = cigar_first;
				cigar_ed = cigar_second;
				st_strand = first_strand;
				ed_strand = second_strand;

				mapFirstStrand.insert(std::pair<std::string, std::string>(first_strand, first_strand));
				mapSecondStrand.insert(std::pair<std::string, std::string>(second_strand, second_strand));
			}
			else {
				st_first = start_second;
				st_second = start_first;
				align_ed_first = alignEnd_second;
				align_ed_second = alignEnd_first;
				cigar_st = cigar_second;
				cigar_ed = cigar_first;
				st_strand = first_strand;
				ed_strand = second_strand;

				mapFirstStrand.insert(std::pair<std::string, std::string>(second_strand, second_strand));
				mapSecondStrand.insert(std::pair<std::string, std::string>(first_strand, first_strand));
			}
		
			if (st_first == st_second) {
				countDisc = countDisc-1;
				if (totalDisc < countDisc) {
					continue; 
				}
			}

			if ( (first_chrom == second_chrom) && (first_chrom_mate == second_chrom_mate) && (second_pos >= first_pos) && (second_pos <= first_pos+upper_limit+450) ) {
				if (mapReadName.count(readName)< 2) {
					adj.supporting ++;
				}

				mapq_v.push_back(mapQ_first);
				mapq_v.push_back(mapQ_first);

				if (bam ==  outdir + "/" +  sampleName + ".FR.bam") {
					start = st_first + align_ed_first;
					end   = st_second;
					adj.align_begins.push_back(start);
					adj.align_ends.push_back(end);
					adjacencies.push_back(adj);
				}
				if (bam ==  outdir + "/" +  sampleName + ".RF.bam") {
					start = st_first;
					end   = st_second + align_ed_second;
					adj.align_begins.push_back(start);
					adj.align_ends.push_back(end);
				}
				if (bam ==  outdir + "/" +  sampleName + ".FF.bam") {
					start = st_first + align_ed_first;
					end   = st_second + align_ed_second;
					adj.align_begins.push_back(start);
					adj.align_ends.push_back(end);
				}
				if (bam ==  outdir + "/" +  sampleName + ".RR.bam") {
					start = st_first;
					end   = st_second;
					adj.align_begins.push_back(start);
					adj.align_ends.push_back(end);
				}
   		       }
		       if ( (first_chrom != second_chrom) || (first_chrom_mate != second_chrom_mate) || (second_pos < first_pos) || (second_pos > first_pos+upper_limit+450) || (totalDisc >= countDisc) ) {

				// Imposing some restrictions
				if (adj.supporting >= minClusterSize && adj.supporting < 200) {
			
					std::sort(std::begin(adj.align_begins), std::end(adj.align_begins));
					std::sort(std::begin(adj.align_ends), std::end(adj.align_ends));

					int anchor1_start = adj.align_begins[0];
					int anchor1_end   = adj.align_begins.back() + r.Length();
					int anchor2_start = adj.align_ends[0];
					int anchor2_end   = adj.align_ends.back() + r.Length();	
					
					int span_anchor1 = anchor1_end - anchor1_start;
					int span_anchor2 = anchor2_end - anchor2_start;	
					int cumulativeSize = span_anchor1 + span_anchor2;

					start = computeMedian(adj.align_begins);
					end = computeMedian(adj.align_ends);

					//double t_cov = getCov(bamFile, adj.chromosome, start, end);
					//double discordant_ratio = adj.supporting/t_cov;
					double t_cov = 50;
					double discordant_ratio = 0.5;

					// Hard filtering to avoid low fraction of discordant reads
					if (discordant_ratio > 0.05) {
						std::string chrom = first_chrom;
						int tmpSt, tmpEd;
						if (start > end) { 
							tmpSt = end;
							tmpEd = start;
							start = tmpSt;
							end = tmpEd;
						}
						chrom  = returnChromLexicoGraphical( chrom );
						chromA = returnChromLexicoGraphical( chromA );
						chromB = returnChromLexicoGraphical( chromB );
													
						int mean_mapq = std::accumulate(mapq_v.begin(), mapq_v.end(), 0.0) / mapq_v.size();
						mapq_v.clear();

						GenomicRange grange;
						grange.chromosomeA = chrom;
						grange.chromosomeB = chrom;
						grange.start= start;
						grange.end  = end;
						grange.t_cov= t_cov;
						grange.supporting = adj.supporting;
						grange.discordant_ratio = discordant_ratio;
						grange.cumulativeSize = cumulativeSize;
						grange.mapq = mean_mapq;
	
						map_of_ranges[chrom].push_back(grange);	
						outfile << chrom << "\t" << start << "\t" << end << "\t" << adj.supporting << "\t" << t_cov << "\t" << discordant_ratio << "\t" << grange.pvalue_discordant << "\t" << cumulativeSize << "\t" << mean_mapq << "\n";
					}
				}

				first_chrom      = myHeader.IDtoName(r.ChrID());
				first_chrom_mate = myHeader.IDtoName(r.MateChrID());
				first_pos = r.Position();
				first_end = r.MatePosition();
				adjacencies.clear();
				mapFirstStrand.clear();
				mapSecondStrand.clear();
				adj.supporting = 1;
				adj.align_begins.clear();
				adj.align_ends.clear();

				adj.chromosome = myHeader.IDtoName(r.ChrID());
				adj.start = r.Position();
				adj.end = r.MatePosition();
				adj.align_begins.push_back(r.Position());
				adj.align_ends.push_back(r.MatePosition());
				adjacencies.push_back(adj);
			}
		}
	}
    mapReadName.clear();
    outfile.close();
}

bool is_file_exist(const char *fileName) {
    std::ifstream infile(fileName);
    return infile.good();
}

void obtainSRinfo ( std::string& bamFile, int& count ) {

	BamReader reader;
	reader.Open(bamFile);
	BamRecord r;

	BamHeader myHeader;
	myHeader = reader.Header();
	reader.Header();
	std::string chrom;
	while (reader.GetNextRecord(r)) {
		count++;		
	}
}

void populatePairInfo( std::string& bamFile, std::map<std::string, ReadPairInfo>& readPair_map, int& count ) {
	BamReader reader;
	reader.Open(bamFile);
	BamRecord r;

	BamHeader myHeader;
	myHeader = reader.Header();
	reader.Header();
	std::string chrom;
	while (reader.GetNextRecord(r)) {

   		BamRecordVector results; 

		chrom = myHeader.IDtoName(r.ChrID());
		chrom = returnChromLexicoGraphical( chrom );

		int alignPos = r.AlignmentPosition();
		int alignEnd = r.AlignmentEndPosition();
		int insert_size = r.InsertSize() < 0 ? r.InsertSize() *-1 : r.InsertSize();

		if (r.AlignmentPosition() < 0) {
			alignPos = r.AlignmentPosition() * -1;
		}
		if (r.AlignmentEndPosition() < 0) {
			alignEnd = r.AlignmentEndPosition() * -1;
		}
		std::string strand = r.ReverseFlag() == true ? "-" : "+";

		count++;

		buildPairInfo( r.Qname(), chrom, r.Position(), r.MatePosition(), r.MapQuality(), r.CigarString(), strand, alignPos, alignEnd, r.FirstFlag(), r.CountBWAChimericAlignments(), readPair_map);
	}
}

void clusterDiscordant::extractDiscordant( std::map<std::string, std::vector<GenomicRange>>& mapOfRegions) {

	// opening bam file	
	BamReader reader;
	reader.Open(bamFile);
	BamRecord r;

	BamHeader myHeader;
	myHeader = reader.Header();

	reader.Header();
	
	std::string FR_file = outdir + "/" +  sampleName + ".FR.bam";
	std::string RF_file = outdir + "/" +  sampleName + ".RF.bam";
	std::string FF_file = outdir + "/" +  sampleName + ".FF.bam";
	std::string RR_file = outdir + "/" +  sampleName + ".RR.bam";
	std::string SR_file = outdir + "/" +  sampleName + ".SR.bam";
	std::string cnvFile = outdir + "/" +  sampleName + ".CNV.bed";


	bool existsFR = false;
	bool existsRF = false;
	bool existsFF = false;
	bool existsRR = false;
	bool existsSR = false;
	bool existsCNV= false;

    if (is_file_exist(FR_file.c_str()) ) {
		existsFR = true;
		cout << " INFO: Populating FR Pair info" << endl;
		populatePairInfo(FR_file, readPair_map, totalFR);
	}
	if ( is_file_exist(RF_file.c_str()) ) { 
		existsRF = true;
		cout << " INFO: Populating RF Pair info" << endl;
		populatePairInfo(RF_file, readPair_map, totalRF);
	}
	if (is_file_exist(FF_file.c_str())  ) {
		existsFF = true;
		cout << " INFO: Populating FF Pair info" << endl;
		populatePairInfo(FF_file, readPair_map, totalFF);
	}
	if (is_file_exist(RR_file.c_str())  ) {
		existsRR = true;
		cout << " INFO: Populating RR Pair info" << endl;
		populatePairInfo(RR_file, readPair_map, totalRR);
	}
	if (is_file_exist(SR_file.c_str())) {
		existsSR = true;
		cout << " INFO: Populating SR info" << endl;
		obtainSRinfo(SR_file, totalSR);
	}
	if (is_file_exist(cnvFile.c_str())) {
		existsCNV = true;
	}

	if ( existsFR && existsFF && existsRR && existsRF && existsSR) {
		return;
	}

	BamWriter FR; 
	if (existsFR == false ) {
		FR.Open( FR_file );
		FR.SetHeader(reader.Header());
		FR.WriteHeader();
	}

	BamWriter RF; 
	if (existsRF == false ) {
		RF.Open(RF_file);
		RF.SetHeader(reader.Header());
		RF.WriteHeader();
	}

	BamWriter FF; 
	if (existsFF == false ) {
		FF.Open(FF_file);
		FF.SetHeader(reader.Header());
		FF.WriteHeader();
	}

	BamWriter RR; 
	if (existsRR == false ) {
		RR.Open(RR_file);
		RR.SetHeader(reader.Header());
		RR.WriteHeader();
	}

	BamWriter SR; 
	if (existsSR == false ) {
		SR.Open(SR_file);
		SR.SetHeader(reader.Header());
		SR.WriteHeader();
	}

	std::tie( mean, stdev) =  calculateInsertSizeLimits();
	upper_limit = mean + (stdev*numSDs); 
	
	std::cout << " INFO: Mean insert size: " << mean << " bp\tStandard Deviation: " << stdev << "\tUpper limit: " << upper_limit << "\n";

	RefGenome ref;
	ref.LoadIndex(genome);

	int reads_window = 0;
	string chr_window;
	int start_window;
	int flag = 0;

	int is_zero_based = 0;

        std::ofstream cov_windows;
	if (getCountsByWindow == "on") {
 		cov_windows.open (outdir + "/" + sampleName + ".counts_window.bed");
	}
	
	std::vector<int> mapq_v;

	totalDiscordants = 0;
	if (mean_coverage > 0 && mean_coverage <= 2.5) {
		binSize = 10000;
	}
	if (mean_coverage > 2.5 && mean_coverage <= 5) {
		binSize = 2500;
	}
	if (mean_coverage > 5 && mean_coverage <= 10) {
		binSize = 1000;
	}
	if (mean_coverage > 10 && mean_coverage <= 20) {
		binSize = 500;
	}
	if (mean_coverage > 20 && mean_coverage <= 30) {
		binSize = 250;
	}	
	if (mean_coverage > 30) {
		binSize = 100;
	}
	int binCount = 0;
	genomeSize = 1;

    if (getCountsByWindow == "on") {
		std::cout << " INFO: Extracting Read Counts on consecutive windows (" << binSize << " bp) and Breakpoint signals " << "\n";
	}
	else {
		std::cout << " INFO: Extracting Breakpoint signals " << "\n";
	}
    
	string chrom;
	while (reader.GetNextRecord(r)) {

   		BamRecordVector results; 

		if (r.ChrID() < 0 || r.ChrID() >= 25 || r.MateChrID()  >= 25 ) { 
			continue; 
		}
		chrom = myHeader.IDtoName(r.ChrID());
		chrom = returnChromLexicoGraphical( chrom );

        int lengthChrom = myHeader.GetSequenceLength(chrom) ;

		if (chrom == "chrM" || chrom == "MT" || chrom == "M" ) { continue; }

		// skipping PCR duplicates
		if (r.DuplicateFlag() == true ) { continue; }

		if (r.CountNBases() > 0) { continue; } 

		if (r.PairedFlag() == false) { continue; }
		
		int alignPos = r.AlignmentPosition();
		int alignEnd = r.AlignmentEndPosition();
		int insert_size = r.InsertSize() < 0 ? r.InsertSize() *-1 : r.InsertSize();
		int readPos = r.Position();
		int readLen = r.Length();
		int readPosEnd = readPos+readLen;

		bool isOnRegion = liesOnBlackRegion( mapOfRegions[chrom], chrom, readPos, readPosEnd);

		if (isOnRegion) {
			continue;
		}
		if (r.NumClip() >= 10 || r.MaxDeletionBases() >= 10 || r.MaxInsertionBases() >= 10 ) {
			totalSR++;

			int numMismatches  = 0;
			std::string tag = "MD";
			std::vector<std::string> MDtag = r.GetSmartStringTag(tag);

			if (r.MaxDeletionBases() == 0 && r.MaxInsertionBases()== 0) {
				numMismatches = getMismatches(MDtag);
			}
			if (numMismatches > 6) {
				continue;
			}
			if (numMismatches >= 1 || r.MaxDeletionBases() > 0 || r.MaxInsertionBases() > 0) {
				alignmentStruct alignment =  rescueSplitReads ( r, genome, bamFile, lengthChrom ); 
				r.SetPosition(alignment.refStart);
				Cigar reCigar(alignment.cigar);
				r.SetCigar(reCigar);
			}
			if (existsSR == false ) {
 				SR.WriteRecord(r);
			}
		}

		if (insert_size > maxInsertSize) {
			continue;
		}

		if (insert_size >= upper_limit) {
			totalDiscordants++;
		}
		if ( getCountsByWindow == "on") {

			// Counting reads on a window
			if (flag == 0) {
				chr_window   = myHeader.IDtoName(r.ChrID());
				start_window = r.Position();
				flag = 1;
			}
			if (flag == 1) {
				if (myHeader.IDtoName(r.ChrID()) == chr_window && r.Position() <= start_window + binSize  && r.CountNBases() == 0) {
					reads_window++;
					mapq_v.push_back(r.MapQuality());
				}
				else if ( (myHeader.IDtoName(r.ChrID()) == chr_window) && (r.Position() > start_window+binSize)) {
				
					int computeMedian_mapq = computeMedian(mapq_v);
					mapq_v.clear();
					chrom = myHeader.IDtoName(r.ChrID());
					chrom = returnChromLexicoGraphical( chrom );
					if ( isChrSomatic( chrom ) ) {
						readCounts_somatic.push_back(reads_window);
					}
					if ( isChrX( chrom ) ) {
						readCounts_germinal.push_back(reads_window); 
					}
					float gc_content = computeGC( chrom, start_window, start_window + binSize, ref);
					if (existsCNV == false ) {
						cov_windows << chrom << "\t" <<  start_window << "\t" << start_window + binSize << "\t" << reads_window << "\t" << gc_content << "\t" << computeMedian_mapq <<  "\n";
					}
					reads_window = 0;
					chr_window   = myHeader.IDtoName(r.ChrID());
					start_window = start_window + binSize +1;
					binCount++;
				}
				else if (myHeader.IDtoName(r.ChrID()) != chr_window ) {
					chr_window   = myHeader.IDtoName(r.ChrID());
					start_window = r.Position();
				}
			}
		}
		// use only uniquely aligned reads
		if (r.CountBWASecondaryAlignments() > 3) { 
			continue; 
		}
		std::string strand = r.ReverseFlag() == true ? "-" : "+";
		// skipping reads with mapq less than 20
		if (r.MapQuality() < 10) { 
			//continue;
		} 		
		// For inter-chromosomal pairs. 
		if (r.Interchromosomal() == true) {
			if (r.MapQuality() < 35) {
				continue; 
			} 
			//buildPairInfo( r.Qname(), r.ChrID(), r.Position(), r.MatePosition(), r.MapQuality(), r.CigarString(), strand, alignPos, alignEnd, r.FirstFlag(), r.CountBWAChimericAlignments(), readPair_map);
			//INT.WriteRecord(r);
			//continue;
		}
		// For intra-chromosomal pairs
		else {

			std::string chrom = myHeader.IDtoName(r.ChrID());
			if (r.AlignmentPosition() < 0) {
				alignPos = r.AlignmentPosition() * -1;
			}
			if (r.AlignmentEndPosition() < 0) {
				alignEnd = r.AlignmentEndPosition() * -1;
			}

			int ori = r.PairOrientation();

			// FR read pair
			if (ori == 0) {
				if (insert_size >= upper_limit) {

					vector<GenomicRange> found_intersect;

					if ( isChrSomatic( chrom ) ) {
						readCounts_somatic.push_back(reads_window);
					}
					else {
						readCounts_germinal.push_back(reads_window); 
					}

					bool isOnRegion = liesOnBlackRegion( mapOfRegions[chrom], chrom, readPos, readPosEnd);

					if (isOnRegion) {
						continue;
					}
					buildPairInfo( r.Qname(), chrom, r.Position(), r.MatePosition(), r.MapQuality(), r.CigarString(), strand, alignPos, alignEnd, r.FirstFlag(), r.CountBWAChimericAlignments(), readPair_map);
					if (existsFR == false ) {
						r.SetPosition(r.Position());
						r.SetCigar(r.CigarString());
			 			FR.WriteRecord(r);
					}				
					totalFR++;

					continue;
				}
			}
			// RR read pair
			if (ori == 3) {

				if (r.MapQuality() < 20) { 
					continue;
				} 
				if ( isChrSomatic( chrom ) ) {
					readCounts_somatic.push_back(reads_window);
				}
				else {
					readCounts_germinal.push_back(reads_window); 
				}
				if (insert_size == 0) {
					continue;
				}

				bool isOnRegion = liesOnBlackRegion( mapOfRegions[chrom], chrom, readPos, readPosEnd);

				if (isOnRegion) {
					continue;
				}
				buildPairInfo( r.Qname(), chrom, r.Position(), r.MatePosition(), r.MapQuality(), r.CigarString(), strand, alignPos, alignEnd, r.FirstFlag(), r.CountBWAChimericAlignments(), readPair_map);
				if (existsRR == false ) {
					r.SetPosition(r.Position());
					r.SetCigar(r.CigarString());
		 			RR.WriteRecord(r);
				}
				totalRR++;

				continue;
			}
			// FF read pair
			if (ori == 1) {

				if (r.MapQuality() < 20) { 
					continue;
				} 

				if ( isChrSomatic( chrom ) ) {
					readCounts_somatic.push_back(reads_window);
				}
				else {
					readCounts_germinal.push_back(reads_window); 
				}
				if (insert_size == 0) {
					continue;
				}
				bool isOnRegion = liesOnBlackRegion( mapOfRegions[chrom], chrom, readPos, readPosEnd);

				if (isOnRegion) {
					continue;
				}
				buildPairInfo( r.Qname(), chrom, r.Position(), r.MatePosition(), r.MapQuality(), r.CigarString(), strand, alignPos, alignEnd, r.FirstFlag(), r.CountBWAChimericAlignments(), readPair_map);
				if (existsFF == false ) {
					r.SetPosition(r.Position());
					r.SetCigar(r.CigarString());
	 				FF.WriteRecord(r);
				}
				totalFF++;

				continue;
			}
			// RF read pair
			if (ori == 2) {

				if (r.MapQuality() < 20) { 
					continue;
				} 

				if ( isChrSomatic( chrom ) ) {
					readCounts_somatic.push_back(reads_window);
				}
				else {
					readCounts_germinal.push_back(reads_window); 
				}
				if (insert_size >= upper_limit) {
					if (insert_size == 0) {
						continue;
					}
					bool isOnRegion = liesOnBlackRegion( mapOfRegions[chrom], chrom, readPos, readPosEnd);

					if (isOnRegion) {
						continue;
					}
					buildPairInfo( r.Qname(), chrom, r.Position(), r.MatePosition(), r.MapQuality(), r.CigarString(), strand, alignPos, alignEnd, r.FirstFlag(), r.CountBWAChimericAlignments(), readPair_map);
					if (existsRF == false ) {
						r.SetPosition(r.Position());
						r.SetCigar(r.CigarString());
	 					RF.WriteRecord(r);
					}
					totalRF++;
				}
				continue;
			}
		}
	}
}

void buildPairInfo ( const std::string& Qname, const std::string& chr, const int& start, const int& mateStart, const int& mapq, const std::string& cigar,  const std::string& strand, const int& align_pos, const int& align_end, const bool& is_first, const int& chimericAlignments,  std::map<std::string, ReadPairInfo>& readPair_map) {


	if ( chimericAlignments > 0 ) {

		// if we have the mate info
		auto it1 = readPair_map.find(Qname);
		if(it1 != readPair_map.end()) {

			if (is_first == true) {
				it1->second.first_mate_chr = chr;
				it1->second.first_mate_start = start;
				it1->second.first_mate_mapq  = mapq;
				it1->second.first_mate_cigar  = cigar;
				it1->second.first_mate_alignPos = align_pos;
				it1->second.first_mate_alignEnd = align_end;
				it1->second.first_mate_strand = strand;
			}
			else {
				it1->second.second_mate_chr =  chr;
				it1->second.second_mate_start = start;
				it1->second.second_mate_mapq  = mapq;
				it1->second.second_mate_cigar  = cigar;
				it1->second.second_mate_alignPos = align_pos;
				it1->second.second_mate_alignEnd = align_end; 
				it1->second.second_mate_strand = strand; 
			}
		}
		// if not read name inside
		else {

			if (is_first == true) {
				ReadPairInfo r1;
				r1.is_chimeric = false;
				r1.first_mate_chr =  chr;
				r1.first_mate_start = start;
				r1.first_mate_mapq  = mapq;
				r1.first_mate_cigar  = cigar;
				r1.first_mate_alignPos = align_pos;
				r1.first_mate_alignEnd = align_end;
				r1.first_mate_strand = strand;  

				readPair_map.insert(std::pair<std::string, ReadPairInfo>(Qname, r1));
			}
			else {
				ReadPairInfo r2;
				r2.second_mate_chr = chr;
				r2.second_mate_start = start;
				r2.second_mate_mapq  = mapq;
				r2.second_mate_cigar  = cigar;
				r2.second_mate_alignPos = align_pos;
				r2.second_mate_alignEnd = align_end;
				r2.second_mate_strand = strand;   
				readPair_map.insert(std::pair<std::string, ReadPairInfo>(Qname, r2));
			}
		}
	}
	else {
		if (is_first == true) {
			// Check if Read Id entry is present
			auto it1 = readPair_map.find(Qname);
			if(it1 != readPair_map.end()) {

				it1->second.first_mate_chr = chr;
				it1->second.first_mate_start = start;
				it1->second.first_mate_mapq  = mapq;
				it1->second.first_mate_cigar  = cigar;
				it1->second.first_mate_alignPos = align_pos;
				it1->second.first_mate_alignEnd = align_end; 
				it1->second.second_mate_strand = strand; 
			}
			else {
				ReadPairInfo r1;
				r1.is_chimeric = false;
				r1.first_mate_chr =  chr;
				r1.first_mate_start = start;
				r1.first_mate_mapq  = mapq;
				r1.first_mate_cigar  = cigar;
				r1.first_mate_alignPos = align_pos;
				r1.first_mate_alignEnd = align_end; 
				r1.first_mate_strand = strand;  

				readPair_map.insert(std::pair<std::string, ReadPairInfo>(Qname, r1));
			}
		}
		else {
			auto it2 = readPair_map.find(Qname);
			if(it2 != readPair_map.end()) {
				it2->second.second_mate_chr = chr;
				it2->second.second_mate_start = start;
				it2->second.second_mate_mapq  = mapq;
				it2->second.second_mate_cigar  = cigar;
				it2->second.second_mate_alignPos = align_pos;
				it2->second.second_mate_alignEnd = align_end; 
				it2 ->second.second_mate_strand = strand; 
			}
			else {
				ReadPairInfo r2;
				r2.second_mate_chr = chr;
				r2.second_mate_start = start;
				r2.second_mate_mapq  = mapq;
				r2.second_mate_cigar  = cigar;
				r2.second_mate_alignPos = align_pos;
				r2.second_mate_alignEnd = align_end; 
				r2.second_mate_strand = strand;   
				readPair_map.insert(std::pair<std::string, ReadPairInfo>(Qname, r2));
			}
		}
	}
}


