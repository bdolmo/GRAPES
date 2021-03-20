#include <iostream>
#include <string.h>
#include <string>
#include <vector>
#include <iostream>
#include <map>
#include <regex>
#include <string.h>
#include <numeric>
#include <algorithm>
#include <math.h>      
#include <omp.h>

#include "varCall.h"
#include "SeqLib/BamReader.h"
#include "SeqLib/BamHeader.h"
#include "SeqLib/RefGenome.h"
#include "SeqLib/BamRecord.h"
#include "SeqLib/BamWriter.h"
#include "Utils.cpp"
#include <boost/progress.hpp>


using namespace SeqLib;
using namespace std;

varCall::varCall( std::string _bamFile, std::string _genome, int _minMapQ, int _minBaseQ, int _minSNV, int _minCov, float _minHomRatio ) {
	bamFile    = _bamFile;
	genome     = _genome;
	minMapQ    = _minMapQ;
	minBaseQ   = _minBaseQ;
	minSNV     = _minSNV;
	minCov     = _minCov;
	minHomRatio= _minHomRatio;
}

double  calculateInsertSizeLimits (std::string bamFile) {

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
	int readLength;
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
			if (insert_size >= maxISize) { continue; }
			sum+=insert_size;
			isizes.push_back(insert_size);
		}
	}

	double mean =  sum / isizes.size();

	int distance = positions[positions.size()-1] - positions[0];

	double mean_coverage = ( (double)read_counts/(double)distance )*readLength;

	for (auto& i : isizes) {
    	accum += (i - mean) * (i - mean);
	}
	double stdev = sqrt(accum / (isizes.size()-1));
	return mean_coverage;
}


int main (int argc, char* argv[]) {
	if (argc < 3) {
		cout << "Usage: ADDINFO <BED> <BAM> <REFERENCE> <OUTBED> <CPU_CORES>" << endl;
		return 0;
	}
	string inBed(argv[1]);
	string bamFile(argv[2]);
	string reference(argv[3]);
	string outBed(argv[4]);
	string threads(argv[5]);

	int numCpus = std::stoi(threads);
	omp_set_num_threads(numCpus);

	RefGenome ref;
	ref.LoadIndex(reference);

	long genomeSize = 2966866909;

	vector<string> strV;
	double meanCoverage = calculateInsertSizeLimits(bamFile);

	ifstream bedFile;
  	bedFile.open (inBed);

        if (bedFile.is_open()) {
		string line;
		while ( std::getline (bedFile, line)) {
			strV.push_back(line);
		}
	}
	bedFile.close();

	/*cout << "#chr" << "\t";
	cout << "start" << "\t";
	cout << "end" << "\t";
	cout << "precision" << "\t";
	cout << "svtype" << "\t";
	cout << "fallsInCNVR" << "\t";
	cout << "svlen" << "\t";
	cout << "mapq" << "\t";
	cout << "kdiv" << "\t";
	cout << "breakreads" << "\t";
	cout << "assembled" << "\t";
	cout << "discordants" << "\t";
	cout << "ratio_RD" << "\t";
	cout << "MAD_RD" << "\t";
	cout << "phred_discordant" << "\t";
	cout << "phred_RD" << "\t";
	cout << "totalSNV" << "\t";
	cout << "Homozygous_ratio" << "\t";
	cout << "gc_content" << endl;*/

  	//boost::progress_display show_progress( strV.size() );

	#pragma omp parallel
	{
		#pragma omp for schedule(dynamic, 1)
		for (int i = 0; i < strV.size(); i++) {

				string line = strV[i];

				vector<string> tmp = split(line, '\t');
				string chr      = tmp[0];
				int start       = stoi(tmp[1]);
				int end         = stoi(tmp[2]);
				string precision= tmp[3];
				string svtype   = tmp[4];
				string fallsInCNVR = tmp[5];
				int svlen       = stoi(tmp[6]);
				int mapq        = stoi(tmp[7]);
				float kdiv      = stof(tmp[8]);
				int breakreads  = stoi(tmp[10]);
				int assembled   = stoi(tmp[11]);
				int discordants = stoi(tmp[12]);

				float ratioRD;
				float madRD;	
	
				if (svtype != "DEL") {
					continue;
				}
				if (tmp[13] != ".") {
					madRD = stof(tmp[14]);	
				}
				else {
					madRD = 0.10;
				}

				string LOHsupp = tmp[15];
				int cumulative = stoi(tmp[16]);
				int numInserts = stoi(tmp[17]);
				string evidence = tmp[18];
	
				double pvalue_upstream, pvalue_downstream, pvalue_twosided;
				int counts_5prime, counts_inner, counts_3prime;
				int coverage_5prime, coverage_inner, coverage_3prime;

				std::tie(pvalue_upstream, pvalue_downstream, pvalue_twosided, counts_5prime, counts_inner, counts_3prime,
				coverage_5prime, coverage_inner, coverage_3prime) = calculatePvalueCoverage(chr, start, end, svtype, precision, bamFile);
				//cout << coverage_5prime << "\t" << coverage_inner << "\t" << coverage_3prime << endl;

				// Checking if the variant is found on a possible segmental duplication
				double ratio5prime = (double)coverage_5prime/(double)meanCoverage;
				double ratio3prime = (double)coverage_3prime/(double)meanCoverage;

				int phred_coverage;
				if ( pvalue_twosided != 0) {
					phred_coverage = ((-10) * (log10(pvalue_twosided)));
				}
				else {
					phred_coverage = 99;
				}
			
				double AF_5prime = (double)counts_inner/counts_5prime;
				double AF_3prime = (double)counts_inner/counts_3prime;
				double AF = (AF_5prime+AF_3prime)/2;

				if (tmp[13] != ".") {
					ratioRD = stof(tmp[13]);	
				}
				else {
					ratioRD = AF;
				}

				double discPvalue = computeDiscordantClusterSignificance (cumulative,  genomeSize, numInserts, discordants);
				int phred_discordant = ((-10) * (log10(discPvalue)));
				if (discordants > 0) {
					if (discPvalue == 0 || !discPvalue) {
						phred_discordant = 99;
					}
					else {
						phred_discordant = phred_discordant < 99 ? ((-10) * (log10(discPvalue))) : 99;
					}
				}
				else {
					phred_discordant = 0;
				}

				varCall SNV(bamFile, reference, 50, 19, 5, 12, 0.8);

				int totalSNV = 0;
				float homRatio = 0.00;

				if (svtype == "DEL" && fallsInCNVR == "no" && svlen > 300) {

					int tmpStart = start;
					int tmpEnd   = end;

					// Setting a padding distance on imprecise calls to avoid calling false positive SNV on breakpoint boundaries
					if (start+100 < end -100) {
						tmpStart = start+100;
						tmpEnd   = end-100;
					}

					SNV.callSNV(chr, tmpStart, tmpEnd, svtype);
					totalSNV = SNV.getLohVars();
					homRatio = SNV.getHomRatio();
				}
				float gc_content = 0.00;
				gc_content = computeGC( chr, start, end, reference);
				#pragma omp critical
				{
					//++show_progress;
					cout << chr << "\t";
					cout << start << "\t";
					cout << end << "\t";
					cout << precision << "\t";
					cout << svtype << "\t";
					cout << fallsInCNVR << "\t";
					cout << svlen << "\t";
					cout << mapq << "\t";
					cout << kdiv << "\t";
					cout << breakreads << "\t";
					cout << assembled << "\t";
					cout << discordants << "\t";
					cout << ratioRD << "\t";
					cout << madRD << "\t";
					cout << phred_discordant << "\t";
					cout << phred_coverage << "\t";
					cout << totalSNV << "\t";
					cout << homRatio << "\t";
					cout << gc_content << "\t";
					cout << ratio5prime << "\t";
					cout << ratio3prime << "\t";
					cout << evidence << endl;
				}
		}
	}

return 0;
}
