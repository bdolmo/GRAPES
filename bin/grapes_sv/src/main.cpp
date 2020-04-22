
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
#include "ssw_cpp.h"

#include "assembler.h"

using namespace std;

assembler::assembler( std::vector<std::string>& _reads, int _minOlapLength, int _maxMismatches, int _maxGaps, int _kSize) {
	reads = _reads;
	minOlapLength = _minOlapLength;
	maxMismatches = _maxMismatches;
	maxGaps  = _maxGaps;
	kSize    = _kSize;
}

//######################################### Calculate median
string inline fillQualString ( std::string& seq ) {

	string outQual;

	for (int i = 0; i < seq.length(); i++) {
		outQual+="I";
	}
	return outQual;
}


int main (int argc, char* argv[]) {
	if (argc < 1) {
		cout << "Usage: sequences.txt>" << endl;
		return 0;
	}
	string Input(argv[1]);
	ifstream inputFile;
  	inputFile.open (Input);

	std::vector<std::string> reads;

        if (inputFile.is_open()) {
		string line;
		while ( std::getline (inputFile, line)) {
			reads.push_back(line);
		}
	}
	assembler A(reads, 10, 3, 2, 20);
	std::vector<std::string> contigs = A.ungappedGreedy();

	assembler B(contigs, 10, 3, 2, 20);
	std::vector<std::string> reContigs = B.gappedGreedy();

	int n = 0;
	for (auto& i : reContigs) {
		cout << "@contig" << n << endl;
		cout <<  i << endl;
		cout << "+" << endl;
		cout << fillQualString(i) << endl;
		n++;
	}
    return 0;
}
