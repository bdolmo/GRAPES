#include <iostream>
#include <vector>
#include <string>
#include <stdlib.h>
#include <string.h>
#include "Trimmer.h"
#include <math.h>

 std::string Trimmer::trim() {

	// Simple trimming by Phred score threshold and a 4-bp sliding window
	// Todo: add universal adapter trimming
	// Trimming 3'
	int i = 0;
	for (i = qualities.length()-1; i>=0; i=i-4) {
		
		std::string codes_4;
		for (int j = i;j >= i-3; j--) {
			codes_4 += qualities[j];
		}
		int sum = 0;
		for (int i = 0; i < codes_4.length(); ++i) {
			// transforming to ASCII 33-offset
			int asci_decimal = codes_4[i]-33 >= 0 ? codes_4[i]-33 : 0;
			sum += asci_decimal;	
		}

		float asci_mean = sum/4;

		if ( i < 4) {	
			break;
		}

		if ( asci_mean > 15 ) {
			break;
		}
	}
	if ( i < 25 ) {
		i = 25;
	}

	trimmed = input_sequence.length()-i-1;
	std::string trimmedSeq  =  input_sequence.substr(0,i);
	std::string trimmedQual =  qualities.substr(0,i);
	//std::vector<std::string> output;

	//output.push_back(finalseq);
	//output.push_back(finalqual);
 
	return trimmedSeq;
 }

 int Trimmer::getTotalTrimmed() {
        return trimmed;
 }
