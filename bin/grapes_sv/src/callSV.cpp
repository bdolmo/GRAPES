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


using namespace std;
using namespace SeqLib;


bool callSV::resolve( BamRecordVector& sam, std::string chr, int refStart, int refEnd ) {

	for (auto& r : sam) {
		cout << r;

		int start = refStart + r.AlignmentPosition();
		int end   = refStart + r.AlignmentEndPosition();
		char strand    = r.ReverseFlag() == true ? '-' : '+';
      		if (!r.SecondaryFlag()) {
			cout << chr << "\t" << start << "\t" << end << "\t" << endl;
		}
	}
	return 1;
}






