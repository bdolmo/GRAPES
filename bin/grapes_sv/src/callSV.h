#include <iostream>
#include <string> // carreguem la llibreria d'strings
#include <fstream> // llibreria per obrir inputs
#include <vector>
#include <list>
#include <algorithm>
#include <omp.h>

#include "GenomicRange.h"
#include "pe_adjacency.h"
#include "ReadPairInfo.h"
#include "sam_t.h"
#include "vcf_t.h"

using namespace SeqLib;

 class callSV {

	// Mètodes
	public:
		// Constructor
		callSV(std::string&, std::string&);

		bool resolve(BamRecordVector&, std::string, int , int );


	private:
		std::string bamFile;
		std::string reference;		
 };
