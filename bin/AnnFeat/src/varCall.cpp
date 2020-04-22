#include <iostream>
#include <vector>
#include <string>
#include <stdlib.h>
#include <string.h>
#include "varCall.h"

#include <boost/algorithm/string/replace.hpp>
#include <boost/algorithm/string.hpp>

#include "SeqLib/BamReader.h"
#include "SeqLib/RefGenome.h"
#include "SeqLib/BamRecord.h"
#include "SeqLib/BFC.h"


using namespace SeqLib;

 struct varInfo {

  std::string chr;
  int position;
  int readSupport;
  int baseQ;
  char Ref;

  int Adenine = 0;
  int Guanine = 0;
  int Cytosine= 0;
  int Thymine = 0;	

 };

 int varCall::getLohVars() {
	return lohVars;
 }
 float varCall::getHomRatio() {
	return homRatio;
 }

 void varCall::callSNV( std::string chr, int start, int end, std::string svType ) {

 	std::string coordinate = chr + ":" +  std::to_string(start) + "-" + std::to_string(end);

	BamReader br;
	br.Open(bamFile);
	BamRecord r;

	std::string tag = "MD";

	GenomicRegion gr(coordinate, br.Header());
	br.SetRegion(gr);

	std::map<int,int> covBase;

	std::map <int, varInfo> varHash;

	while (br.GetNextRecord(r)) {
		if (r.MapQuality() < minMapQ) { 
			continue;
		} 
		if (r.CountNBases() > 0) {
			continue;
		}
		if (r.MaxDeletionBases() > 0 || r.MaxInsertionBases() > 0) {
			continue;
		}
		if (r.CountBWAChimericAlignments() > 0) { 
			continue; 
		}
		if (r.CountBWASecondaryAlignments() > 0) { 
			continue; 
		}

		// Skipping reads with more than 2 mismatches
		int numMismatches = r.Length() - r.NumClip() -r.NumMatchBases();
		if (numMismatches > 2) {
			continue;
		}
		int readLen = r.Length();
		int j = 0;
		for (int pos = r.Position() + r.AlignmentPosition()+1; pos <= r.Position() + r.AlignmentEndPosition(); pos++ ) {

			//std::string seq = r.Sequence();
			
			std::string qual = r.Qualities();
			int baseQuality = qual[j]-33;

			//Avoid calling SNVs at read edges
			if (j >10 && j < readLen-10) {
				if ( covBase.count(pos) == 0 ) {
					if (baseQuality > minBaseQ) {
						covBase.insert(std::pair<int, int>(pos, 1));
					}					
				} 
				else {
					if (baseQuality > minBaseQ) {				
						covBase[pos]++;
					}
				}
			}
			j++;
		}
				
		std::string tag = "MD";

		std::vector<std::string> MDtag = r.GetSmartStringTag(tag);

		for (auto&element : MDtag) {

			std::string tmp = "";
			int readPos = 0;

			for (int i=0; i<element.length(); i++) {
				if (element[i] == 'A' || element[i] == 'C' || element[i] == 'T' || element[i] == 'G') {
					if (tmp == "") {
						tmp = '0';
					}
					int tmpInt = std::stoi(tmp);
					readPos+= std::stoi(tmp);

					std::string seq = r.Sequence();
					int varPos = r.Position() + readPos+1;
					
					std::string qual = r.Qualities();

					char Alt = seq[readPos];
					char Ref = element[i];

					int baseQuality = qual[readPos]-33;
					
					if ( varHash.count(varPos) == 0) {

						varInfo var;
						var.chr = r.ChrID();
						var.position = varPos;
						var.readSupport = 1;
						var.Ref = Ref;
						if (baseQuality > minBaseQ) {
							if (Alt == 'A') {
								var.Adenine++;
							}
							else if (Alt == 'C') {
								var.Cytosine++;
							}
							else if (Alt == 'T') {
								var.Thymine++;
							}
							else if (Alt == 'G') {
								var.Guanine++;
							}
							varHash.insert(std::pair<int, varInfo>(varPos, var));
						}					
					} 
					else {
						if (baseQuality > minBaseQ) {
							varHash[varPos].readSupport++;
							if (Alt == 'A') {
								varHash[varPos].Adenine++;
							}
							else if (Alt == 'C') {
								varHash[varPos].Cytosine++;
							}
							else if (Alt == 'T') {
								varHash[varPos].Thymine++;
							}
							else if (Alt == 'G') {
								varHash[varPos].Guanine++;
							}
						}
					}
					tmp = "";
					readPos++;
				}
				else {
					tmp = tmp + element[i];
				}
			}
		}
	}
	
	// Greedy genotyping
	int totalVars = 0;
	int totalHom  = 0;
	int totalHet  = 0;

	for (auto const& var : varHash) {
		
		if (covBase[var.first] < minCov) {
			continue;
		}
		int max = 0;
		std::string alternative;

		if (varHash[var.first].Adenine > max) {
			max = varHash[var.first].Adenine;
			alternative = 'A';
		}
		if (varHash[var.first].Guanine > max) {
			max = varHash[var.first].Guanine;
			alternative = 'G';
		}		
		if (varHash[var.first].Thymine > max) {
			max = varHash[var.first].Thymine;
			alternative = 'T';
		}
		if (varHash[var.first].Cytosine > max) {
			max = varHash[var.first].Cytosine;
			alternative = 'C';
		}
	
		float alleleFreq = (double) max/covBase[var.first];	
		std::string genotype;

		if (alleleFreq >= 0.2 && covBase[var.first] > minCov && varHash[var.first].position >= start && varHash[var.first].position <= end  ) {
			totalVars++;

			if (alleleFreq > minHomRatio) {
				genotype = "1/1";
				totalHom++;
			}
			else {
				genotype = "0/1";
				totalHet++;
			}
			//std::cout << chr << "\t" << varHash[var.first].position << "\t" << alternative << "\t" << alleleFreq << "\t" << covBase[var.first] << std::endl;
		}
	}

	homRatio = totalVars > 0 ? (double)totalHom/totalVars : 0.00;
	lohVars = totalVars;

	//if (totalVars > 2 && homRatio < minHomRatio) {
	//	return 1;
	//} 
	//else if ( totalVars > 2 && homRatio >= minHomRatio) {
	//	return 2;
	//}	
	//else if ( totalVars <= 2) {
	//	return 3;
	//}
	// return 1: if more than 2 variants found and LOH does not support the deletion call. LOH filter
	// return 2: if more than 2 variants found and LOH is present. HQ deletion
	// return 3: if less than 2 variants found and LOH filter cannot be applied
}

