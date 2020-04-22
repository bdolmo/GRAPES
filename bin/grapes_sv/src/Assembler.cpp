#include <iostream>
#include <string> // carreguem la llibreria d'strings
#include <fstream> // llibreria per obrir inputs
#include <vector>
#include <list>
#include <algorithm>
#include <map>
#include "Assembler.h"
#include "ssw_cpp.h"
#include <unordered_map>

using namespace std;

 struct node_t {
     	std::string lab;
        int indegree  = 0;
        int outdegree = 0;

	bool operator<(const node_t& x) const {
		return (lab < x.lab);
	}
 };

string inline revComp(string seq) {

	string rev_seq(seq);
	reverse(rev_seq.begin(),rev_seq.end());
	for (int i = 0; i < rev_seq.length(); i++) {
		switch (rev_seq[i]) 
		{
			case 'A': rev_seq[i] = 'T'; break;
			case 'T': rev_seq[i] = 'A'; break;
			case 'C': rev_seq[i] = 'G'; break;
			case 'G': rev_seq[i] = 'C'; break;
			case 'a': rev_seq[i] = 't'; break;
			case 't': rev_seq[i] = 'a'; break;
			case 'c': rev_seq[i] = 'g'; break;
			case 'g': rev_seq[i] = 'c'; break;
		}
	}
	return rev_seq; 
}

std::vector<Nucleotide> revertConsensus( vector<Nucleotide>& vec) {

	std::vector<Nucleotide> outVec;
	for (auto& b : vec) {
		Nucleotide ntd;
		ntd.A = b.T;
		ntd.C = b.G;
		ntd.G = b.C;
		ntd.T = b.A;
		ntd.totalBases = b.totalBases;
		outVec.push_back(ntd);
	}
	reverse(outVec.begin(),outVec.end());
	return outVec;
}

std::multimap<std::string, std::string> Assembler::createReadTable() {

	std::multimap<std::string, std::string> readHash;

	for (auto read : reads) {
		readHash.insert(std::pair<std::string, std::string>(read, read));
	}
	return readHash;
}


std::unordered_multimap<std::string, Contig> Assembler::createPrefixTable() {

	std::unordered_multimap<std::string, Contig> prefixHash;

	for (auto contig : ContigList) {

		// Equal than SSAKE, VCAKE, popuating the first 11 bases on a hash table with the value being the full sequence
		std::string prefix = contig.seq.substr(0, kSize);
		std::string revCompRead = revComp(contig.seq);

		Contig revCompContig = contig;	 // Create revcomp contig struct
		revCompContig.seq = revCompRead; // RevComp of the contig
		std::string revCompPrefix = revCompRead.substr(0, kSize); // Prefix of revcomp

		prefixHash.insert(std::pair<std::string, Contig>(prefix, contig));
		prefixHash.insert(std::pair<std::string, Contig>(revCompPrefix, contig));
	}
	return prefixHash;
}

std::pair<std::string, int> mostVotedBase(Nucleotide Base) {
	int max = 0;
	std::string ntd;

	if (Base.A > max) {
		max = Base.A;
		ntd = "A";
	}
	if (Base.C > max) {
		max = Base.C;
		ntd = "C";
	}
	if (Base.T > max) {
		max = Base.T;
		ntd = "T";
	}
	if (Base.G > max) {
		max = Base.G;
		ntd = "G";
	}
	//cout << endl;
	return std::make_pair(ntd, max);	
}

int Assembler::getNumAssembled() {
	return totalAssembled;
}
int Assembler::getTotalReads() {
	return totalReads;
}

std::vector<std::string> Assembler::ungappedGreedy() {

	for (auto& r: reads) {
		Contig c;
		c.seq = r;

		for (auto& base :r) {
			Nucleotide ntd;
			ntd.totalBases++;

			if (base == 'A') {
				ntd.A=0;
			}
			if (base == 'C') {
				ntd.C=0;
			}
			if (base == 'T') {
				ntd.T=0;
			}
			if (base == 'G') {
				ntd.G=0;
			}
			c.Consensus.push_back(ntd);
		}
		ContigList.push_back(c);
	}

	int notExtended = 0;
	std::map<std::string, int> seenMap;
	std::map<std::string, int> seenPair;
	std::multimap<std::string, std::string> readsHash = createReadTable();
	std::unordered_multimap<std::string, Contig> prefixTable = createPrefixTable();

	//int totalReads = 100;
	int totalReads = readsHash.size();
	//cout << "Total reads: " << totalReads << endl;
	int round = 0;

	while (1) {

		Contig read = ContigList[0];
		round++;
		seenPair[read.seq]++;
		if (debug) {
			cout << "Ronda" << round << "\t" << read.seq << "\t" << seenPair[read.seq] << endl;
		}
		if (seenPair[read.seq] > 1) {
			bool allSeen = 0;
			for (auto& i : ContigList) {
				if (seenPair[i.seq] > 1) {
					allSeen = 1;
				}
				else {
					allSeen = 0;
				}
			}
			if (allSeen) {
				break;
			}
		}	

		if (reads.size() <= 1 ) {
			break;
		}

		Contig contig;

		//cout << "3'Extension..."<< endl;
		// 3' Extension
		bool isRevComp = false;
	    std::tie(notExtended, contig) = Extend( prefixTable, readsHash, seenPair, read, notExtended, isRevComp );
		// 5' Extension (through reverse complement)
		//cout << contig.seq << endl;

		isRevComp = true;
	    std::tie(notExtended, contig) = Extend( prefixTable, readsHash, seenPair, contig, notExtended, isRevComp );
	}
	for (auto& contig : reads) {
		//cout << "Contig:" << contig << endl;
	}
	// Remove non-assembled reads
	if (reads.size() > 1) {
		for (auto& r : readsHash) {
			//cout << "Eliminem " << r.first << endl;
			reads.erase(std::remove(reads.begin(), reads.end(), r.first), reads.end());
		}
	}
	totalAssembled = totalReads-reads.size();
	//cout << "Total assembled: " << totalAssembled << endl;
	std::vector<std::string> outReads;
	for (auto& r:reads) {
		if (r.length() >= kSize) {
			outReads.push_back(r);
		}
	}

	return outReads;
}

std::pair<int, Contig> Assembler::Extend( std::unordered_multimap<std::string, Contig>& prefixTable, std::multimap<std::string, std::string>& readHash, std::map<std::string, int>& seenPair, Contig& read, int notExtended, bool isRevComp ) {

	std::vector<candidateRead> arrayHits;
	std::vector<Nucleotide> consensusVector;
	std::string contig;
	bool isFound = false;
	std::string SEQ = read.seq;

	Contig outContig;
	if (isRevComp) {
		//cout << "Anem a assemblar: " << read.seq << "\t" << isRevComp << endl;
	}
	else {
		//cout << "Anem a assemblar: " << read.seq << "\t" << isRevComp << endl;
	}
	if (SEQ.length() < kSize) {

		ContigList.erase(std::remove(ContigList.begin(), ContigList.end(), read), ContigList.end());
		reads.erase(std::remove(reads.begin(), reads.end(), SEQ), reads.end());
		reads.push_back(read.seq);
		ContigList.push_back(read);
		return std::make_pair(notExtended, outContig);
	}

	std::string revCompRead;
	if (isRevComp) {
		revCompRead = revComp(read.seq);
	}

	int min = 100000;
	int max = 0;
	int startOffset;
	int endOffset;
	int lastPos = read.seq.length();
	int span;
	if (isRevComp) {
		SEQ = revComp(read.seq);
	}

	for (int i = 0; i < SEQ.length()-kSize; i++) {

		std::string seed = SEQ.substr(i, kSize);
			
		std::pair <std::unordered_multimap<std::string, Contig>::iterator, std::unordered_multimap<std::string, Contig>::iterator> ret;
		ret = prefixTable.equal_range(seed);

		for (std::unordered_multimap<std::string, Contig>::iterator it=ret.first; it!=ret.second; ++it) {

			int matches = 0;
			int mismatches = 0;
			int m = 0;
			std::string SEQ2 = it->second.seq;

			if (isRevComp) {
				SEQ2 = revComp(SEQ2);
			}

			for (int n = i+kSize; n < SEQ.length(); n++) {
				if (m+kSize > SEQ2.length()-1) {
					break;
				}
				if (SEQ[n] == SEQ2[m+kSize]) {
					matches++;
				}
				else {
					mismatches++;
				}
				m++;
			}
			int mmAllowd = (m * maxMismatches)/100;
			//cout << "Mismatchos:" << mismatches << "\t" << "Allowd:" << mmAllowd <<  "\t" << it->second.seq << endl;

			if (mismatches <= mmAllowd ) {
				if (isRevComp) {
					//seenPair[revComp(SEQ2)]++;
				}
				else {
					//seenPair[SEQ2]++;
				}
				//if (seenPair[SEQ2] > 1 || SEQ2 == SEQ) {
				if (SEQ2 == SEQ) {					
					continue;
				}

				if (!isRevComp) {
					//cout << " L'afegim a la llista de hits" << "\t" << seed << "\t" << it->second.seq << endl;
				}
				else {
					//cout << "  L'afegim a la llista de hits" << "\t" << revComp(seed) << "\t" << it->second.seq << endl;
				}

				isFound = true;
				candidateRead cr;

				if (isRevComp) {
					cr.Seed = revComp(it->first);
					cr.Sequence = SEQ2;
					cr.ContigStruct   = it->second;
					cr.ContigStruct.seq = it->second.seq;
					std::vector<Nucleotide> tmpConsensus = it->second.Consensus;
					reverse(tmpConsensus.begin(),tmpConsensus.end());
					cr.Consensus = tmpConsensus;
				}
				else {				
					cr.Seed = it->first;
					cr.Sequence = SEQ2;
					cr.ContigStruct   = it->second;
					cr.Consensus =  it->second.Consensus;
				}
				cr.offset   = i;
				arrayHits.push_back(cr);
			}
		}
	}
	if (isFound==false) {

		notExtended++;
		contig = SEQ;

		if (isRevComp) {
			if (debug) {
				cout << "found false revcomp" << endl;
			}

			contig = revComp(contig);
			reads.erase(std::remove(reads.begin(), reads.end(), contig), reads.end());
			reads.erase(std::remove(reads.begin(), reads.end(), revComp(contig)), reads.end());

			ContigList.erase(std::remove(ContigList.begin(), ContigList.end(), read), ContigList.end());
			//cout << endl;
			std::string prefixRead = read.seq.substr(0, kSize);
			std::string refCompPrefixRead = revComp(prefixRead);

			prefixTable.erase(prefixRead);
			prefixTable.erase(revComp(prefixRead));

		}
		else {
			//reads.erase(std::remove(reads.begin(), reads.end(), contig), reads.end());
			//reads.erase(std::remove(reads.begin(), reads.end(), revComp(contig)), reads.end());
			if (debug) {
				cout << "found false normal" << endl;
			}
		}

		ContigList.erase(std::remove(ContigList.begin(), ContigList.end(), read), ContigList.end());
		reads.erase(std::remove(reads.begin(), reads.end(), contig), reads.end());

		// Do not remove the read itself if it cannot be extended, instead, append it at the end
		ContigList.push_back(read);
		reads.push_back(contig);
		std::string prefix = contig.substr(0, kSize);

		std::string revCompContig = revComp(contig);
		std::string prefRevComp = revCompContig.substr(0, kSize);

		//prefixTable.insert(std::pair<std::string, Contig>(prefix, read));
		//prefixTable.insert(std::pair<std::string, Contig>(prefRevComp, read));
		//seenPair[contig]++;

		//cout << "Tenim disponibles:\n";
		//for (auto& l : prefixTable) {
			//cout << l.first << endl;
		//}
		//cout << endl << endl;
		outContig = read;
	}
	else {

		// Adding read offsets and others
		candidateRead cr;
		cr.Sequence  = SEQ;

		if (isRevComp) {
			if (debug) {
				cout << "found true revcomp\n";
			}
			std::vector<Nucleotide> tmpConsensus = read.Consensus;
			reverse(tmpConsensus.begin(),tmpConsensus.end());
			cr.Consensus = tmpConsensus;
		}
		else {
			if (debug) {
				cout << "found true normal\n";
			}
			cr.Consensus = read.Consensus;
		}
		cr.offset    = 0;
		cr.Seed      = read.seq.substr(0, kSize);
		arrayHits.push_back(cr);

		for (auto& h : arrayHits) {

			if (h.offset < min) {
				startOffset = h.offset;
				min = startOffset;
			}
			if (h.offset+ h.Sequence.length() > max) {
				endOffset = h.offset+ h.Sequence.length();
				max = endOffset;
			}
			std::string revCompRead = revComp(h.Sequence);

			if (h.Sequence != SEQ) {
				readHash.erase(h.Sequence);
				readHash.erase(revCompRead);
			}

			std::string prefixInitRead = h.Sequence.substr(0, kSize);
			std::string revCompPrefix = revCompRead.substr(0, kSize);

			//prefixTable.erase(prefixInitRead);
			//prefixTable.erase(revCompPrefix);	
			//cout << "buscant " << h.Sequence.substr(0, kSize) << " " << revCompPrefix << endl;

			std::pair <std::unordered_multimap<std::string, Contig>::iterator, std::unordered_multimap<std::string, Contig>::iterator> ret2;
			ret2 = prefixTable.equal_range(prefixInitRead);

			for (std::unordered_multimap<std::string, Contig>::iterator it=ret2.first; it!=ret2.second; ++it) {
				if (h.Sequence == it->second.seq or h.Sequence == revComp(it->second.seq)) {
					prefixTable.erase(it);
					//cout << " 1Prefixhash removing " << prefixInitRead << "\t" << it->second.seq << endl;
					break;
				}
			}

			std::pair <std::unordered_multimap<std::string, Contig>::iterator, std::unordered_multimap<std::string, Contig>::iterator> ret3;
			ret3 = prefixTable.equal_range(revCompPrefix);
			//cout << "buscant " << revCompPrefix << endl;

			for (std::unordered_multimap<std::string, Contig>::iterator it=ret3.first; it!=ret3.second; ++it) {
				if (h.Sequence == it->second.seq or h.Sequence == revComp(it->second.seq)) {
					prefixTable.erase(it);		
					//cout << " 2Prefixhash removing " << revCompPrefix << "\t" << it->second.seq << endl;
					break;
				}
			}

			reads.erase(std::remove(reads.begin(), reads.end(), h.Sequence), reads.end());
			reads.erase(std::remove(reads.begin(), reads.end(), revCompRead), reads.end());
			//cout << " removing " << read.seq <<" "<< h.ContigStruct.seq << endl;

			ContigList.erase(std::remove(ContigList.begin(), ContigList.end(), read), ContigList.end());
			ContigList.erase(std::remove(ContigList.begin(), ContigList.end(), h.ContigStruct), ContigList.end());
		}

		// Initializing consesnsusVector
		for (int pos = 0; pos < endOffset; pos++) {
			Nucleotide base;
			consensusVector.push_back(base);
		}	
		for(auto& r : arrayHits) {
			Nucleotide base;
			int j = r.offset;

			int i = 0;
			//cout << "iterem per r.Sequence " << r.Sequence << endl;

  			for (auto& b : r.Sequence){ 
				//cout << b << " [A:" << r.Consensus[i].A << ",C:" << r.Consensus[i].C << ",T:" << r.Consensus[i].T << ",G:" << r.Consensus[i].G << "] ";

				consensusVector[j].totalBases+= r.Consensus[i].totalBases;

				if (isRevComp) {
					consensusVector[j].A += r.Consensus[i].T;
					consensusVector[j].C += r.Consensus[i].G;
					consensusVector[j].T += r.Consensus[i].A;
					consensusVector[j].G += r.Consensus[i].C;
				}	
				else {
					consensusVector[j].A += r.Consensus[i].A;
					consensusVector[j].C += r.Consensus[i].C;
					consensusVector[j].T += r.Consensus[i].T;
					consensusVector[j].G += r.Consensus[i].G;
				}
				switch(b) { 		
					case 'A':
						consensusVector[j].A++;
						break;
					case 'C':
						consensusVector[j].C++;
						break;		
					case 'T':
						consensusVector[j].T++;
						break;
					case 'G':
						consensusVector[j].G++;
						break;															
				}
				j++;
				i++;
			}
			//cout << endl;
		}
		int x = 0;
		std::vector<Nucleotide> NewConsensus;

		for (auto& ntd : consensusVector) {

			std::string b;
			int totalb;

			if (ntd.totalBases > 1) {
				std::tie(b, totalb) = mostVotedBase(ntd);
				contig+=b;
				NewConsensus.push_back(ntd);
			}
		}
		if (isRevComp) {
			contig = revComp(contig);
		}

		ContigList.erase(std::remove(ContigList.begin(), ContigList.end(), read), ContigList.end());

		readHash.erase(read.seq);
		reads.erase(std::remove(reads.begin(), reads.end(), read.seq), reads.end());

		reads.erase(std::remove(reads.begin(), reads.end(), revCompRead), reads.end());
		readHash.erase(revCompRead);

		auto it = reads.insert(reads.begin(), contig);
		std::string prefixContig = contig.substr(0, kSize);

		Contig NewContig;
		NewContig.seq = contig;
		NewContig.Consensus = NewConsensus;

		if (debug) {
			cout << "Contig => " << contig << endl;
		}

		if (isRevComp) {
			NewContig.Consensus = revertConsensus(NewConsensus);
		}
		if (contig.length() >= kSize) {
			if (debug) {
				cout << "size is ok " << contig.length() << "\t" << kSize << "\t" << reads.size() << endl;
			}
			ContigList.insert(ContigList.begin(), NewContig);
			for(auto& r:ContigList) {
				//cout << "disponible: " << r.seq << endl;
			}
			for(auto& r:prefixTable) {
				//cout << "PrefixHash disponible: " << r.first << " " << r.second.seq << endl;
			}

			//cout << endl;
			std::string prefixRead = read.seq.substr(0, kSize);
			std::string refCompPrefixRead = revComp(prefixRead);

			prefixTable.erase(prefixRead);
			prefixTable.erase(revComp(prefixRead));

			prefixTable.insert(std::pair<std::string, Contig>(prefixContig, NewContig));
			prefixTable.insert(std::pair<std::string, Contig>(revComp(prefixContig), NewContig));

			outContig = NewContig;
		}
		else {
			if (debug) {
				cout << "size is notok " << contig.length() << endl;
			}
			std::string prefixRead        = read.seq.substr(0, kSize);
			std::string refCompPrefixRead = revComp(prefixRead);

			prefixTable.erase(prefixRead);
			prefixTable.erase(revComp(prefixRead));

			reads.erase(std::remove(reads.begin(), reads.end(), read.seq), reads.end());
			reads.erase(std::remove(reads.begin(), reads.end(), revCompRead), reads.end());

			readHash.erase(revComp(contig));
			readHash.erase(contig);

			/*prefixTable.erase(contig);
			prefixTable.erase(revComp(contig));
			
			ContigList.erase(std::remove(ContigList.begin(), ContigList.end(), contig), ContigList.end());
			ContigList.erase(std::remove(ContigList.begin(), ContigList.end(), read), ContigList.end());

			reads.erase(std::remove(reads.begin(), reads.end(), read.seq), reads.end());
			reads.erase(std::remove(reads.begin(), reads.end(), revCompRead), reads.end());
			readHash.erase(revComp(contig));
			readHash.erase(contig);*/
		}

	}
	return std::make_pair(notExtended, outContig);
}

int getScore (const StripedSmithWaterman::Alignment& alignment){
	return alignment.sw_score;
}
int getRefBegin (const StripedSmithWaterman::Alignment& alignment){
	return alignment.ref_begin;
}

int getRefEnd (const StripedSmithWaterman::Alignment& alignment){
	return alignment.ref_end;
}

int getQueryBegin (const StripedSmithWaterman::Alignment& alignment){
	return alignment.query_begin;
}

int getQueryEnd (const StripedSmithWaterman::Alignment& alignment){
	return alignment.query_end;
}
int getMismatches (const StripedSmithWaterman::Alignment& alignment){
	return alignment.mismatches;
}

std::string getCigarString (const StripedSmithWaterman::Alignment& alignment){
	return alignment.cigar_string;
}

std::vector<uint32_t> getCigar (const StripedSmithWaterman::Alignment& alignment){
	return alignment.cigar;
}

std::pair<int, int> cigar2Gaps (std::string& cigar ) {

	std::string concat;
	int gaps = 0;
	int matches = 0;
	for (auto& i : cigar)  {
		if (i == 'D' or i == 'I' or i == 'M' or i == 'X' or i == 'S' or i == '=') {
			if (i == 'D' or i == 'I') {
				int n = std::stoi(concat);
				gaps +=n;	
			}
			if (i == '=') {
				int n = std::stoi(concat);
				matches +=n;	
			}
			concat = "";				
		}
		else  {
			concat += i;
		}
	}
	return std::make_pair(gaps, matches);
}


 std::string classifyOverlap(int query_begin, int query_end, int ref_begin, int ref_end, int l_query, int l_ref, std::string query_read, std::string ref_read) {

	std::string contig;

	//Case1        Prefix(query) - Suffix(Ref)
        // query         T A A C T C T A G C A C A T T T A C A
	//	         | | | | | | | | | | | | | | | |	
   	// ref   G A C C T A A C T C T A G C A C A T T T
	if (ref_begin >= query_begin && ref_end > query_end ) {

		std::string a = query_read.substr(query_end+1);
		contig = ref_read + a;
		return contig;	
	}
	//Case2        Suffix(query) - Prefix(Ref)
   	// query   G A C C T A A C T C T A G C A C A T T T
	//	           | | | | | | | | | | | | | | | |	
        // ref             T A A C T C T A G C A C A T T T A C A
	else if (query_begin >= ref_begin && query_end > ref_end ) {
     	std::string a = ref_read.substr(ref_end+1);
		contig = query_read + a;
		return contig;
	}
	//Case3 if whole query inside ref and ref is larger than query on both ends
        // query        T G A A T G C C C T T T A C T T T T A A
        //		| | | | | | | | | | | | | | | | | | | |        
        // ref      T T T G A A T G C C C T T T A C T T T T A A G T 

	else if (query_begin <= ref_begin && query_end <= ref_end ) {
		contig = ref_read;
		return contig;
	}
	//Case4 if whole ref inside query and ref is larger than query on both ends
        // query      T T T G A A T G C C C T T T A C T T T T A G      
	//		  | | | | | | | | | | | | | | | | | |
        // ref            T G A A T G C C C T T T A C T T T T 
	else if (ref_begin <= query_begin &&  ref_end <= query_end ) {
		contig = query_read;
		return contig;
	}
	//Case5  No accepted overlaps
	else {
		contig = "N";
		return contig;
	}
 }

overlapData Assembler::findBestOverlap() {

	std::vector<std::string> reads2 = reads;

	int bestScore = 0;
	std::string bestPair;	
	std::string contig = "N";

	int matchesOut;
	int mismatchesOut;
	int gapsOut;
	std::string refOut;
	std::string queryOut;
	bool isRevComp = false;
	for (auto& i : reads ) {

    	reads2.erase (reads2.begin());
		for (auto& j : reads2) {

			// Declares a default Aligner
			StripedSmithWaterman::Aligner ssw_aligner;

			// Declares a default filter
			StripedSmithWaterman::Filter filter;

			// Declares an alignment that stores the result
			StripedSmithWaterman::Alignment alignment;
			ssw_aligner.Align(j.c_str(), i.c_str(), i.length(), filter, &alignment);

			std::string jrevComp = revComp(j);
			
			StripedSmithWaterman::Alignment alignmentRev;
			ssw_aligner.Align(jrevComp.c_str(), i.c_str(), i.length(), filter, &alignmentRev);

			int score, mismatches, matches, gaps, refBegin, refEnd, queryBegin, queryEnd;

			if ( getScore(alignment) >= getScore(alignmentRev) ) {
				score      = getScore(alignment);
				mismatches = getMismatches(alignment);
				refBegin   = getRefBegin(alignment);
				refEnd     = getRefEnd(alignment);
				queryBegin = getQueryBegin(alignment);
				queryEnd   = getQueryEnd(alignment);
				std::string cigarStr = getCigarString(alignment);
	  			std::tie(gaps, matches) = cigar2Gaps(cigarStr);
			}
			else {
				score      = getScore(alignmentRev);
				mismatches = getMismatches(alignmentRev);
				refBegin   = getRefBegin(alignmentRev);
				refEnd     = getRefEnd(alignmentRev);
				queryBegin = getQueryBegin(alignmentRev);
				queryEnd   = getQueryEnd(alignmentRev);
				std::string cigarStr = getCigarString(alignmentRev);
	  			std::tie(gaps, matches) = cigar2Gaps(cigarStr);
				isRevComp = true;
			}

			if (score >= bestScore && matches >= minOlapLength && gaps < maxGaps && mismatches < maxMismatches) {
				bestScore = score;
				refOut = i;
				queryOut = j;
				matchesOut = matches;
				mismatchesOut = mismatches;
				gapsOut = gaps;
 				contig = classifyOverlap(queryBegin, queryEnd, refBegin, refEnd, j.length(), i.length(), j, i);
			}
		}
	}
	
	overlapData Contig;

	Contig.seq1      = refOut;
	Contig.seq2      = queryOut;
	Contig.matches   = matchesOut;
	Contig.mismatches = mismatchesOut;
	Contig.gaps      = gapsOut;
	Contig.contigSeq = contig;
	Contig.isRevComp = isRevComp;

	return Contig;
}

std::vector<std::string> Assembler::gappedGreedy() {

	// Greedy assembly consists on two simple steps:
	// Compute all overlaps between all reads (aka build Overlap graph)
	// Pick best overlapping pairs and build a contig
	// iterate process until no reads can be incorporated to the contig
	
	if (reads.size() == 1) {
		return reads;
	}

	overlapData contigData = Assembler::findBestOverlap();
	
	int matches        = contigData.matches;
	int mismatches     = contigData.mismatches;
	int gaps           = contigData.gaps;
	std::string contig = contigData.contigSeq;
	std::string seq1   = contigData.seq1;
	std::string seq2   = contigData.seq2;
	bool isRevComp     = false;

	if (reads.size() == 2 && contig != "N") {
		std::vector<std::string> out;
		out.push_back(contig);
		totalAssembled+=2;
		return out;
	}
	if ( matches >= minOlapLength && mismatches < maxMismatches && gaps < maxGaps) {
		totalAssembled++;
	}
	int iteration = 0;

	while (matches >= minOlapLength && gaps < maxGaps && mismatches < maxMismatches && contig != "N") {

		reads.erase(std::remove(reads.begin(), reads.end(), seq1), reads.end());
		reads.erase(std::remove(reads.begin(), reads.end(), seq2), reads.end());
		reads.push_back(contig);

		if (reads.size() == 1) {
			return reads;
		}

		contigData = Assembler::findBestOverlap();
		contig    = contigData.contigSeq;
		seq1      = contigData.seq1;
		seq2      = contigData.seq2;
		matches   = contigData.matches;
		mismatches= contigData.mismatches;
		gaps      = contigData.gaps;
		isRevComp = contigData.isRevComp;
		iteration++;
	}
	return reads;
}

std::vector<std::string> Assembler::deBruijnGraph() {

	int k_size = 21;
	std::vector<std::string> tmp;

	std::multimap<node_t, node_t> edges;
	std::map<node_t, node_t> vertices;

	// 1. Chop string into k-mers
	for (int i = 0; i < reads.size(); i++) {

		for (int j = 0; j < reads[i].length()-k_size+1; j++) {
			std::string subseq = reads[i].substr(j, k_size);

			std::string v1 = subseq.substr(0, k_size-1);
			std::string v2 = subseq.substr(1, k_size);

           		node_t nodeL;
			nodeL.lab = v1;
			
			node_t nodeR;
			nodeR.lab = v2;

		    // Node not found, so we will append a new one
			if ( edges.count(nodeL) > 0 ) {
				vertices[nodeL].outdegree++;
				edges.insert(std::pair<node_t, node_t>(nodeL, nodeR ));
			} 
			else {	 	
				vertices[nodeL].outdegree++;
				edges.insert(std::pair<node_t, node_t>(nodeL, nodeR ));
			}
		    // Edge not found, so we will append a new one
			if ( edges.count(nodeR) > 0) {
				vertices[nodeR].indegree++;
			}
			else  {
				vertices[nodeR].indegree++;
				edges.erase(nodeR);
			}
		 }
	}

	// Pick starting node
	node_t First_node       = vertices.begin()->first;
	node_t First_indegree   = vertices.begin()->second;

	std::vector<node_t> initial_nodes;
	for (std::map<node_t, node_t>::iterator it=vertices.begin(); it!=vertices.end(); ++it) {

		if (it->second.indegree < First_indegree.indegree) {
			First_node = it->first;
			First_indegree.indegree = it->second.indegree;
		}
	}
	for (std::map<node_t, node_t>::iterator it=vertices.begin(); it!=vertices.end(); ++it) {
		if (it->second.indegree == First_node.indegree) {
			initial_nodes.push_back(it->first);
		}
	}

	std::multimap<node_t, node_t> old_edges;
	old_edges = edges;

	int contigL = 0;
	std::string longestContig;

	std::vector<std::string> contigs;

	for ( int i=0; i < initial_nodes.size(); i++) {

		First_node = initial_nodes[i];

		std::string contig = First_node.lab;
		node_t current = First_node;

		edges = old_edges;

		while (edges.count(First_node) > 0) {

			node_t next_node = edges.find(First_node)->second;
			contig += next_node.lab.back();

			std::multimap<node_t,node_t>::iterator it = edges.find(First_node);
			edges.erase(it);
			First_node = next_node; 
		}
		contigs.push_back(contig);
		if (contig.length() > contigL) {

			contigL= contig.length();
			longestContig = contig;
		}
        }
	return contigs;
 }


