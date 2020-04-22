#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <algorithm>
#include <tuple>
#include <sstream>
#include <regex>
#include <fstream>
#include "seed_t.h"

 class hashAligner {

	// Mètodes
	public:
		// Constructor
		hashAligner(int, int, int);

		std::multimap<std::string, int> indexContig(std::string&);

		std::vector<seed_t> findSeeds(std::string&, std::string&, std::multimap<std::string, int>&);

		std::vector<seed_t> mergeIntervals( std::vector<seed_t>&);

		void createDAG(std::vector<seed_t>& intervalVector);


	private:
		int kSize;
		int minSeeds;
		int maxMismatches;

 };



