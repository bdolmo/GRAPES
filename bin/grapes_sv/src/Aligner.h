#include <iostream>
#include <string.h>
#include <string> // carreguem la llibreria d'strings
#include <fstream> // llibreria per obrir inputs
#include <map>
#include <vector>

class Aligner {

	public:
		// Constructor
		Aligner(const std::string&, const std::string&, int);

		// Standard smith waterman alignment
		std::vector<int> localAligner();

		// getters
		int getScore();

		int getMismatches();
		
		int getGaps();		

	private:
		std::string query;
		std::string reference;
		int max_mismatches;
		int max_gaps;
		int mismatches;
		int gaps;
		int score;
};
