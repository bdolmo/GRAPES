#include <iostream>
#include <string.h>
#include <string>


class Trimmer {

	public:
		// Constructor
		Trimmer( std::string, std::string );

		// Read trimming by base quality
		std::string trim();

		int getTotalTrimmed();

	private:
	    std::string input_sequence;
	    std::string qualities;
	    int trimmed;
};


