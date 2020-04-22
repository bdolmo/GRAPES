#include <iostream>
#include <string.h>



class varCall {

	public:
		// Constructor
		varCall( std::string, std::string, int, int, int, int, float );

		// Determine LOH
		int isLOH( std::string, int, int, std::string );
		
	private:
		std::string bamFile;
		std::string genome;
		int minMapQ;
		int minBaseQ;
		int minSNV;
		int minCov;
		float minHomRatio;
};
