#include <iostream>
#include <string.h>

class varCall {

	public:
		// Constructor
		varCall( std::string, std::string, int, int, int, int, float );

		// Determine LOH
		void callSNV( std::string, int, int, std::string );

		int getLohVars();
		float getHomRatio();
		
	private:
		std::string bamFile;
		std::string genome;
		int minMapQ;
		int minBaseQ;
		int minSNV;
		int minCov;
		float minHomRatio;
		int lohVars = 0;
		float homRatio = 0.00;

};
