#ifndef PE_ADJACENCY_H
#define PE_ADJACENCY_H


#include <iostream>
#include <string>

struct PE_adjacency {

	std::string chromosome;
	int start;
	int end;

	std::vector<int> align_begins;
	std::vector<int> align_ends;
	std::vector<std::string> first_mate_strands;
	std::vector<std::string> second_mate_strands;

	int supporting = 0;
	bool operator<(const PE_adjacency& x) const {
		return (start < x.start);
	}
};

#endif
