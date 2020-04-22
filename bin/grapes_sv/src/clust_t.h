#ifndef CLUST_T_H
#define CLUST_T_H


#include <iostream>
#include <string>

std::string IntToString ( int& );

 struct clust_t {

	int32_t chr;
	int32_t pos;
	bool operator <(const clust_t& x)const {
	    return std::tie(chr, pos) < std::tie(x.chr, x.pos);
	}

		
 };


#endif
