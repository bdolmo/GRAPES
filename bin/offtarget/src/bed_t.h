#ifndef BED_T_H
#define BED_T_H

#include <string>
#include <iostream>

 struct bed_t  {

	std::string chr;
	int start;
	int end;
	std::string info;
  std::vector<std::string> fields;
	float gc;

	bool operator<(const bed_t& x) const {
		return (chr < x.chr);
	}
	inline bool comp(const bed_t& lhs, const bed_t& rhs){
	  return std::tie(lhs.chr, lhs.start, lhs.end, lhs.info, lhs.gc, lhs.fields) < std::tie(rhs.chr, rhs.start, rhs.end, rhs.info, rhs.gc, rhs.fields);
	}
 };

#endif
