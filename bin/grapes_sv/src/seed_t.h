#ifndef SEED_T_H
#define SEED_T_H

#include <iostream>
#include <string>

struct seed_t {
  int start;
  int end;
  int length;
  //std::string seq;
  int refStart;
  int refEnd;

   bool operator<(const seed_t& x) const {
    return (start < x.start);
   }

};

#endif
