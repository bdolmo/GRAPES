#include <iostream>
#include <vector>
#include <string> 

#include "Aligner.h"

 int Aligner::getScore() {
	return score;
 }

 int Aligner::getMismatches() {
	return mismatches;
 }

 int Aligner::getGaps() {
	return gaps;
 }

 double mu    = 1.33;
 double delta = 1.33;
 int ind;

 ////////////////////////////////// Auxiliary functions /////////////////////////////////

 double similarity_score(char a,char b){

	  double result;
	  if(a==b){
	      result=1.;
	  }
	  else{
	      result=-mu;
	  }
	  return result;
 }

 /////////////////////////////////////////////////////////////////////////////

 double find_array_max(double array[],int length) {

	  double max = array[0];            // start with max = first element
	  int ind = 0;
	  for(int i = 1; i<length; i++){
	      if(array[i] > max){
		max = array[i];
		ind = i; 
	      }
	  }
	  return max;                    // return highest value in array
 }


 std::vector<int> Aligner::localAligner() {


 mismatches = 0;
 gaps = 0;

 std::vector<int> output;

 int N_max = 1000;
 std::string header;
	
 std::string seq_a,seq_b; 
	  
 seq_a = query;  
 seq_b = reference;

 int N_a = seq_a.length();                     // get the actual lengths of the sequences
 int N_b = seq_b.length();
	 
 // initialize H
 double H[N_a+1][N_b+1];     
  for(int i=0;i<=N_a;i++){
    for(int j=0;j<=N_b;j++){
      H[i][j]=0.;
    }
  } 

  double temp[4];
  int I_i[N_a+1][N_b+1],I_j[N_a+1][N_b+1];     // Index matrices to remember the 'path' for backtracking
	  
  // here comes the actual algorithm

  int MATCH = 1;
  int MISMATCH = -1;
  int GAP = -2;

  for(int i=1;i<=N_a;i++){
        for(int j=1;j<=N_b;j++){
	       if ( seq_a[i-1] == seq_b[j-1]) {
	      		temp[0] = H[i-1][j-1] + MATCH;
	       }
	       else {
			temp[0] = H[i-1][j-1] + MISMATCH;
	       }
	       temp[1] = H[i-1][j]+GAP;                  
	       temp[2] = H[i][j-1]+GAP;                 
	       temp[3] = 0.;
	       double max = temp[0];
	       int ind = 0;
	       for(int i = 1; i<4; i++){
		    if(temp[i] > max){
			max = temp[i];
			ind = i; 
		     }
	       }
	       H[i][j] = max;	
	       if (ind == 0) {                                  // score in (i,j) stems from a match/mismatch
	   		I_i[i][j] = i-1;
			I_j[i][j] = j-1;
	       }		
	       if (ind == 1) {                                  // score in (i,j) stems from a deletion in sequence A
	     		I_i[i][j] = i-1;
			I_j[i][j] = j;
	       }		
	       if (ind == 2) {                                  // score in (i,j) stems from a deletion in sequence B
	      		I_i[i][j] = i;
			I_j[i][j] = j-1;
	       }
	       if (ind == 3) {                                  // (i,j) is the beginning of a subsequence
	      		I_i[i][j] = i;
			I_j[i][j] = j;	
	       }
        }
  }

  // search H for the maximal score
  double H_max = 0.;
  int i_max=0,j_max=0;
  for(int i=1;i<=N_a;i++){
    for(int j=1;j<=N_b;j++){
      if(H[i][j]>H_max){
	H_max = H[i][j];
	i_max = i;
	j_max = j;
      }
    }
  }

  score = H_max;

  // Backtracking from H_max
  int current_i=i_max,  current_j=j_max;
  int next_i=I_i[current_i][current_j];
  int next_j=I_j[current_i][current_j];
  int tick=0;
  char consensus_a[N_a+N_b+2], consensus_b[N_a+N_b+2];

  int offset;
  int query_end = current_i;

  int ref_end = j_max;

  while(((current_i!=next_i) || (current_j!=next_j)) && (next_j>=0) && (next_i>=0)){

        if(next_i==current_i)  consensus_a[tick] = '-';                  // deletion in A
        else                   consensus_a[tick] = seq_a[current_i-1];   // match/mismatch in A
    
        if(next_j==current_j)  consensus_b[tick] = '-';                  // deletion in B
        else                   consensus_b[tick] = seq_b[current_j-1];   // match/mismatch in B
    
        // mkroon: fix for adding first character of the alignment.
        if (next_i == 0) {
            next_i = -1;
        } else if (next_j == 0) {
            next_j = -1;
        } else {
            current_i = next_i;
            current_j = next_j;
            next_i = I_i[current_i][current_j];
            next_j = I_j[current_i][current_j];
        }
        tick++;
   }

   offset = tick;
 
 //Output of the consensus motif to the console
/*std::cout<<"\n"<<"***********************************************"<<"\n";
std::cout<<"The alignment of the sequences"<<"\n"<<"\n";
for(int i=0;i<N_a;i++){std::cout<<seq_a[i];}; std::cout<<"  and"<<"\n";
for(int i=0;i<N_b;i++){std::cout<<seq_b[i];}; std::cout<<"\n" <<"\n";
for(int i=tick-1;i>=0;i--) std::cout<<consensus_a[i]; std::cout<<"\n";
for(int j=tick-1;j>=0;j--) std::cout<<consensus_b[j]; std::cout<<"\n";*/

    for (int i=tick-1; i>=0; i--) {
	if (consensus_a[i] == '-' || consensus_b[i] == '-') {
	    gaps++;
	}

        if ((consensus_a[i] != consensus_b[i]) && (consensus_a[i] != '-' && consensus_b[i] != '-')) {
            mismatches++;
        }
    }


 int query_start = query_end - offset;
 int ref_start   = ref_end - offset;
 
 output.push_back(query_start);
 output.push_back(query_end-1);
 output.push_back(ref_start);
 output.push_back(ref_end-1);

 return output;

}





