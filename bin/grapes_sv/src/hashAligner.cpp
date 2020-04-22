#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <algorithm>
#include <tuple>
#include <sstream>
#include <regex>
#include <fstream>
#include "hashAligner.h"
#include "seed_t.h"
#include <set>
#include <list>

using namespace std;

string ToStr( char c ) {
   return string( 1, c );
}

struct node {
	int start;
	int end;
	int refStart;
	int refEnd;	

	bool operator<(const node& x) const {
		return (start < x.start);
	}
};


// data structure to store graph edges
struct Edge {
	node src;
	node dest;
	int weight;
};


// class to represent a graph object
class Graph {
	public:
		// construct a vector of vectors of Pairs to represent an adjacency list
		map<node, vector<node> > adjList;

		// Graph Constructor
		Graph(vector<Edge> const &edges, int N) {

			// add edges to the directed graph
			for (auto &edge: edges) {
				node src    = edge.src;
				node dest   = edge.dest;

				// insert at the end
				adjList[src].push_back(dest);
			}

			/*for (auto& n : adjList) {

				for (auto& m : n.second) {
					cout << n.first.start << "\t" << m.start << endl;
				}
			}*/
		}
};



void hashAligner::createDAG(std::vector<seed_t>& intervalVector ) {

	std::set< pair<node,node> > Nodes;
	vector<Edge> edges;

	for (auto& vertex1 : intervalVector) {
		cout << "Vertex start: " << vertex1.start << "-" << vertex1.end << "\t "<< "RefSequence: " << vertex1.refStart << "-" << vertex1.refEnd << endl;
		for (auto& vertex2 : intervalVector) {
			node v1;
			v1.start = vertex1.start;
			v1.end   = vertex1.end;
			v1.refStart = vertex1.refStart;			
			v1.refEnd = vertex1.refEnd;			

            		node v2;
			v2.start = vertex2.start;
			v2.end   = vertex2.end;
			v2.refStart = vertex2.refStart;			
			v2.refEnd = vertex2.refEnd;

			if ( abs(vertex2.start - vertex1.end) < 8 ) {
				Edge edge;
				edge.src  = v1;
				edge.dest = v2;
				edge.weight = v1.end-v1.start;
				edges.push_back(edge);
				//Nodes.insert(make_pair(v1,v2));
			}
		}
	}

	// Number of nodes in the graph
	int N = intervalVector.size();

	// construct graph
	Graph graph(edges, N);

	// print adjacency list representation of graph
	//printGraph(graph, N);


	// Filling DAG
	/*for (auto& vertex1 : intervalVector) {

		cout << "Vertex start: " << vertex1.start << "-" << vertex1.end << "\t "<< "RefSequence: " << vertex1.refStart << "-" << vertex1.refEnd << endl;
		
		for (auto& vertex2 : intervalVector) {
			
           		node v1;
			v1.start = vertex1.start;
			v1.end   = vertex1.end;
			v1.refStart = vertex1.refStart;			
			v1.refEnd = vertex1.refEnd;			

            		node v2;
			v2.start = vertex2.start;
			v2.end   = vertex2.end;
			v2.refStart = vertex2.refStart;			
			v2.refEnd = vertex2.refEnd;	

			if ( abs(vertex2.start - vertex1.end) < 8 ) {

				cout << "Vertex is compatible " << abs(vertex2.start - vertex1.end) << endl;
				Nodes.insert(make_pair(v1,v2));
			}
		}
	}

	for (auto& vertex : intervalVector) {
		Add source and destination vertexs
 		node source;
		node destination;
		Nodes.insert(make_pair(source,vertex));
		Nodes.insert(make_pair(vertex,destination));
	}
		

	for(auto f : Nodes) {
	   cout << f.first.start << "-" << f.first.end << "\t" <<  f.second.start << "-" << f.second.end << endl;
		
	  // use f here
	}*/    
}


std::multimap<std::string, int> hashAligner::indexContig(std::string& contig) {

	std::multimap<std::string, int> contigKmerHash;

	for (int i=0; i<contig.length()-kSize+1; i++) {

		std::string kmer = contig.substr (i, kSize);

		contigKmerHash.insert(std::pair<std::string, int>(kmer, i)); 
	}
	return contigKmerHash;
}

vector<seed_t> hashAligner::mergeIntervals( std::vector<seed_t>& seedVector ) {

	// Create an empty vector of intervals 
    	std::vector<seed_t> s; 
	  
	//std::map<seed_t, seed_t> seen;
		
	s.push_back(seedVector[0]);

	for (int i = 0; i < seedVector.size(); i++) {

		seed_t backSeed = s.back(); 

		// if current interval is not overlapping with stack top, 
		// push it to the stack 
		if (backSeed.refEnd < seedVector[i].refStart-2) {
		    s.push_back(seedVector[i]); 
		}

		// Otherwise update the ending time of top if ending of current 
		// interval is more 
		else if (backSeed.refEnd < seedVector[i].refEnd) {
		    backSeed.end = seedVector[i].end;
		    backSeed.refEnd = seedVector[i].refEnd;
		    backSeed.length = backSeed.end-backSeed.start;	
	
		    s.pop_back(); 
		    s.push_back(backSeed); 
		} 
	}
	return s;
}

vector<seed_t> hashAligner::findSeeds(std::string& refSequence, std::string& contig, std::multimap<std::string, int>& contigKmerHash){

	vector<seed_t> seedVector;

	for (int i = 0; i < refSequence.length()-kSize+1; i++) {

		std::string refKmer = refSequence.substr(i, kSize);

		auto const& r = contigKmerHash.equal_range(refKmer);

		for (std::multimap<std::string,int>::iterator it=r.first; it!=r.second; ++it) {

			int seedStart = it->second;
			int seedEnd   = it->second+kSize-1;

			std::string extSeed = it->first;
			int mismatches = 0;

			seed_t seed;
			seed.start = seedStart;
			seed.end   = seedEnd;
			seed.length = seed.end-seed.start;
			seed.refStart = i;
			seed.refEnd   = i+kSize-1;
			seedVector.push_back(seed);
		}
	}

	// Sorting by contig-start position
	std::sort(seedVector.begin(), seedVector.end());
	
	for (auto& s : seedVector ) {
		//cout << "Ordered " << s.start << "-" << s.end << "\t" << s.refStart << "-" << s.refEnd << endl;
	}
	

	return seedVector;
}

