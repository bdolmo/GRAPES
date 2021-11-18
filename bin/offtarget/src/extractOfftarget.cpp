
#include <iostream>
#include <string.h>
#include <string>
#include <vector>
#include <iostream>
#include <map>
#include <regex>
#include <string.h>
#include <numeric>
#include <algorithm>
#include "bed_t.h"

#include <math.h>       /* log10 */
#include "SeqLib/BamHeader.h"
#include "SeqLib/RefGenome.h"
#include "SeqLib/BamReader.h"
#include "SeqLib/BamRecord.h"
#include "SeqLib/BamWriter.h"
#include "SeqLib/BWAWrapper.h"

using namespace SeqLib;
using namespace std;

float computeGC ( std::string chr, int start, int end, RefGenome& reference ) {

	std::string sequence = reference.QueryRegion(chr, start, end);
	int count = 0;
	for (auto& ntd : sequence) {
		if (ntd == 'C' || ntd == 'c' || ntd == 'G' || ntd == 'g' ) {
			count++;
		}
	}
	float gc_content = (count/ (float) sequence.length()) * (float) 100;
	return gc_content;
}


vector<string> inline split( std::string const& original, char separator )
{
    std::vector<std::string> results;
    std::string::const_iterator start = original.begin();
    std::string::const_iterator end = original.end();
    string::const_iterator next = std::find( start, end, separator );
    while ( next != end ) {
        results.push_back( string( start, next ) );
        start = next + 1;
        next = std::find( start, end, separator );
    }
    results.push_back( std::string( start, next ) );
    return results;
}

void bed2Vector( std::string& bed, vector<bed_t>& RegionsVector, RefGenome& ref) {
	  ifstream bedFile;
  	bedFile.open (bed);
    if (bedFile.is_open()) {
			string line;
			while ( std::getline (bedFile, line)) {
				vector<string> tmp = split(line, '\t');
				string chr = tmp[0];
				int start  = atoi(tmp[1].c_str());
				int end    = atoi(tmp[2].c_str());

				float gc;
				bed_t region;
				region.chr  = chr;
				region.start= start;
				region.end  = end;
				region.gc   = gc;

				if (tmp.size() > 3) {
					for (int i =3; i < tmp.size();i++) {
						region.fields.push_back(tmp[i]);
					}
					string info = tmp[3];
					region.info = info;
				}
				RegionsVector.push_back(region);
			}
		}
}


int main (int argc, char* argv[]) {

	string bamFile(argv[1]);
	string reference(argv[2]);
	string bed(argv[3]);

	BamReader reader;
	reader.Open(bamFile);
	BamRecord r;

	BamHeader myHeader;
	myHeader = reader.Header();

	RefGenome ref;
	ref.LoadIndex(reference);

	string chrA;
	int startA;
	int endA;

	string chrB;
	int startB;
	int endB;

	int first_start;
	int last_end;

	int i = 0;
	int reads_in_peak = 0;
	long genomeSize = 2700000000;

	vector<bed_t> RegionsVector;
	bed2Vector(bed, RegionsVector, ref);

	int count = 0;
	int tag = 0;

	string chr_region;
	int start_region;
	int end_region;
	string info_region;
	string fields_region;
	while (reader.GetNextRecord(r)) {

		int chr = r.ChrID();
		if (chr < 0 || chr > 24) {
			continue;
		}

		string chrom  = myHeader.IDtoName(r.ChrID());
    //cout << chrom << "\t" << r.Position() << endl;
		if (i == 0) {
			chr_region  = RegionsVector[i].chr;
			start_region= RegionsVector[i].start;
			end_region  = RegionsVector[i].end;
			fields_region = "";
			for (auto& field : RegionsVector[i].fields) {
				fields_region+= "\t" + field;
			}
		}

		// if ( !reader.GetNextRecord(r)) {
		// 	cout <<  i << " " << r.Position() << " " << r.Position()+r.Length() << " "<<RegionsVector.size()-1 << " "<< chr_region << "\t" << start_region << "\t" << end_region  << fields_region << "\t" << count << endl;
		// }


		if ( (r.MapQuality() > 50) && (chrom == chr_region) && ( ( r.Position() <= start_region) && (r.Position()+r.Length() > start_region )
		|| ( r.Position() < end_region && r.Position()+r.Length() >= end_region)
		|| ( r.Position() >= start_region && r.Position()+r.Length() <= end_region))) {
			count++;
		}
		if ( chrom == chr_region && r.Position() > end_region && i < RegionsVector.size()-1) {
			cout << chr_region << "\t" << start_region << "\t" << end_region  << fields_region << "\t" << count << endl;
			i++;
			count = 0;
			chr_region  = RegionsVector[i].chr;
			start_region= RegionsVector[i].start;
			end_region  = RegionsVector[i].end;
			fields_region = "";
			for (auto& field : RegionsVector[i].fields) {
				fields_region+= "\t" + field;
			}
			continue;
		}
		if (chrom != chr_region && i < RegionsVector.size()-1 ) {

			cout << chr_region << "\t" << start_region << "\t" << end_region  << fields_region << "\t" << count << endl;
			i++;
			count = 0;
			chr_region  = RegionsVector[i].chr;
			start_region= RegionsVector[i].start;
			end_region  = RegionsVector[i].end;
			fields_region = "";
			for (auto& field : RegionsVector[i].fields) {
				fields_region+= "\t" + field;
			}

			if ( (r.MapQuality() > 50) && (chrom == chr_region) && ( ( r.Position() <= start_region) && (r.Position()+r.Length() > start_region )
			|| ( r.Position() < end_region && r.Position()+r.Length() >= end_region)
			|| ( r.Position() >= start_region && r.Position()+r.Length() <= end_region))) {
				count++;
			}
			continue;
		}

		if ( i ==  RegionsVector.size()-1 || !reader.GetNextRecord(r)) {
			if ( (r.MapQuality() > 50) && (chrom == chr_region) && ( ( r.Position() <= start_region) && (r.Position()+r.Length() > start_region )
			|| ( r.Position() < end_region && r.Position()+r.Length() >= end_region)
			|| ( r.Position() >= start_region && r.Position()+r.Length() <= end_region))) {
				count++;
			}
			// if ( !reader.GetNextRecord(r)) {
			// 	cout <<  i << " " <<RegionsVector.size()-1 << " "<< chr_region << "\t" << start_region << "\t" << end_region  << fields_region << "\t" << count << endl;
			// }
		}

	}
	for (int j=i;j<=RegionsVector.size()-1;j++)  {
		count = 0;
		chr_region  = RegionsVector[j].chr;
		start_region= RegionsVector[j].start;
		end_region  = RegionsVector[j].end;
		fields_region = "";
		for (auto& field : RegionsVector[j].fields) {
			fields_region+= "\t" + field;
		}
		cout << chr_region << "\t" << start_region << "\t" << end_region  << fields_region << "\t" << count << endl;

	}

	return 0;
}
