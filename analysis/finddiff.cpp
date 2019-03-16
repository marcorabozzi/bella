//=======================================================================
// Title:  C++ program to compare BELLA and MECAT output
// Author: G. Guidi
// Date:   14 Mar 2019
//=======================================================================

// Compile: g++ -std=c++14 -g -march=native -fopenmp -fpermissive -O3 finddiff.cpp -o diff
// Run: ./diff <BELLA-alignments> <BELLA-overlaps> <MECAT-overlaps> <MECAT-indexes> <output-name-1> <output-name-2>

// file1 contains BELLA's alignments and file2 contains BELLA's overlaps (alignment threashold at 0)
// MECAT's extra overlaps (compared to BELLA's alignments) are then searched amongst BELLA's overlaps
// The goal is to identify potential good overlaps we lose while filtering 			

// TODO: to be generalized for whatever overlapper		

#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <istream>
#include <iostream>
#include <vector>
#include <string>
#include <stdlib.h>
#include <algorithm>
#include <utility>
#include <array>
#include <tuple>
#include <queue>
#include <memory>
#include <stack>
#include <functional>
#include <cstring>
#include <string.h>
#include <math.h>
#include <sstream>
#include <iterator>
#include <cassert>
#include <numeric>
#include <ios>
#include <sys/stat.h>
#include <sys/types.h>
#include <sys/sysctl.h>
#include <map>
#include <unordered_map>
#include <omp.h>

//=======================================================================
// 
// Common functions
// 
//=======================================================================

typedef std::pair<std::string, std::string> _myOverlap;
typedef std::vector<_myOverlap> _myOverlapVector;

int estimate (int begpV, int endpV, int lenV, int begpH, int endpH, int lenH)
{
	int diffV = endpV - begpV;
	int diffH = endpH - begpH;
	int minL  = std::min(begpV, begpH);
	int minR  = std::min(lenV - endpV, lenH - endpH);
	int ovlen = minL + minR + (diffV + diffH)/2;

	return ovlen;
}

std::vector<std::string> split (const std::string &s, char delim)
{
	std::vector<std::string> result;
	std::stringstream ss (s);
	std::string item;

	while (std::getline (ss, item, delim))
	{
		result.push_back (item);
	}

	return result;
}

/* from mecat index file to a map<uint32_t,string> : map<index,read-name> */
void mecatidx (std::ifstream& idx2read, std::map<uint32_t, std::string>& names)
{
	std::string num, name, seq;
	uint32_t idx;

	if(idx2read.is_open())
	{   
		std::string line;
		while(getline(idx2read, line))
		{
			std::stringstream linestream(line);

			getline(linestream, num, ' ');
			getline(linestream, name, ' ' );

			/* sequence on new line */
			getline(idx2read, seq); 

			idx = stoi(num);
			/* remove first char'>' */
			name.erase(0, 1); 
			names.insert(std::make_pair(idx, name));
		}
		std::cout << "MECAT idx2read table created" << std::endl;
	}
	else
	{
		std::cout << "Error creating names table from idx2read" << std::endl;
		exit(1);
	}
}

/* from mecat numeric id to read name */
std::string idx2read(uint32_t idx, std::map<uint32_t, std::string>& names)
{
	std::map<uint32_t, std::string>::iterator it;
	std::string name;

	it = names.find(idx);
	if(it != names.end())
		return names[idx];
	else
	{
		std::cout << "Read " << idx << " not present in MECAT output" << std::endl;
		exit(1);
	}
}
//=======================================================================
// 
// main (TODO: break down in functions)
// 
//=======================================================================

int main(int argc, char const *argv[])
{
	std::ifstream alignb(argv[1]);	// BELLA alignment file
	std::ifstream overlb(argv[2]);	// BELLA overlap file
	std::ifstream alignm(argv[3]);	// MECAT alignment file
	std::ifstream indexm(argv[4]);	// MECAT index
	const char* out_b = argv[5];	// MECAT extras found in BELLA -- MECAT side
	const char* out_m = argv[6];	// MECAT extras found in BELLA -- BELLA side

	std::map<uint32_t, std::string> names;
	mecatidx (indexm, names);

	int maxt = 1;
#pragma omp parallel
	{
		maxt = omp_get_num_threads();
	}
	// count alignments from bella
	uint64_t n_alignb = std::count(std::istreambuf_iterator<char>(alignb), std::istreambuf_iterator<char>(), '\n');
	uint64_t n_overlb = std::count(std::istreambuf_iterator<char>(overlb), std::istreambuf_iterator<char>(), '\n');
	uint64_t n_alignm = std::count(std::istreambuf_iterator<char>(alignm), std::istreambuf_iterator<char>(), '\n');

	alignb.seekg(0,std::ios_base::beg);
	alignm.seekg(0,std::ios_base::beg);
	overlb.seekg(0,std::ios_base::beg); 		

	std::vector<std::string> lines_alignb;
	std::vector<std::string> lines_alignm;
	std::vector<std::string> lines_overlb;

	std::vector<std::stringstream> local_out_b(maxt);      
	std::vector<std::stringstream> local_out_m(maxt);      

	// read BELLA alignment file
	if(alignb)
		for (int i = 0; i < n_alignb; ++i)
		{
			std::string line;
			std::getline(alignb, line);
			lines_alignb.push_back(line);
		}
	alignb.close();

	// read MECAT alignment file
	if(alignm)
		for (int i = 0; i < n_alignm; ++i)
		{
			std::string line;
			std::getline(alignm, line);
			lines_alignm.push_back(line);
		}
	alignm.close();

	// read BELLA overlap file
	if(overlb)
		for (int i = 0; i < n_overlb; ++i)
		{
			std::string line;
			std::getline(overlb, line);
			lines_overlb.push_back(line);
		}
	overlb.close();

	_myOverlapVector global_alignb(n_alignb*2);	
	_myOverlapVector global_alignm(n_alignm*2);	
	_myOverlapVector global_overlb(n_overlb*2);	
	std::vector<_myOverlapVector> local_alignb(maxt);   
	std::vector<_myOverlapVector> local_alignm(maxt);  
	std::vector<_myOverlapVector> local_overlb(maxt);  

#pragma omp parallel
	{

	#pragma omp for nowait
		for(uint64_t i = 0; i < n_alignb; i++) 
		{
			int MYTHREAD = omp_get_thread_num();

			std::vector<std::string> v = split (lines_alignb[i], '\t');
			std::string id1 = v[0]; 	// BELLA has names in position 0, 1
			std::string id2 = v[1];	 	// BELLA has names in position 0, 1

			local_alignb[MYTHREAD].push_back(std::make_pair(id1,id2));
			local_alignb[MYTHREAD].push_back(std::make_pair(id2,id1));
		}

	#pragma omp for nowait
		for(uint64_t j = 0; j < n_alignm; j++) 
		{
			int MYTHREAD = omp_get_thread_num();

			std::vector<std::string> v = split (lines_alignm[j], '\t');
			// mecat idx to nametag translation 	
			std::string id1 = idx2read (stoi(v[0]), names);
			std::string id2 = idx2read (stoi(v[1]), names);

			local_alignm[MYTHREAD].push_back(std::make_pair(id1,id2));
			local_alignm[MYTHREAD].push_back(std::make_pair(id2,id1));
		}

	#pragma omp for nowait
		for(uint64_t j = 0; j < n_overlb; j++) 
		{
			int MYTHREAD = omp_get_thread_num();

			std::vector<std::string> v = split (lines_overlb[j], '\t');
			// mecat idx to nametag translation 	
			std::string id1 = v[0];
			std::string id2 = v[1];

			local_overlb[MYTHREAD].push_back(std::make_pair(id1,id2));
			local_overlb[MYTHREAD].push_back(std::make_pair(id2,id1));
		}
	}

	uint64_t alignmentssofar_b = 0;
	uint64_t alignmentssofar_m = 0;
	uint64_t overlapssofar_b = 0;
	for(int t = 0; t < maxt; ++t)
	{
		std::copy(local_alignb[t].begin(), local_alignb[t].end(), global_alignb.begin() + alignmentssofar_b);
		std::copy(local_alignm[t].begin(), local_alignm[t].end(), global_alignm.begin() + alignmentssofar_m);
		std::copy(local_overlb[t].begin(), local_overlb[t].end(), global_overlb.begin() + overlapssofar_b);
		alignmentssofar_b += local_alignb[t].size();
		alignmentssofar_m += local_alignm[t].size();
		overlapssofar_b += local_overlb[t].size();
	}

	std::sort(global_alignb.begin(), global_alignb.end());
	std::sort(global_alignm.begin(), global_alignm.end());
	std::sort(global_overlb.begin(), global_overlb.end());

	// the difference is the number of elements that are present in the first set, but not in the second one 
	_myOverlapVector diff1(n_alignb*2);		
	_myOverlapVector::iterator it1;	

	it1 = std::set_difference (global_alignb.begin(), global_alignb.end(), 
							   global_alignm.begin(), global_alignm.end(), 
							   diff1.begin());

	diff1.resize(it1-diff1.begin());

	// the difference is the number of elements that are present in the first set, but not in the second one 
	_myOverlapVector extras(n_alignm*2);	
	_myOverlapVector::iterator it2;	

	it2 = std::set_difference (global_alignm.begin(), global_alignm.end(), 
							   global_alignb.begin(), global_alignb.end(), 
							   extras.begin());

	extras.resize(it2-extras.begin()); 	// these are extras in MECAT (overlaps in MECAT's output but not in BELLA's)

	std::cout << "Alignments in BELLA : " << global_alignb.size() << std::endl;
	std::cout << "Overlaps in BELLA (alignment threshold at 0) : " << global_overlb.size() << std::endl;
	std::cout << "Alignemnts in MECAT : " << global_alignm.size() << std::endl;
	std::cout << "BELLA overlaps missing in MECAT : " << diff1.size() << std::endl;
	std::cout << "MECAT overlaps missing in BELLA : " << extras.size() << std::endl;

	// now we're gonna look for these overlaps in BELLA's overlap outputs (alignment threshold at 0)
	// elements of extras with are in overlb
	std::sort(extras.begin(), extras.end());
	_myOverlapVector extrasset(extras.size());
	_myOverlapVector::iterator it3;

	it3 = std::set_intersection(global_overlb.begin(), global_overlb.end(),
						  extras.begin(), extras.end(),
						  extrasset.begin()); //	container that supports a push_back operation

	extrasset.resize(it3-extrasset.begin());
	// extrasset cannot be greater than extras
	std::cout << "MECAT overlaps found in BELLA overlaps : " << extrasset.size() << std::endl;

#pragma omp parallel
	{
	#pragma omp for nowait
		for(uint64_t i = 0; i < n_overlb; i++) 
		{
			_myOverlapVector::iterator it;		
			it = std::find (extrasset.begin(), extrasset.end(), global_overlb[i]);

			if (it == extrasset.end())
				continue;

			int MYTHREAD = omp_get_thread_num();
			std::vector<std::string> v = split (lines_overlb[i], '\t');
			// mecat idx to nametag translation 	
			std::string id1 	= v[0];
			std::string id2 	= v[1];
			std::string score 	= v[3];
			std::string overlap = v[4];
			std::string start1 	= v[6];
			std::string end1 	= v[7];
			std::string len1 	= v[8];
			std::string start2 	= v[9];
			std::string end2 	= v[10];
			std::string len2 	= v[11];

			local_out_b[MYTHREAD] << id1 << "\t" << id2 << "\t" << start1 << "\t" << end1 << "\t" << len1
				<< "\t" << start2 << "\t" << end2 << "\t" << len2 << "\t" << score << "\t" << overlap << std::endl;
		}
	}

	// write BELLA output to new file
	int64_t * bytes1 = new int64_t[maxt];
	for(int i = 0; i < maxt; ++i)
	{
		local_out_b[i].seekg(0, std::ios::end);
		bytes1[i] = local_out_b[i].tellg();
		local_out_b[i].seekg(0, std::ios::beg);
	}
	int64_t bytestotal1 = std::accumulate(bytes1, bytes1 + maxt, static_cast<int64_t>(0));

	std::ofstream output1(out_b, std::ios::binary | std::ios::app);
	cout << "Creating or appending to output file with " << (double)bytestotal1/(double)(1024 * 1024) << " MB" << endl;
	output1.seekp(bytestotal1 - 1);
	/* this will likely create a sparse file so the actual disks won't spin yet */
	output1.write("", 1); 
	output1.close();

	#pragma omp parallel
	{
		int ithread = omp_get_thread_num(); 

		FILE *ffinal1;
		/* then everyone fills it */
		if ((ffinal1 = fopen(out_b, "rb+")) == NULL) 
		{
			fprintf(stderr, "File %s failed to open at thread %d\n", out_b, ithread);
		}
		int64_t bytesuntil1 = std::accumulate(bytes1, bytes1 + ithread, static_cast<int64_t>(0));
		fseek (ffinal1 , bytesuntil1 , SEEK_SET);
		std::string text1 = local_out_b[ithread].str();
		fwrite(text1.c_str(),1, bytes1[ithread], ffinal1);
		fflush(ffinal1);
		fclose(ffinal1);
	}
	delete [] bytes1;

#pragma omp parallel
	{
	#pragma omp for nowait
		for(uint64_t i = 0; i < n_alignm; i++) 
		{
			_myOverlapVector::iterator it;		
			it = std::find (extrasset.begin(), extrasset.end(), global_alignm[i]);

			if (it == extrasset.end())
				continue;

			int MYTHREAD = omp_get_thread_num();
			std::vector<std::string> v = split (lines_overlb[i], '\t');

			// mecat idx to nametag translation 	
			std::string id1 = idx2read (stoi(v[0]), names);
			std::string id2 = idx2read (stoi(v[1]), names);
			std::string ident 	= v[2];
			std::string start1 	= v[5];
			std::string end1 	= v[6];
			std::string len1 	= v[7];
			std::string start2 	= v[9];
			std::string end2	= v[10];
			std::string len2 	= v[11];

			int ovlen = estimate (stoi(start1), stoi(end1), stoi(len1), stoi(start2), stoi(end2), stoi(len2));
			int score = floor((stod(ident)*ovlen) / 100);

			local_out_m[MYTHREAD] << id1 << "\t" << id2 << "\t" << start1 << "\t" << end1 << "\t" << len1
				<< "\t" << start2 << "\t" << end2 << "\t" << len2 << "\t" << score << "\t" << ovlen << std::endl;
		}
	}

	// write MECAT output to new file 	
	int64_t * bytes2 = new int64_t[maxt];
	for(int i = 0; i < maxt; ++i)
	{
		local_out_m[i].seekg(0, std::ios::end);
		bytes2[i] = local_out_m[i].tellg();
		local_out_m[i].seekg(0, std::ios::beg);
	}
	int64_t bytestotal2 = std::accumulate(bytes2, bytes2 + maxt, static_cast<int64_t>(0));

	std::ofstream output2(out_m, std::ios::binary | std::ios::app);
	std::cout << "Creating or appending to output file with " << (double)bytestotal2/(double)(1024 * 1024) << " MB" << std::endl;
	output2.seekp(bytestotal2 - 1);
	/* this will likely create a sparse file so the actual disks won't spin yet */
	output2.write("", 1); 
	output2.close();

	#pragma omp parallel
	{
		int ithread = omp_get_thread_num(); 

		FILE *ffinal2;
		/* then everyone fills it */
		if ((ffinal2 = fopen(out_m, "rb+")) == NULL) 
		{
			fprintf(stderr, "File %s failed to open at thread %d\n", out_m, ithread);
		}
		int64_t bytesuntil2 = std::accumulate(bytes2, bytes2 + ithread, static_cast<int64_t>(0));
		fseek (ffinal2, bytesuntil2 , SEEK_SET);
		std::string text2 = local_out_m[ithread].str();
		fwrite(text2.c_str(),1, bytes2[ithread], ffinal2);
		fflush(ffinal2);
		fclose(ffinal2);
	}
	delete [] bytes2;

	return 0;
}