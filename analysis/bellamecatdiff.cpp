//=======================================================================
// Title:  C++ program to compare BELLA and MECAT output
// Author: G. Guidi
// Date:   14 Mar 2019
//=======================================================================

// Compile: g++ -std=c++14 -g -march=native -fopenmp -fpermissive -O3 bellamecatdiff.cpp -o diff
// Run: ./diff <file1> <file2> <file3>
// file1 is supposed to be the reference (BELLA), file2 the query (MECAT), and file3 MECAT index file

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
#include <ios>
#include <sys/stat.h>
#include <sys/types.h>
#include <sys/sysctl.h>
#include <map>
#include <unordered_map>
#include <omp.h>

typedef std::pair<std::string, std::string> _myOverlap;
typedef std::vector<_myOverlap> _myOverlapVector;

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

int main(int argc, char const *argv[])
{
	std::ifstream file1(argv[1]);	// BELLA overlap
	std::ifstream file2(argv[2]);	// MECAT overlap
	std::ifstream file3(argv[3]);	// MECAT index

	std::map<uint32_t, std::string> names;
	mecatidx (file3, names);

	//	map<uint32_t, std::string> names;
	//	mecatidx (index, names);

	int maxt = 1;
#pragma omp parallel
	{
		maxt = omp_get_num_threads();
	}

	uint64_t numoverlap1 = std::count(std::istreambuf_iterator<char>(file1), std::istreambuf_iterator<char>(), '\n');
	uint64_t numoverlap2 = std::count(std::istreambuf_iterator<char>(file2), std::istreambuf_iterator<char>(), '\n');

	file1.seekg(0,std::ios_base::beg);
	file2.seekg(0,std::ios_base::beg); 

	std::vector<std::string> lines1;
	std::vector<std::string> lines2;

	/* Read file1 */
	if(file1)
		for (int i = 0; i < numoverlap1; ++i)
		{
			std::string line;
			std::getline(file1, line);
			lines1.push_back(line);
		}
	file1.close();

	/* Read file1 */
	if(file2)
		for (int i = 0; i < numoverlap2; ++i)
		{
			std::string line;
			std::getline(file2, line);
			lines2.push_back(line);
		}
	file2.close();

	_myOverlapVector global1(numoverlap1*2);	
	_myOverlapVector global2(numoverlap2*2);	
	std::vector<_myOverlapVector> local1(maxt);   
	std::vector<_myOverlapVector> local2(maxt);  

#pragma omp parallel
	{

	#pragma omp for nowait
		for(uint64_t i = 0; i < numoverlap1; i++) 
		{
			std::stringstream linestream(lines1[i]);
			int tid = omp_get_thread_num();

			std::vector<std::string> v = split (lines1[i], '\t');
			std::string id1 = v[0]; 	// BELLA has names in position 0, 1
			std::string id2 = v[1];	 	// BELLA has names in position 0, 1

			local1[tid].push_back(std::make_pair(id1,id2));
			local1[tid].push_back(std::make_pair(id2,id1));
		}

	#pragma omp for nowait
		for(uint64_t j = 0; j < numoverlap2; j++) 
		{
			int MYTHREAD = omp_get_thread_num();

			std::vector<std::string> v = split (lines2[j], '\t');
			// mecat idx to nametag translation 	
			std::string id1 = idx2read (stoi(v[0]), names);
			std::string id2 = idx2read (stoi(v[1]), names);

			local2[MYTHREAD].push_back(std::make_pair(id1,id2));
			local2[MYTHREAD].push_back(std::make_pair(id2,id1));
		}
	}

	uint64_t overlapssofar1 = 0;
	uint64_t overlapssofar2 = 0;
	for(int t = 0; t < maxt; ++t)
	{
		std::copy(local1[t].begin(), local1[t].end(), global1.begin() + overlapssofar1);
		std::copy(local2[t].begin(), local2[t].end(), global2.begin() + overlapssofar2);
		overlapssofar1 +=	local1[t].size();
		overlapssofar2 +=	local2[t].size();
	}

	std::sort(global1.begin(), global1.end());
	std::sort(global2.begin(), global2.end());

	// the difference is the number of elements that are present in the first set, but not in the second one 
	_myOverlapVector diff1(numoverlap1*2);	
	_myOverlapVector::iterator it1;	

	it1 = std::set_difference (global1.begin(), global1.end(), 
							   global2.begin(), global2.end(), 
							   diff1.begin());

	diff1.resize(it1-diff1.begin());

	// the difference is the number of elements that are present in the first set, but not in the second one 
	_myOverlapVector diff2(numoverlap2*2);	
	_myOverlapVector::iterator it2;	

	it2 = std::set_difference (global2.begin(), global2.end(), 
							   global1.begin(), global1.end(), 
							   diff2.begin());

	diff2.resize(it2-diff2.begin());

	std::cout << "Overlaps in file1 : " << global1.size() << std::endl;
	std::cout << "Overlaps in file2 : " << global2.size() << std::endl;
	std::cout << "Overlaps missing in MECAT : " << diff1.size() << std::endl;
	std::cout << "Overlaps extra in MECAT   : " << diff2.size() << std::endl;

	return 0;
}