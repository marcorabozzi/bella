//=======================================================================
// Title:  C++ program to classify spanning entries
// Author: G. Guidi
// Date:   31 Mar 2019
//=======================================================================

// Compile: g++ -std=c++14 -g -march=native -fopenmp -fpermissive -O3 spanningRead.cpp -o span
// Run: ./span <read-to-contig> 

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

typedef std::map<std::string, std::vector< std::string >> mymap;

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
	std::ifstream file1(argv[1]);	// read to contig

	int maxt = 1;
#pragma omp parallel
	{
		maxt = omp_get_num_threads();
	}

	uint64_t numentry = std::count(std::istreambuf_iterator<char>(file1), std::istreambuf_iterator<char>(), '\n');
	file1.seekg(0,std::ios_base::beg);
	std::vector<std::string> lines1;

	/* Read file1 */
	if(file1)
		for (int i = 0; i < numentry; ++i)
		{
			std::string line;
			std::getline(file1, line);
			lines1.push_back(line);
		}
	file1.close();

	mymap readkey;	
	mymap contigkey;	
	std::vector<mymap> local1(maxt);   
	std::vector<mymap> local2(maxt);   

#pragma omp parallel
	{
	#pragma omp for
		for(uint64_t i = 0; i < numentry; i++) 
		{
			int tid = omp_get_thread_num();

			std::vector<std::string> v = split (lines1[i], '\t');
			std::string read = v[0]; 	
			std::string contig = v[1];	

			std::vector<std::string> temp1;
			std::vector<std::string> temp2;
			temp1.push_back(contig);
			temp2.push_back(read);

			mymap::iterator it1 = local1[tid].find(read);
			mymap::iterator it2 = local2[tid].find(contig);

			// read = key
			if(it1 == local1[tid].end())
				local1[tid].insert( std::make_pair ( read, temp1 ) );
			else
				local1[tid][read].push_back( contig );

			// contig = key 
			if(it2 == local2[tid].end())
				local2[tid].insert( std::make_pair ( contig, temp2 ) );
			else
				local2[tid][contig].push_back( read );
		}
	}

	for(int t = 0; t < maxt; ++t)
	{
		readkey.insert(local1[t].begin(), local1[t].end());
		contigkey.insert(local2[t].begin(), local2[t].end());
	}

	std::cout << "num contigs : " << contigkey.size() << std::endl;
	std::cout << "num reads supplementary alignment : " << readkey.size() << std::endl;

	for (mymap::iterator it = readkey.begin(); it != readkey.end(); it++) 	// iterating over reads
	{
		if(it->second.size() == 1)   	// if read aligns to a single contig
		{
			mymap::iterator it3 = contigkey.find(it->second.front()); 	// it does exist
			std::vector<std::string>& vec = it3->second;			// vector of read of this contig
			//std::cout << "before " << vec.size() << std::endl;
			vec.erase(std::remove(vec.begin(), vec.end(), it->first), vec.end()); 	// remove this read from the vector of read of this contig
			//std::cout << "after " << vec.size() << std::endl;
			readkey.erase(it->first); 	// remove the read from the read map
		}
	}

	std::cout << "num spanning reads supplementary alignment : " << readkey.size() << std::endl;
	std::map<std::pair<std::string, std::string>, std::pair<std::string, std::string>> span;

	for (auto const& i : readkey) 			// iterating over reads
	{
		std::string origin = i.second.front();
		std::next(&i.second, 1);

		for (auto const& j : i.second) 		// iterating over contigs (from the second one)
		{
			auto c = contigkey.find(j);
			for(auto const& k : c->second)	// iterating over reads of the second contig
				if(k != i.first)
				{
					std::pair<std::string, std::string> pair1 = make_pair(i.first, k);		// pair of spanning reads
					std::pair<std::string, std::string> pair2 = make_pair(origin, j);	// pair of spanned contigs
					span.insert(make_pair( pair1, pair2));
				}
		}
	}

	std::cout << "num spanning pair : " << span.size() << std::endl;

	return 0;
}