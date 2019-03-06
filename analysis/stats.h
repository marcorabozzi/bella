#include <omp.h>
#include <fstream>
#include <iostream>
#include <string>
#include <cstdlib>
#include <cstdio>
#include <algorithm>
#include <cmath>
#include <numeric>
#include <vector>
#include <sys/types.h> 
#include <sys/stat.h> 
#include <math.h>
#include <limits.h>
#include <bitset>
#include <map>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <unistd.h>
#include <ctype.h> 
#include <sstream>
#include <set>
#include <memory>
#include <typeinfo>

using namespace std;

/* PAF format: https://github.com/lh3/miniasm/blob/master/PAF.md
 * column seq name
 * column seq length
 * column seq start
 * column seq end
 * strand (+/-)
 * row seq name
 * row seq length
 * row seq start
 * row seq end
 * number of residue matches (alignment score)
 * alignment block length (overlap length)
 * mapping quality (0-255; 255 for missing) */

//=======================================================================
// 
// Common functions
// 
//=======================================================================

typedef std::vector< std::pair< std::string,std::string > > pairv;
typedef std::tuple< pairv, double > alignstats;

int estimate (int begpV, int endpV, int lenV, int begpH, int endpH, int lenH)
{
	int diffV = endpV - begpV;
	int diffH = endpH - begpH;
	int minL  = min(begpV, begpH);
	int minR  = min(lenV - endpV, lenH - endpH);
	int ovlen = minL + minR + (diffV + diffH)/2;

	return ovlen;
}

vector<std::string> split (const std::string &s, char delim)
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
void mecatidx (ifstream& idx2read, map<uint32_t, std::string>& names)
{
	string num, name, seq;
	uint32_t idx;

	if(idx2read.is_open())
	{   
		string line;
		while(getline(idx2read, line))
		{
			stringstream linestream(line);

			getline(linestream, num, ' ');
			getline(linestream, name, ' ' );

			/* sequence on new line */
			getline(idx2read, seq); 

			idx = stoi(num);
			/* remove first char'>' */
			name.erase(0, 1); 
			names.insert(std::make_pair(idx, name));
		}
		cout << "MECAT idx2read table created" << endl;
	}
	else
	{
		cout << "Error creating names table from idx2read" << endl;
		exit(1);
	}
}

/* from mecat numeric id to read name */
std::string idx2read(uint32_t idx, map<uint32_t, std::string>& names)
{
	map<uint32_t, std::string>::iterator it;
	string name;

	it = names.find(idx);
	if(it != names.end())
		return names[idx];
	else
	{
		cout << "Read " << idx << " not present in MECAT output" << endl;
		exit(1);
	}
}

//=======================================================================
// 
// BELLA statistics
// 
//=======================================================================

alignstats BELLA2PAF(ifstream& input, char* filename)
{
	int maxt = 1;
#pragma omp parallel
	{
		maxt = omp_get_num_threads();
	}

	uint64_t numoverlap = std::count(std::istreambuf_iterator<char>(input), std::istreambuf_iterator<char>(), '\n');
	input.seekg(0, std::ios_base::beg);

	vector<std::string> entries;
	vector<pairv> local(maxt);      
	pairv pairs;
	alignstats result;
	vector<std::vector<double>> identity(maxt);      
	vector<double> avgsvec(maxt);
	double avg;

	/* read input file */
	if(input)
		for (int i = 0; i < numoverlap; ++i)
		{
			std::string line;
			std::getline(input, line);
			entries.push_back(line);
		}
	input.close();

	/* get stats from BELLA output in PAF format */
#pragma omp parallel
{
	#pragma omp for
	for(uint64_t i = 0; i < numoverlap; i++) 
	{
		int ithread = omp_get_thread_num();
		std::vector<std::string> v = split (entries[i], '\t');

		/* improve readability */
		std::string& nameV = v[0];
		std::string& lengV = v[1];
		std::string& begpV = v[2];
		std::string& endpV = v[3];
		std::string& isRev = v[4];
		std::string& nameH = v[5];
		std::string& lengH = v[6];
		std::string& begpH = v[7];
		std::string& endpH = v[8];
		std::string& score = v[9];
		std::string& ovlen = v[10];
		std::string& mapQV = v[11];

		/* compute pair identity */
		int matches = floor((stoi(ovlen) + stoi(score)) / 2); // we might need to lower down this a bit    	
		double id = (float(matches)/stof(score)) * 100;

		identity[ithread].push_back(id);

		local[ithread].push_back(std::make_pair(nameV, nameH));
		local[ithread].push_back(std::make_pair(nameH, nameV));		// complementaty pair 	
	}

	avgsvec[ithread] = std::reduce(std::execution::par, identity[ithread].begin(), identity[ithread].end());
	avgsvec[ithread] = avgsvec[ithread] / identity[ithread].size();
}

	avg = std::reduce(std::execution::par, avgsvec.begin(), avgsvec.end());
	avg = avg / maxt;

	size_t pairssofar = 0;
	for(int t = 0; t < maxt; ++t)
	{
		std::copy(local[t].begin(), local[t].end(), pairs.begin()+pairssofar);
		pairssofar += local[t].size();
	}

	return make_tuple(pairs, avg);
}

//=======================================================================
// 
// MHAP statistics
// 
//=======================================================================

//void MHAP2PAF(ifstream& input, char* filename)
//{
//	int maxt = 1;
//#pragma omp parallel
//	{
//		maxt = omp_get_num_threads();
//	}
//
//	uint64_t numoverlap = std::count(std::istreambuf_iterator<char>(input), std::istreambuf_iterator<char>(), '\n');
//	input.seekg(0, std::ios_base::beg);
//
//	vector<std::string> entries;
//	vector<std::stringstream> local(maxt);      
//
//	/* read input file */
//	if(input)
//		for (int i = 0; i < numoverlap; ++i)
//		{
//			std::string line;
//			std::getline(input, line);
//			entries.push_back(line);
//		}
//	input.close();
//
//	/* transform MHAP output in PAF format */
//#pragma omp parallel for
//	for(uint64_t i = 0; i < numoverlap; i++) 
//	{
//		int ithread = omp_get_thread_num();
//		/* MHAP format: cname, rname, err, nkmer, cstrand, cstart, cend, clen, rstrand, rstart, rend, rlen */
//		std::vector<std::string> v = split (entries[i], ' ');
//		/* improve readability */
//		std::string& nameV = v[0];
//		std::string& nameH = v[1];
//		std::string& begpV = v[5];
//		std::string& endpV = v[6];
//		std::string& lengV = v[7];
//		std::string& isRev = v[8];
//		std::string& begpH = v[9];
//		std::string& endpH = v[10];
//		std::string& lengH = v[11];
//
//		/* change strand formatting */
//		if(isRev == "0") isRev = "+";         
//			else isRev = "-";
//
//		/* compute overlap length if missing (begpV, endpV, lenV, begpH, endpH, lenH) */
//		int ovlen = estimate (stoi(begpV), stoi(endpV), stoi(lengV), stoi(begpH), stoi(endpH), stoi(lengH));
//
//		/* GGGG: If alignment is missing I estimate it as % of the overlap length and I determine that % using the error rate */
//		// GGGG: Error rate is now hard-coded, need to be an input parameter
//
//		float error = 0.15;
//		float identity = (1-error)*(1-error);
//		int score = floor(identity*ovlen);
//
//		local[ithread] << nameV << "\t" << lengV << "\t" << begpV << "\t" << endpV << "\t" << isRev 
//			<< "\t" << nameH << "\t" << lengH << "\t" << begpH << "\t" << endpH << "\t" << score 
//				<< "\t" << ovlen << "\t255" << endl;
//	}
//
//	/* write to a new file */
//	int64_t * bytes = new int64_t[maxt];
//	for(int i = 0; i < maxt; ++i)
//	{
//		local[i].seekg(0, ios::end);
//		bytes[i] = local[i].tellg();
//		local[i].seekg(0, ios::beg);
//	}
//	int64_t bytestotal = std::accumulate(bytes, bytes + maxt, static_cast<int64_t>(0));
//
//	std::ofstream output(filename, std::ios::binary | std::ios::app);
//#ifdef PRINT
//	cout << "Creating or appending to output file with " << (double)bytestotal/(double)(1024 * 1024) << " MB" << endl;
//#endif
//	output.seekp(bytestotal - 1);
//	/* this will likely create a sparse file so the actual disks won't spin yet */
//	output.write("", 1); 
//	output.close();
//
//	#pragma omp parallel
//	{
//		int ithread = omp_get_thread_num(); 
//
//		FILE *ffinal;
//		/* then everyone fills it */
//		if ((ffinal = fopen(filename, "rb+")) == NULL) 
//		{
//			fprintf(stderr, "File %s failed to open at thread %d\n", filename, ithread);
//		}
//		int64_t bytesuntil = std::accumulate(bytes, bytes + ithread, static_cast<int64_t>(0));
//		fseek (ffinal , bytesuntil , SEEK_SET);
//		std::string text = local[ithread].str();
//		fwrite(text.c_str(),1, bytes[ithread], ffinal);
//		fflush(ffinal);
//		fclose(ffinal);
//	}
//	delete [] bytes;
//}
//
////=======================================================================
//// 
//// MECAT statistics
//// 
////=======================================================================
//
//void MECAT2PAF(ifstream& input, char* filename, ifstream& index)
//{
//	map<uint32_t, std::string> names;
//	mecatidx (index, names);
//
//	int maxt = 1;
//#pragma omp parallel
//	{
//		maxt = omp_get_num_threads();
//	}
//
//	uint64_t numoverlap = std::count(std::istreambuf_iterator<char>(input), std::istreambuf_iterator<char>(), '\n');
//	input.seekg(0, std::ios_base::beg);
//
//	vector<std::string> entries;
//	vector<std::stringstream> local(maxt);      
//
//	/* read input file */
//	if(input)
//		for (int i = 0; i < numoverlap; ++i)
//		{
//			std::string line;
//			std::getline(input, line);
//			entries.push_back(line);
//		}
//	input.close();
//
//	/* transform MECAT output in PAF format */
//#pragma omp parallel for
//	for(uint64_t i = 0; i < numoverlap; i++) 
//	{
//		std::stringstream linestream(entries[i]);
//		int ithread = omp_get_thread_num();
//
//		/* MECAT format: cid, rid, score, id, cstr, cstart, cend, clen, rstr, rstart, rend, rlen */
//		std::vector<std::string> v = split (entries[i], '\t');
//
//		/* mecat idx to nametag translation */
//		std::string nameV = idx2read (stoi(v[0]), names);
//		std::string nameH = idx2read (stoi(v[1]), names);
//
//		std::string& ident = v[2];
//		std::string& begpV = v[5];
//		std::string& endpV = v[6];
//		std::string& lengV = v[7];
//		std::string& isRev = v[8];
//		std::string& begpH = v[9];
//		std::string& endpH = v[10];
//		std::string& lengH = v[11];
//
//		/* change strand formatting */
//		if(isRev == "0") isRev = "+";         
//			else isRev = "-";
//
//		/* compute overlap length if missing (begpV, endpV, lenV, begpH, endpH, lenH) */
//		int ovlen = estimate (stoi(begpV), stoi(endpV), stoi(lengV), stoi(begpH), stoi(endpH), stoi(lengH));
//		/* If alignment is missing I estimate it as (ident) *ovlen */
//		//int score = floor((stod(ident)*ovlen) / 100);
//		int score = floor((stod(ident)*ovlen) / 100);
//
//		/* GGGG: I might need to translate back idx to original names ---> YES, I NEED THE ORIGINAL NAMES. */
//		local[ithread] << nameV << "\t" << lengV << "\t" << begpV << "\t" << endpV << "\t" << isRev 
//			<< "\t" << nameH << "\t" << lengH << "\t" << begpH << "\t" << endpH << "\t" << score 
//				<< "\t" << ovlen << "\t255" << endl;
//	}
//
//	/* write to a new file */
//	int64_t * bytes = new int64_t[maxt];
//	for(int i = 0; i < maxt; ++i)
//	{
//		local[i].seekg(0, ios::end);
//		bytes[i] = local[i].tellg();
//		local[i].seekg(0, ios::beg);
//	}
//	int64_t bytestotal = std::accumulate(bytes, bytes + maxt, static_cast<int64_t>(0));
//
//	std::ofstream output(filename, std::ios::binary | std::ios::app);
//#ifdef PRINT
//	cout << "Creating or appending to output file with " << (double)bytestotal/(double)(1024 * 1024) << " MB" << endl;
//#endif
//	output.seekp(bytestotal - 1);
//	/* this will likely create a sparse file so the actual disks won't spin yet */
//	output.write("", 1); 
//	output.close();
//
//	#pragma omp parallel
//	{
//		int ithread = omp_get_thread_num(); 
//
//		FILE *ffinal;
//		/* then everyone fills it */
//		if ((ffinal = fopen(filename, "rb+")) == NULL) 
//		{
//			fprintf(stderr, "File %s failed to open at thread %d\n", filename, ithread);
//		}
//		int64_t bytesuntil = std::accumulate(bytes, bytes + ithread, static_cast<int64_t>(0));
//		fseek (ffinal , bytesuntil , SEEK_SET);
//		std::string text = local[ithread].str();
//		fwrite(text.c_str(),1, bytes[ithread], ffinal);
//		fflush(ffinal);
//		fclose(ffinal);
//	}
//	delete [] bytes;
//}
//
////=======================================================================
//// 
//// BLASR statistics
//// 
////=======================================================================
//
//void BLASR2PAF(ifstream& input, char* filename)
//{
//	int maxt = 1;
//#pragma omp parallel
//	{
//		maxt = omp_get_num_threads();
//	}
//
//	uint64_t numoverlap = std::count(std::istreambuf_iterator<char>(input), std::istreambuf_iterator<char>(), '\n');
//	input.seekg(0, std::ios_base::beg);
//
//	vector<std::string> entries;
//	vector<std::stringstream> local(maxt);      
//
//	/* read input file */
//	if(input)
//		for (int i = 0; i < numoverlap; ++i)
//		{
//			std::string line;
//			std::getline(input, line);
//			entries.push_back(line);
//		}
//	input.close();
//
//	/* transform BLASR output in PAF format */
//#pragma omp parallel for
//	for(uint64_t i = 0; i < numoverlap; i++) 
//	{
//		int ithread = omp_get_thread_num();
//		/* BLASR format: cname, rname, score, id, cstr, cstart, cend, clen, rstr, rstart, rend, rlen, qv */
//		std::vector<std::string> v = split (entries[i], ' ');
//		/* improve readability */
//		std::string& nameV = v[0];
//		std::string& nameH = v[1];
//		std::string& score = v[2];
//		std::string& strnV = v[4];
//		std::string& begpV = v[5];
//		std::string& endpV = v[6];
//		std::string& lengV = v[7];
//		std::string& strnH = v[8];
//		std::string& begpH = v[9];
//		std::string& endpH = v[10];
//		std::string& lengH = v[11];
//		std::string& mapQV = v[12];
//
//		/* change strand formatting */
//		std::string isRev;
//		if(strnH == strnV) isRev = "+";         
//			else isRev = "-";
//
//		// GGGG: BLSR scores are negatives? Dig into this.
//		score.erase(0, 1); 
//		/* compute overlap length if missing (begpV, endpV, lenV, begpH, endpH, lenH) */
//		int ovlen = estimate (stoi(begpV), stoi(endpV), stoi(lengV), stoi(begpH), stoi(endpH), stoi(lengH));
//
//		local[ithread] << nameV << "\t" << lengV << "\t" << begpV << "\t" << endpV << "\t" << isRev 
//			<< "\t" << nameH << "\t" << lengH << "\t" << begpH << "\t" << endpH << "\t" << score 
//				<< "\t" << ovlen << "\t" << mapQV << endl;
//	}
//
//	/* write to a new file */
//	int64_t * bytes = new int64_t[maxt];
//	for(int i = 0; i < maxt; ++i)
//	{
//		local[i].seekg(0, ios::end);
//		bytes[i] = local[i].tellg();
//		local[i].seekg(0, ios::beg);
//	}
//	int64_t bytestotal = std::accumulate(bytes, bytes + maxt, static_cast<int64_t>(0));
//
//	std::ofstream output(filename, std::ios::binary | std::ios::app);
//#ifdef PRINT
//	cout << "Creating or appending to output file with " << (double)bytestotal/(double)(1024 * 1024) << " MB" << endl;
//#endif
//	output.seekp(bytestotal - 1);
//	/* this will likely create a sparse file so the actual disks won't spin yet */
//	output.write("", 1); 
//	output.close();
//
//	#pragma omp parallel
//	{
//		int ithread = omp_get_thread_num(); 
//
//		FILE *ffinal;
//		/* then everyone fills it */
//		if ((ffinal = fopen(filename, "rb+")) == NULL) 
//		{
//			fprintf(stderr, "File %s failed to open at thread %d\n", filename, ithread);
//		}
//		int64_t bytesuntil = std::accumulate(bytes, bytes + ithread, static_cast<int64_t>(0));
//		fseek (ffinal , bytesuntil , SEEK_SET);
//		std::string text = local[ithread].str();
//		fwrite(text.c_str(),1, bytes[ithread], ffinal);
//		fflush(ffinal);
//		fclose(ffinal);
//	}
//	delete [] bytes;
//}
//
////=======================================================================
//// 
//// DALIGNER statistics
//// 
////=======================================================================
//
//void DALIGNER2PAF(ifstream& input, char* filename)
//{
//	int maxt = 1;
//#pragma omp parallel
//	{
//		maxt = omp_get_num_threads();
//	}
//
//	uint64_t numoverlap = std::count(std::istreambuf_iterator<char>(input), std::istreambuf_iterator<char>(), '\n');
//	input.seekg(0, std::ios_base::beg);
//
//	vector<std::string> entries;
//	vector<std::stringstream> local(maxt);      
//
//	/* read input file */
//	if(input)
//		for (int i = 0; i < numoverlap; ++i)
//		{
//			std::string line;
//			std::getline(input, line);
//			entries.push_back(line);
//		}
//	input.close();
//
//	/* transform DALIGNER output in PAF format */
//#pragma omp parallel for
//	for(uint64_t i = 0; i < numoverlap; i++) 
//	{
//		int ithread = omp_get_thread_num();
//		/* DALIGNER format: cname, rname, rev, cstart, cend, clen, rstart, rend, rlen */
//		std::vector<std::string> v = split (entries[i], ' ');
//		/* improve readability */
//		std::string& nameV = v[0];
//		std::string& nameH = v[1];
//		std::string& isRev = v[2];
//		std::string& begpV = v[3];
//		std::string& endpV = v[4];
//		std::string& lengV = v[5];
//		std::string& begpH = v[6];
//		std::string& endpH = v[7];
//		std::string& lengH = v[8];
//
//		/* change strand formatting */
//		if(isRev == "n") isRev = "+";         
//			else isRev = "-";
//
//		/* compute overlap length if missing (begpV, endpV, lenV, begpH, endpH, lenH) */
//		int ovlen = estimate (stoi(begpV), stoi(endpV), stoi(lengV), stoi(begpH), stoi(endpH), stoi(lengH));
//
//		/* GGGG: If alignment is missing I estimate it as % of the overlap length and I determine that % using the error rate */
//		// GGGG: Error rate is now hard-coded, need to be an input parameter
//
//		float error = 0.15;
//		float identity = (1-error)*(1-error);
//		int score = floor(identity*ovlen);
//
//		local[ithread] << nameV << "\t" << lengV << "\t" << begpV << "\t" << endpV << "\t" << isRev 
//			<< "\t" << nameH << "\t" << lengH << "\t" << begpH << "\t" << endpH << "\t" << score 
//				<< "\t" << ovlen << "\t255" << endl;
//	}
//
//	/* write to a new file */
//	int64_t * bytes = new int64_t[maxt];
//	for(int i = 0; i < maxt; ++i)
//	{
//		local[i].seekg(0, ios::end);
//		bytes[i] = local[i].tellg();
//		local[i].seekg(0, ios::beg);
//	}
//	int64_t bytestotal = std::accumulate(bytes, bytes + maxt, static_cast<int64_t>(0));
//
//	std::ofstream output(filename, std::ios::binary | std::ios::app);
//#ifdef PRINT
//	cout << "Creating or appending to output file with " << (double)bytestotal/(double)(1024 * 1024) << " MB" << endl;
//#endif
//	output.seekp(bytestotal - 1);
//	/* this will likely create a sparse file so the actual disks won't spin yet */
//	output.write("", 1); 
//	output.close();
//
//	#pragma omp parallel
//	{
//		int ithread = omp_get_thread_num(); 
//
//		FILE *ffinal;
//		/* then everyone fills it */
//		if ((ffinal = fopen(filename, "rb+")) == NULL) 
//		{
//			fprintf(stderr, "File %s failed to open at thread %d\n", filename, ithread);
//		}
//		int64_t bytesuntil = std::accumulate(bytes, bytes + ithread, static_cast<int64_t>(0));
//		fseek (ffinal , bytesuntil , SEEK_SET);
//		std::string text = local[ithread].str();
//		fwrite(text.c_str(),1, bytes[ithread], ffinal);
//		fflush(ffinal);
//		fclose(ffinal);
//	}
//	delete [] bytes;
//}
