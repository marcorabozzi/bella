#ifndef _ALIGNMENT_H_
#define _ALIGNMENT_H_

#include "common.h"
#include "../logan/src/simd/score.h"
#include "../logan/src/simd/simd_utils.h"
#include "../logan/src/simd/logan_xa.h"
#include <omp.h>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <cstdlib>
#include <algorithm>
#include <ctype.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <dirent.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <stdint.h>

using namespace std;

double adaptiveSlope(double error)
{
	double p_mat = pow(1-error,2);  // match
	double p_mis = 1-p_mat;         // mismatch/gap
	double alpha = 1;               // match penalty
	double beta = 1;                // mismatch/gap penalty

	return alpha*p_mat - beta*p_mis;
}

bool toEnd(int colStart, int colEnd, int colLen, int rowStart, int rowEnd, int rowLen, int relaxMargin)
{
	int minLeft = min(colStart, rowStart);
	int minRight = min(colLen-colEnd, rowLen-rowEnd);

	 if(minLeft-relaxMargin <= 0)
		minLeft = 0;
	 if(minRight-relaxMargin <= 0)
		minRight = 0;

	 if((minLeft == 0 || minRight == 0))
		return true;
	else
		return false;
}

char revComplement(char n)
{
	switch(n)
	{
	case 'A':
		return 'T';
	case 'T':
		return 'A';
	case 'G':
		return 'C';
	case 'C':
		return 'G';
	}
	assert(false);
	return ' ';
}

/**
 * @brief alignSeqAn does the seed-and-extend alignment
 * @param row
 * @param col
 * @param rlen is the length of the row sequence
 * @param i is the starting position of the k-mer on the first read
 * @param j is the starting position of the k-mer on the second read
 * @param xdrop
 * @return alignment score and extended seed
 */
loganResult alignLogan(const std::string & row, const std::string & col, int rlen, int i, int j, int xdrop, int kmer_len) {

	ScoringSchemeL scoringScheme(1,-1,-1);

	std::string seedH = row.substr(i, kmer_len);
	std::string seedV = col.substr(j, kmer_len);
	std::string strand;

	std::pair<int, int> temp;
	loganResult result;
	TSeed seed(i, j, kmer_len);

	if(seedH != seedV)
	{
		strand = 'c';
		// reverse complement horizonatal sequence and update its seed position
		std::string cpyrow = row;

		std::transform(
			std::begin(cpyrow),
			std::end(cpyrow),
			std::begin(cpyrow),
		revComplement);
		i = rlen - i - kmer_len;

		setBeginPositionH(seed, i);
		setBeginPositionV(seed, j);
		setEndPositionH(seed, i+kmer_len);
		setEndPositionV(seed, j+kmer_len);

		#pragma omp critical
		{
			printf("s %d %d\n", i, j);
			printf("c %d %d %d %d %d %d\n", cpyrow.length(), col.length(), getBeginPositionH(seed), getEndPositionH(seed), getBeginPositionV(seed), getEndPositionV(seed));
		}

		// perform alignment
		temp = LoganXDrop(seed, LOGAN_EXTEND_RIGHT, cpyrow, col, scoringScheme, xdrop, kmer_len);

		#pragma omp critical
		{
			printf("c %d %d %d %d %d %d %d\n\n", temp.second, cpyrow.length(), col.length(), getBeginPositionH(seed), getEndPositionH(seed), getBeginPositionV(seed), getEndPositionV(seed));
		}
	}
	else
	{
		strand = 'n';
		// perform alignment
		#pragma omp critical
		{
			printf("s %d %d\n", i, j);
			printf("n %d %d %d %d %d %d\n", row.length(), col.length(), getBeginPositionH(seed), getEndPositionH(seed), getBeginPositionV(seed), getEndPositionV(seed));
		}
		temp = LoganXDrop(seed, LOGAN_EXTEND_RIGHT, row, col, scoringScheme, xdrop, kmer_len);

		#pragma omp critical
		{
			printf("n %d %d %d %d %d %d %d\n\n", temp.second, row.length(), col.length(), getBeginPositionH(seed), getEndPositionH(seed), getBeginPositionV(seed), getEndPositionV(seed));
		}
	
	} 

	result.score = temp;
	result.seed = seed;
	result.strand = strand;
	return result;
}

#endif
