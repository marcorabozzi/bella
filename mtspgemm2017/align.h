#ifndef _ALIGNMENT_H_
#define _ALIGNMENT_H_

#include <seqan/sequence.h>
#include <seqan/align.h>
#include <seqan/score.h>
#include <seqan/modifier.h>
#include <seqan/seeds.h>
#include "../logan/src/simd/score.h"
#include "../logan/src/simd/simd_utils.h"
#include "../logan/src/simd/logan_xa.h"
#include "common.h"
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

using namespace seqan;
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

seqAnResult alignSeqAn(const std::string & row, const std::string & col, int rlen, int i, int j, int xdrop, int kmer_len) {

	ScoringSchemeL scoringSchemeLogan(1,-1,-1);
	Score<int, Simple> scoringScheme(1,-1,-1);

	Dna5String seqH(row); 
	Dna5String seqV(col); 
	Dna5String seedH;
	Dna5String seedV;
	string strand;
	std::pair<int, int> temp;
	int longestExtensionTemp;
	seqAnResult longestExtensionScore;
	loganResult result;

	TSeed seed(i, j, i+kmer_len, j+kmer_len);
	TSeedLogan seedLogan(i, j, kmer_len);
	seedH = infix(seqH, beginPositionH(seed), endPositionH(seed));
	seedV = infix(seqV, beginPositionV(seed), endPositionV(seed));

	//std::string seedH = row.substr(i, kmer_len);
	//std::string seedV = col.substr(j, kmer_len);
	//std::string strand;

	/* we are reversing the "row", "col" is always on the forward strand */
	Dna5StringReverseComplement twin(seedH);

	if(twin == seedV)
	{
		strand = 'c';
		Dna5StringReverseComplement twinRead(seqH);

		std::string cpyrow = row;
		std::reverse(cpyrow.begin(), cpyrow.end());

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

		LoganSetBeginPositionH(seedLogan, i);
		LoganSetBeginPositionV(seedLogan, j);
		LoganSetEndPositionH(seedLogan, i + kmer_len);
		LoganSetEndPositionV(seedLogan, j + kmer_len);

		/* Perform match extension */
		//std::cout << seqH << std::endl;
		//std::cout << twinRead << std::endl;
		longestExtensionTemp = extendSeed(seed, twinRead, seqV, EXTEND_BOTH, scoringScheme, xdrop, kmer_len, GappedXDrop());
		std::cout << longestExtensionTemp << std::endl;

		//std::cout << row << std::endl;
		//std::cout << cpyrow << std::endl;

		temp = LoganXDrop(seedLogan, LOGAN_EXTEND_BOTH, cpyrow, col, scoringSchemeLogan, xdrop, kmer_len);
		std::cout << temp.second << std::endl;

	}
	 else
	{
		longestExtensionTemp = extendSeed(seed, seqH, seqV, EXTEND_BOTH, scoringScheme, xdrop, kmer_len, GappedXDrop());
		std::cout << longestExtensionTemp << std::endl;

		temp = LoganXDrop(seedLogan, LOGAN_EXTEND_BOTH, row, col, scoringSchemeLogan, xdrop, kmer_len);
		std::cout << temp.second << std::endl;

	} 

	longestExtensionScore.score = longestExtensionTemp;
	longestExtensionScore.seed = seed;
	longestExtensionScore.strand = strand;

	result.score = temp;
	result.seed = seedLogan;
	result.strand = strand;

	return longestExtensionScore;
}

#endif
