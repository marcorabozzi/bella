//==================================================================
// Title:  LOGAN: X-Drop Adaptive Banded Alignment
// Author: G. Guidi, E. Younis
// Date:   29 April 2019
//==================================================================

#include <vector>
#include <iostream>
#include <string.h>
#include <algorithm>
#include <chrono>
#include <stdlib.h>
#include <inttypes.h>
#include <assert.h>
#include <iterator>
#include <x86intrin.h>
#include "logan_xa_affine_int8.h"
//#include "logan_ga.h"

//======================================================================================
// SEQUENCE GENERATION (source: https://github.com/ocxtal/libgaba)
//======================================================================================

char random_base(void)
{
	char const table[4] = {'A', 'C', 'G', 'T'};
	return(table[rand() % 4]);
}

void generate_random_sequence(std::string& seq, unsigned int const& len)
{
	for(int i = 0; i < len; i++)
		seq.append(1, random_base());
}

std::string generate_mutated_sequence(std::string& seq, unsigned int const& len, 
	unsigned short const& bw, double const& indels, double const& substs)
{
	int i, j, wave = 0;	// wave is q-coordinate of the alignment path
	std::string mutated;

	for(i = 0, j = 0; i < len; i++)
	{
		if(((double)rand()/(double)RAND_MAX) < substs)
		{
			mutated.append(1, random_base());
			j++; // mismatch
		}
		else if(((double)rand()/(double)RAND_MAX) < indels)
		{
			if(rand() & 0x01 && wave > (-bw + 1))
			{
				char tmp = (j < len) ? seq[j++] : random_base();
				mutated.append(1, tmp);
				j++; wave--; // deletion
			}
			else if(wave < (bw-2))
			{
				mutated.append(1, random_base());
				wave++; // insertion
			}
			else
			{
				char tmp = (j < len) ? seq[j++] : random_base();
				mutated.append(1, tmp);
			}
		}
		else
		{
			char tmp = (j < len) ? seq[j++] : random_base();
			mutated.append(1, tmp);
		}
	}
	return mutated;
}

//======================================================================================
// DEMO
//======================================================================================

int main(int argc, char const *argv[])
{
	/* Declarations */
	std::string seq1, seq2;

	/* Sequences length */
	unsigned int len1 = 10000;
	unsigned int len2 = 11500;

	/* Bandwidth (the alignment path of the input sequence and the result does not go out of the band) */
	unsigned short bw = 32;

	/* Error rate composition */
	double indels = 0.12; // indels probability
	double substs = 0.03; // substitution probability

	/* Penalties (LOGAN temporarily supports only linear gap penalty) */
	int8_t match    =  1;
	int8_t mismatch = -1;
	int8_t gap 	   = -1;

	/* Initialize scoring scheme */
	ScoringSchemeL penalties(match, mismatch, gap);

	/* Generate pair of sequences */
	generate_random_sequence(seq1, len1);
	seq2 = generate_mutated_sequence(seq1, len2, bw, indels, substs);

	/* x-drop value */
	unsigned short x = 100;

	/* seed/k-mer length */
	unsigned short k = 17;

	/* seed starting position on seq1, seed starting position on seq2, k-mer lenght */
	SeedL seed(0, 0, k);

	/* result.first = best score, result.second = exit score when (if) x-drop termination is satified */
	std::pair<int, int> result_x; 

	//======================================================================================
	// LOGAN (X-Drop Adaptive Banded Alignment)
	//======================================================================================

	result_x = LoganXDrop(seed, LOGAN_EXTEND_BOTH, seq1, seq2, penalties, x);
	std::cout << "Best score : " << result_x.first << "\tExit score : " << result_x.second << std::endl;

	//======================================================================================
	// LOGAN (Global Adaptive Banded Alignment)
	//======================================================================================

	//std::pair<short, short> result_g; 
	//result_g = LoganGlobal(seq1, seq2, penalties);
	//std::cout << "Best score : " << result_g.first << "\tGlobal score : " << result_g.second << std::endl;

	return 0;
}
