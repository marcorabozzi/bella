//==================================================================
// Title:  LOGAN: Adaptive Banded Global Alignment
// Author: G. Guidi, E. Younis
// Date:   30 April 2019
//==================================================================

#ifndef LOGANGA_H
#define LOGANGA_H

#include<vector>
#include<iostream>
#include<omp.h>
#include<algorithm>
#include<inttypes.h>
#include<assert.h>
#include<iterator>
#include<x86intrin.h>
#include"utils.h"
#include"simd_utils.h"
#include"score.h"

//======================================================================================
// GLOBAL ADAPTIVE BANDED ALIGNMENT
//======================================================================================

std::pair<short, short>
LoganGlobal
(
	std::string const& targetSeg,
	std::string const& querySeg,
	ScoringSchemeL& scoringScheme
)
{
	unsigned int hlength = targetSeg.length() + 1;
	unsigned int vlength = querySeg.length()  + 1;

	if (hlength <= 1 || vlength <= 1)
		return std::make_pair(0, 0);

	// Convert from string to int array
	// This is the entire sequences
	short* queryh = new short[hlength];
	short* queryv = new short[vlength];
	std::copy(targetSeg.begin(), targetSeg.end(), queryh);
	std::copy(querySeg.begin(), querySeg.end(), queryv);

	short matchCost    = scoreMatch(scoringScheme   );
	short mismatchCost = scoreMismatch(scoringScheme);
	short gapCost      = scoreGap(scoringScheme     );

	vector_t vmatchCost    = set1_func (matchCost   );
	vector_t vmismatchCost = set1_func (mismatchCost);
	vector_t vgapCost      = set1_func (gapCost     );

	//======================================================================================
	// PHASE I (initial values load using dynamic programming)
	//======================================================================================

#ifdef DEBUG
	printf("Phase I\n");
#endif
	// we need one more space for the off-grid values and one more space for antiDiag2
	short phase1_data[LOGICALWIDTH + 2][LOGICALWIDTH + 2];

	// phase1_data initialization
	phase1_data[0][0] = 0;
	for (int i = 1; i < LOGICALWIDTH + 2; i++)
	{
		phase1_data[0][i] = -i;
		phase1_data[i][0] = -i;
	}

	// dynamic programming loop to fill phase1_data[][]
	for(int i = 1; i < LOGICALWIDTH + 2; i++)
		for(int j = 1; j < LOGICALWIDTH + 2; j++)
		{
			short onef = phase1_data[i-1][j-1];
			if(queryh[i-1] == queryv[j-1])
				onef += matchCost;
			else
				onef += mismatchCost;

			short twof = std::max(phase1_data[i-1][j], phase1_data[i][j-1]) + gapCost;
			phase1_data[i][j] = std::max(onef, twof);
		}

#ifdef DEBUG
	// print phase1_data[][]
	for(int i = 1; i < LOGICALWIDTH + 2; i++)
	{
		for(int j = 1; j < LOGICALWIDTH + 2; j++)
			std::cout << phase1_data[i][j] << '\t';
		std::cout << std::endl;
	}
#endif

	vector_union_t antiDiag1; 	// 16 (vector width) 16-bit integers
	vector_union_t antiDiag2; 	// 16 (vector width) 16-bit integers
	vector_union_t antiDiag3; 	// 16 (vector width) 16-bit integers

	vector_union_t vqueryh;
	vector_union_t vqueryv;

	// Initialize vqueryh and vqueryv
	for ( int i = 0; i < LOGICALWIDTH; ++i )
	{
		vqueryh.elem[i] = queryh[i + 1];
		vqueryv.elem[i] = queryv[LOGICALWIDTH - i];
	}

	vqueryh.elem[LOGICALWIDTH] = NINF;
	vqueryv.elem[LOGICALWIDTH] = NINF;

	// this should point to the next value to be loaded into vqueryh and vqueryv
	int hoffset = LOGICALWIDTH;
	int voffset = LOGICALWIDTH;

	// load phase1_data into antiDiag1 vector
	for (int i = 1; i <= LOGICALWIDTH; ++i)
		antiDiag1.elem[i-1] = phase1_data[i][LOGICALWIDTH - i + 1];
	antiDiag1.elem[LOGICALWIDTH] = NINF;

	// load phase1_data into antiDiag2 vector going myRIGHT (our arbitrary decision)
	// the first antiDiag3 computation is going myDOWN
	// shift to the right on updated vector 2 (This places the left-aligned vector 3 as a right-aligned vector 2)
	for (int i = 1; i <= LOGICALWIDTH; ++i)
		antiDiag2.elem[i] = phase1_data[i + 1][LOGICALWIDTH - i + 1];
	antiDiag2.elem[0] = NINF;

	// initialize antiDia3 to -inf
	antiDiag3.simd = set1_func (NINF);

	//======================================================================================
	// PHASE II (core vectorized computation)
	//======================================================================================

	short antiDiagNo = 1;
	short antiDiagBest = antiDiagNo * gapCost;
	short best = 0;

#ifdef DEBUG
	printf("Phase II\n");
#endif

	while(hoffset < hlength && voffset < vlength)
	{

#ifdef DEBUG
	printf("\n");
	print_vector_c(vqueryh.simd);
	print_vector_c(vqueryv.simd);
#endif

		// antiDiagBest initialization
		antiDiagNo++;
		antiDiagBest = antiDiagNo * gapCost;

		// antiDiag1F (final)
		// POST-IT: -1 for a match and 0 for a mismatch
		vector_t m = cmpeq_func (vqueryh.simd, vqueryv.simd);
		m = blendv_func (vmismatchCost, vmatchCost, m);
		vector_t antiDiag1F = add_func (m, antiDiag1.simd);

	#ifdef DEBUG
		printf("antiDiag1: ");
		print_vector_d(antiDiag1.simd);
		printf("antiDiag1F: ");
		print_vector_d(antiDiag1F);
	#endif

		// antiDiag2S (shift)
		vector_union_t antiDiag2S = leftShift (antiDiag2);
	#ifdef DEBUG
		printf("antiDiag2S: ");
		print_vector_d(antiDiag2S.simd);
	#endif
		// antiDiag2M (pairwise max)
		vector_t antiDiag2M = max_func (antiDiag2S.simd, antiDiag2.simd);
	#ifdef DEBUG
		printf("antiDiag2M: ");
		print_vector_d(antiDiag2M);
	#endif
		// antiDiag2F (final)
		vector_t antiDiag2F = add_func (antiDiag2M, vgapCost);
	#ifdef DEBUG
		printf("antiDiag2F: ");
		print_vector_d(antiDiag2F);
	#endif
	#ifdef DEBUG
		printf("antiDiag2: ");
		print_vector_d(antiDiag2.simd);
	#endif
		// Compute antiDiag3
		antiDiag3.simd = max_func (antiDiag1F, antiDiag2F);
		// we need to have always antiDiag3 left-aligned
		antiDiag3.elem[LOGICALWIDTH] = NINF;
	#ifdef DEBUG
		printf("antiDiag3: ");
		print_vector_d(antiDiag3.simd);
	#endif

		// update best
		antiDiagBest = *std::max_element(antiDiag3.elem, antiDiag3.elem + VECTORWIDTH);
		best = (best > antiDiagBest) ? best : antiDiagBest;

		// antiDiag swap, offset updates, and new base load
		// TODO : optimize this
		int maxpos, max = 0;
		for(int i = 0; i < VECTORWIDTH; ++i)
			if(antiDiag3.elem[i] > max)
			{
				maxpos = i;
				max = antiDiag3.elem[i];
			}

		//if(antiDiag3.elem[MIDDLE] < antiDiag3.elem[MIDDLE + 1])
		if(maxpos > MIDDLE)
			moveRight (antiDiag1, antiDiag2, antiDiag3, hoffset, voffset, vqueryh, vqueryv, queryh, queryv);
		else
			moveDown (antiDiag1, antiDiag2, antiDiag3, hoffset, voffset, vqueryh, vqueryv, queryh, queryv);
	}

	//======================================================================================
	// PHASE III (we are one edge)
	//======================================================================================

	int dir = hoffset >= hlength ? myDOWN : myRIGHT;

#ifdef DEBUG
	printf("Phase III\n");
#endif

	while(hoffset < hlength || voffset < vlength)
	{

	#ifdef DEBUG
		printf("\n");
		print_vector_c(vqueryh.simd);
		print_vector_c(vqueryv.simd);
	#endif
		// antiDiagBest initialization
		antiDiagNo++;
		antiDiagBest = antiDiagNo * gapCost;

		// antiDiag1F (final)
		// POST-IT: -1 for a match and 0 for a mismatch
		vector_t m = cmpeq_func (vqueryh.simd, vqueryv.simd);
		m = blendv_func (vmismatchCost, vmatchCost, m);
		vector_t antiDiag1F = add_func (m, antiDiag1.simd);

	#ifdef DEBUG
		printf("antiDiag1: ");
		print_vector_d(antiDiag1.simd);
		printf("antiDiag1F: ");
		print_vector_d(antiDiag1F);
	#endif

		// antiDiag2S (shift)
		vector_union_t antiDiag2S = leftShift (antiDiag2);
	#ifdef DEBUG
		printf("antiDiag2S: ");
		print_vector_d(antiDiag2S.simd);
	#endif
		// antiDiag2M (pairwise max)
		vector_t antiDiag2M = max_func (antiDiag2S.simd, antiDiag2.simd);
	#ifdef DEBUG
		printf("antiDiag2M: ");
		print_vector_d(antiDiag2M);
	#endif
		// antiDiag2F (final)
		vector_t antiDiag2F = add_func (antiDiag2M, vgapCost);
	#ifdef DEBUG
		printf("antiDiag2F: ");
		print_vector_d(antiDiag2F);
	#endif
	#ifdef DEBUG
		printf("antiDiag2: ");
		print_vector_d(antiDiag2.simd);
	#endif
		// Compute antiDiag3
		antiDiag3.simd = max_func (antiDiag1F, antiDiag2F);
		// we need to have always antiDiag3 left-aligned
		antiDiag3.elem[LOGICALWIDTH] = NINF;
	#ifdef DEBUG
		printf("antiDiag3: ");
		print_vector_d(antiDiag3.simd);
	#endif

		// update best
		antiDiagBest = *std::max_element(antiDiag3.elem, antiDiag3.elem + VECTORWIDTH);
		best = (best > antiDiagBest) ? best : antiDiagBest;

		// antiDiag swap, offset updates, and new base load
		if (dir == myRIGHT)
			moveRight (antiDiag1, antiDiag2, antiDiag3, hoffset, voffset, vqueryh, vqueryv, queryh, queryv);
		else
			moveDown (antiDiag1, antiDiag2, antiDiag3, hoffset, voffset, vqueryh, vqueryv, queryh, queryv);
	}

	//======================================================================================
	// PHASE IV (reaching end of sequences)
	//======================================================================================

#ifdef DEBUG
	printf("Phase IV\n");
#endif
	for (int i = 0; i < (LOGICALWIDTH - 3); i++)
	{
		// antiDiag1F (final)
		// POST-IT: -1 for a match and 0 for a mismatch
		vector_t m = cmpeq_func (vqueryh.simd, vqueryv.simd);
		m = blendv_func (vmismatchCost, vmatchCost, m);
		vector_t antiDiag1F = add_func (m, antiDiag1.simd);

	#ifdef DEBUG
		printf("\n");
		printf("antiDiag1: ");
		print_vector_d(antiDiag1.simd);
		printf("antiDiag1F: ");
		print_vector_d(antiDiag1F);
	#endif

		// antiDiag2S (shift)
		vector_union_t antiDiag2S = leftShift (antiDiag2);
	#ifdef DEBUG
		printf("antiDiag2S: ");
		print_vector_d(antiDiag2S.simd);
	#endif
		// antiDiag2M (pairwise max)
		vector_t antiDiag2M = max_func (antiDiag2S.simd, antiDiag2.simd);
	#ifdef DEBUG
		printf("antiDiag2M: ");
		print_vector_d(antiDiag2M);
	#endif

		// antiDiag2F (final)
		vector_t antiDiag2F = add_func (antiDiag2M, vgapCost);
	#ifdef DEBUG
		printf("antiDiag2F: ");
		print_vector_d(antiDiag2F);
	#endif
		// Compute antiDiag3
		antiDiag3.simd = max_func (antiDiag1F, antiDiag2F);
		// we need to have always antiDiag3 left-aligned
		antiDiag3.elem[LOGICALWIDTH] = NINF;
	#ifdef DEBUG
		printf("antiDiag2: ");
		print_vector_d(antiDiag2.simd);
	#endif
	#ifdef DEBUG
		printf("antiDiag3: ");
		print_vector_d(antiDiag3.simd);
	#endif

		// update best
		antiDiagBest = *std::max_element(antiDiag3.elem, antiDiag3.elem + VECTORWIDTH);
		best = (best > antiDiagBest) ? best : antiDiagBest;

		// antiDiag swap, offset updates, and new base load
		short nextDir = dir ^ 1;
		// antiDiag swap, offset updates, and new base load
		if (nextDir == myRIGHT)
			moveRight (antiDiag1, antiDiag2, antiDiag3, hoffset, voffset, vqueryh, vqueryv, queryh, queryv);
		else
			moveDown (antiDiag1, antiDiag2, antiDiag3, hoffset, voffset, vqueryh, vqueryv, queryh, queryv);
		// direction update
		dir = nextDir;
	}

	delete [] queryh;
	delete [] queryv;

	return std::make_pair(best, antiDiagBest);
}

#endif