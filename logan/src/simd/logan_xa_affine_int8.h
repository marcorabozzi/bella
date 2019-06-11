//==================================================================
// Title:  LOGAN: X-Drop Adaptive Banded Alignment
// Author: G. Guidi, E. Younis
// Date:   22 April 2019
//==================================================================

#ifndef LOGANXAA_H
#define LOGANXAA_H
//#define HISTORY
#include<vector>
#include<iostream>
#include<omp.h>
#include<algorithm>
#include<inttypes.h>
#include<assert.h>
#include<iterator>
#include<x86intrin.h>
#include"utils_affine_int8.h"
#include"simd_utils_affine_int8.h"
#include"score.h"

// #define DEBUGLogan

//======================================================================================
// X-DROP ADAPTIVE BANDED ALIGNMENT
//======================================================================================

std::pair<int, int>
LoganOneDirection
(
	SeedL & seed,
	std::string const& targetSeg,
	std::string const& querySeg,
	ScoringSchemeL& scoringScheme,
	int const &scoreDropOff
)
{
	unsigned int hlength = targetSeg.length() + 1;
	unsigned int vlength = querySeg.length()  + 1;
	int offset = 0;

	if (hlength <= 1 || vlength <= 1)
		return std::make_pair(0, 0);

	// Convert from string to int array
	// This is the entire sequences
	int8_t* queryh = new int8_t[hlength];
	int8_t* queryv = new int8_t[vlength];
	std::copy(targetSeg.begin(), targetSeg.end(), queryh);
	std::copy(querySeg.begin(), querySeg.end(), queryv);

	int8_t matchCost    = scoreMatch(scoringScheme   );
	int8_t mismatchCost = scoreMismatch(scoringScheme);
	int8_t gapCost      = scoreGap(scoringScheme     );
	int8_t gapOpening   = 0; // scoreGapOpening()

	vector_t vmatchCost    = set1_func (matchCost   );
	vector_t vmismatchCost = set1_func (mismatchCost);
	vector_t vgapCost      = set1_func (gapCost     );
	vector_t vgapopening   = set1_func (gapOpening  );

	//======================================================================================
	// PHASE I (initial values load using dynamic programming)
	//======================================================================================

#ifdef DEBUGLogan
	printf("Phase I\n");
#endif
	// we need one more space for the off-grid values and one more space for antiDiag2
	// TODO: worry about overflow in phase 1 with int8_t
	int8_t phase1_data[LOGICALWIDTH + 2][LOGICALWIDTH + 2];
	int8_t phase1_gaphistup[ LOGICALWIDTH + 2 ][ LOGICALWIDTH + 2 ];
	int8_t phase1_gaphistleft[ LOGICALWIDTH + 2 ][ LOGICALWIDTH + 2 ];

	// phase1_data initialization
	phase1_data[0][0] = 0;
	for (int i = 1; i < LOGICALWIDTH + 2; i++)
	{
		phase1_data[0][i] = -i + gapOpening;
		phase1_data[i][0] = -i + gapOpening;
	}

	phase1_gaphistleft[0][0] = gapOpening;
	phase1_gaphistup[0][0]   = gapOpening;

	for (int i = 1; i < LOGICALWIDTH + 2; i++)
	{
		phase1_gaphistleft[0][i] = 0;
		phase1_gaphistup[0][i]   = gapOpening;

		phase1_gaphistleft[i][0] = gapOpening;
		phase1_gaphistup[i][0]   = 0;
	}

	// dynamic programming loop to fill phase1_data[][]
	for(int i = 1; i < LOGICALWIDTH + 2; i++)
	{
		for(int j = 1; j < LOGICALWIDTH + 2; j++)
		{
			int8_t onef = phase1_data[i-1][j-1];
			if(queryh[i-1] == queryv[j-1])
				onef += matchCost;
			else
				onef += mismatchCost;

			int8_t twou = phase1_data[i-1][j] + gapCost + phase1_gaphistup[i-1][j];
			int8_t twol = phase1_data[i][j-1] + gapCost + phase1_gaphistleft[i][j-1];
			int8_t twof = std::max( twou, twol );

			phase1_data[i][j] = std::max(onef, twof);

			if ( phase1_data[i][j] == onef )
			{
				phase1_gaphistup[i][j]   = gapOpening;
				phase1_gaphistleft[i][j] = gapOpening;
			}
			else if ( phase1_data[i][j] == twou )
			{
				phase1_gaphistup[i][j]   = 0;
				phase1_gaphistleft[i][j] = gapOpening;
			}
			else
			{
				phase1_gaphistup[i][j]   = gapOpening;
				phase1_gaphistleft[i][j] = 0;
			}
		}
	}

#ifdef DEBUGLogan
	// print phase1_data[][]
	for(int i = 1; i < LOGICALWIDTH + 2; i++)
	{
		for(int j = 1; j < LOGICALWIDTH + 2; j++)
			std::cout << phase1_data[i][j] << '\t';
		std::cout << std::endl;
	}
#endif

	vector_union_t antiDiag1; 	// 32 (vector width) 8-bit integers
	vector_union_t antiDiag2; 	// 32 (vector width) 8-bit integers
	vector_union_t antiDiag3; 	// 32 (vector width) 8-bit integers

	vector_union_t gaphist_up;
	vector_union_t gaphist_left;

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

	// initialize antiDiag3 to -inf
	antiDiag3.simd = set1_func (NINF);

	// initialize gaphist
	// GGGG: why 1?
	for (int i = 0; i <= LOGICALWIDTH; ++i)
	{
		gaphist_left.elem[i] = phase1_gaphistleft[i + 1][LOGICALWIDTH - i + 1];
	}
	//gaphist_left.elem[0] = NINF;

	for (int i = 0; i <= LOGICALWIDTH; ++i)
	{
		gaphist_up.elem[i] = phase1_gaphistup[i + 1][LOGICALWIDTH - i + 1];
	}
	//gaphist_up.elem[0] = NINF;

	//======================================================================================
	// PHASE II (core vectorized computation)
	//======================================================================================

	int antiDiagNo   = 1;
	int best   = 0;
	int dir;
	int8_t antiDiagBest = (int8_t)antiDiagNo * gapCost;

#ifdef DEBUGLogan
	printf("Phase II\n");
#endif

	while(hoffset < hlength && voffset < vlength)
	{
	printf("\n");
	print_vector_d(gaphist_up.simd);
	print_vector_d(gaphist_left.simd);
#ifdef DEBUGLogan
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
		m = blend_func (vmismatchCost, vmatchCost, m);
		vector_t antiDiag1F = add_func (m, antiDiag1.simd);

	#ifdef DEBUGLogan
		printf("antiDiag1: ");
		print_vector_d(antiDiag1.simd);
		printf("antiDiag1F: ");
		print_vector_d(antiDiag1F);
	#endif

	// antiDiag2U (Up)
		vector_union_t antiDiag2U = leftShift (antiDiag2);
		antiDiag2U.simd = add_func (antiDiag2U.simd, vgapCost);
		antiDiag2U.simd = add_func (antiDiag2U.simd, gaphist_up.simd);
	#ifdef DEBUGLogan
		printf("antiDiag2U: ");
		print_vector_d(antiDiag2U.simd);
	#endif

	// antiDiag2L (Left)
		vector_union_t antiDiag2L = antiDiag2;
		antiDiag2L.elem[VECTORWIDTH - 1] = NINF;
		antiDiag2L.simd = add_func (antiDiag2L.simd, vgapCost);
		antiDiag2L.simd = add_func (antiDiag2L.simd, gaphist_left.simd);
	#ifdef DEBUGLogan
		printf("antiDiag2L: ");
		print_vector_d(antiDiag2L.simd);
	#endif

		// antiDiag2M (pairwise max)
		vector_t antiDiag2M = max_func (antiDiag2L.simd, antiDiag2U.simd);
	#ifdef DEBUGLogan
		printf("antiDiag2M: ");
		print_vector_d(antiDiag2M);
	#endif
	#ifdef DEBUGLogan
		printf("antiDiag2: ");
		print_vector_d(antiDiag2.simd);
	#endif
		// Compute antiDiag3
		antiDiag3.simd = max_func (antiDiag1F, antiDiag2M);
		// we need to have always antiDiag3 left-aligned
		antiDiag3.elem[LOGICALWIDTH] = NINF;
	#ifdef DEBUGLogan
		printf("antiDiag3: ");
		print_vector_d(antiDiag3.simd);
	#endif

		// Gap History Computation
	#ifdef HISTORY
		printf("antiDiag2M: ");
		print_vector_d(antiDiag2M);
		printf("antiDiag2U.simd: ");
		print_vector_d(antiDiag2U.simd);
	#endif
		vector_t twogap 	= cmpeq_func	(antiDiag2M, antiDiag2U.simd);
	#ifdef HISTORY
		printf("twogap: ");
		print_vector_d(twogap);
	#endif
		gaphist_up.simd 	= blend_func	(_mm256_setzero_si256(), vgapopening, twogap); 	// TODO: Check values
	#ifdef HISTORY
		printf("gaphist_up.simd: ");
		print_vector_d(gaphist_up.simd);
	#endif
		gaphist_left.simd 	= blend_func 	(vgapopening, _mm256_setzero_si256(), twogap); 	// TODO: Check values
	#ifdef HISTORY
		printf("gaphist_left.simd: ");
		print_vector_d(gaphist_left.simd);
	#endif
		vector_t threegap 	= cmpeq_func 	(antiDiag3.simd, antiDiag1F);
	#ifdef HISTORY
		printf("threegap: ");
		print_vector_d(threegap);
	#endif
		gaphist_up.simd 	= blend_func 	(vgapopening, gaphist_up.simd, threegap); 	// TODO: Check values
	#ifdef HISTORY
		printf("gaphist_up.simd: ");
		print_vector_d(gaphist_up.simd);
	#endif
		gaphist_left.simd 	= blend_func 	(vgapopening, gaphist_left.simd, threegap); // TODO: Check values
	#ifdef HISTORY
		printf("gaphist_left.simd: ");
		print_vector_d(gaphist_left.simd);
	#endif
		// TODO: x-drop termination
		// Note: Don't need to check x drop every time
		antiDiagBest = *std::max_element(antiDiag3.elem, antiDiag3.elem + VECTORWIDTH);
		if(antiDiagBest + offset < best - scoreDropOff)
		{
			setEndPositionH(seed, hoffset);
			setEndPositionV(seed, voffset);

			delete [] queryh;
			delete [] queryv;

			return std::make_pair(best, antiDiagBest + offset);
		}

		if(antiDiagBest > INTLIMIT)
		{
			int8_t antiDiagMin = *std::min_element(antiDiag3.elem, antiDiag3.elem + LOGICALWIDTH);
			antiDiag2.simd = sub_func(antiDiag2.simd, set1_func(antiDiagMin));
			antiDiag3.simd = sub_func(antiDiag3.simd, set1_func(antiDiagMin));
			// std::cout << "offset prev " << offset << std::endl;
			offset += (int)antiDiagMin;
			// std::cout << "antiDiagMin " << (int)antiDiagMin << std::endl;
			// std::cout << "offset post " << offset << std::endl;
		}

		// update best
		best = (best > antiDiagBest + offset) ? best : antiDiagBest + offset;

		// antiDiag swap, offset updates, and new base load
		// TODO : optimize this
		int maxpos, max = 0;
		for(int i = 0; i < VECTORWIDTH; ++i)
			if(antiDiag3.elem[i] > max)
			{
				maxpos = i;
				max = antiDiag3.elem[i];
			}

		if(maxpos > MIDDLE)
		//if(antiDiag3.elem[MIDDLE] < antiDiag3.elem[MIDDLE + 1])
		{
			#ifdef DEBUGLogan
			printf("myRIGHT\n");
			#endif
			moveRight (antiDiag1, antiDiag2, antiDiag3, hoffset, voffset, vqueryh, vqueryv, queryh, queryv);
			dir = myRIGHT;
		}
		else
		{
			#ifdef DEBUGLogan
			printf("myDOWN\n");
			#endif
			moveDown (antiDiag1, antiDiag2, antiDiag3, hoffset, voffset, vqueryh, vqueryv, queryh, queryv);
			dir = myDOWN;
		}
	}

	// at least one sequence has achieved its max length
	int hextension = hoffset - 1;
	int vextension = voffset - 1;

	//======================================================================================
	// PHASE III (we are one edge)
	//======================================================================================

	// REMOVED (NOT NEEDED IN BELLA)

	//======================================================================================
	// PHASE IV (reaching end of sequences)
	//======================================================================================

#ifdef DEBUGLogan
	printf("Phase IV\n");
#endif
	for (int i = 0; i < (LOGICALWIDTH - 3); i++)
	{

#ifdef DEBUGLogan
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
		m = blend_func (vmismatchCost, vmatchCost, m);
		vector_t antiDiag1F = add_func (m, antiDiag1.simd);

	#ifdef DEBUGLogan
		printf("antiDiag1: ");
		print_vector_d(antiDiag1.simd);
		printf("antiDiag1F: ");
		print_vector_d(antiDiag1F);
	#endif

		// antiDiag2U (Up)
		vector_union_t antiDiag2U = leftShift (antiDiag2);
		antiDiag2U.simd = add_func (antiDiag2U.simd, vgapCost);
		antiDiag2U.simd = add_func (antiDiag2U.simd, gaphist_up.simd);
	#ifdef DEBUGLogan
		printf("antiDiag2U: ");
		print_vector_d(antiDiag2U.simd);
	#endif

	// antiDiag2L (Left)
		vector_union_t antiDiag2L = antiDiag2;
		antiDiag2L.elem[VECTORWIDTH - 1] = NINF;
		antiDiag2L.simd = add_func (antiDiag2L.simd, vgapCost);
		antiDiag2L.simd = add_func (antiDiag2L.simd, gaphist_left.simd);
	#ifdef DEBUGLogan
		printf("antiDiag2L: ");
		print_vector_d(antiDiag2L.simd);
	#endif

		// antiDiag2M (pairwise max)
		vector_t antiDiag2M = max_func (antiDiag2L.simd, antiDiag2U.simd);
	#ifdef DEBUGLogan
		printf("antiDiag2M: ");
		print_vector_d(antiDiag2M);
	#endif
	#ifdef DEBUGLogan
		printf("antiDiag2: ");
		print_vector_d(antiDiag2.simd);
	#endif
		// Compute antiDiag3
		antiDiag3.simd = max_func (antiDiag1F, antiDiag2M);
		// we need to have always antiDiag3 left-aligned
		antiDiag3.elem[LOGICALWIDTH] = NINF;
	#ifdef DEBUGLogan
		printf("antiDiag3: ");
		print_vector_d(antiDiag3.simd);
	#endif

		// Gap History Computation
		vector_t twogap 	= cmpeq_func (antiDiag2M, antiDiag2U.simd);
		gaphist_up.simd 	= blend_func (_mm256_setzero_si256(), vgapopening, twogap); // TODO: Check values
		gaphist_left.simd 	= blend_func (vgapopening, _mm256_setzero_si256(), twogap); // TODO: Check values
		vector_t threegap 	= cmpeq_func (antiDiag3.simd, antiDiag1F);
		gaphist_up.simd 	= blend_func (vgapopening, gaphist_up.simd, threegap); 		// TODO: Check values
		gaphist_left.simd 	= blend_func (vgapopening, gaphist_left.simd, threegap); 	// TODO: Check values

		// x-drop termination
		antiDiagBest = *std::max_element(antiDiag3.elem, antiDiag3.elem + VECTORWIDTH);
		if(antiDiagBest + offset < best - scoreDropOff)
		{
			// Longest extension and update seed
			setEndPositionH(seed, hextension);
			setEndPositionV(seed, vextension);

			delete [] queryh;
			delete [] queryv;

			return std::make_pair(best, antiDiagBest + offset);
		}

		if (antiDiagBest > INTLIMIT)
		{
			int8_t antiDiagMin = *std::min_element(antiDiag3.elem, antiDiag3.elem + LOGICALWIDTH);
			antiDiag2.simd = sub_func (antiDiag2.simd, set1_func(antiDiagMin));
			antiDiag3.simd = sub_func (antiDiag3.simd, set1_func(antiDiagMin));
			// std::cout << "offset prev " << offset << std::endl;
			offset += (int)antiDiagMin;
			// std::cout << "antiDiagMin " << (int)antiDiagMin << std::endl;
			// std::cout << "offset post " << offset << std::endl;
		}

		// update best
		best = (best > antiDiagBest + offset) ? best : antiDiagBest + offset;

		// antiDiag swap, offset updates, and new base load
		short nextDir = dir ^ 1;
		// antiDiag swap, offset updates, and new base load
		if (nextDir == myRIGHT)
		{
		#ifdef DEBUGLogan
			printf("myRIGHT\n");
		#endif
			moveRight (antiDiag1, antiDiag2, antiDiag3, hoffset, voffset, vqueryh, vqueryv, queryh, queryv);
		}
		else
		{
		#ifdef DEBUGLogan
			printf("myDOWN\n");
		#endif
			moveDown (antiDiag1, antiDiag2, antiDiag3, hoffset, voffset, vqueryh, vqueryv, queryh, queryv);
		}
		// direction update
		dir = nextDir;
	}

	setEndPositionH(seed, hextension);
	setEndPositionV(seed, vextension);

	delete [] queryh;
	delete [] queryv;

	return std::make_pair(best, antiDiagBest + offset);
}

std::pair<int, int>
LoganXDrop
(
	SeedL& seed,
	ExtDirectionL direction,
	std::string const& target,
	std::string const& query,
	ScoringSchemeL& scoringScheme,
	int const &scoreDropOff
)
{
	// TODO: check scoring scheme correctness/input parameters
	if (direction == LOGAN_EXTEND_LEFT)
	{
		SeedL seed1 = seed;
		std::pair<int, int> extLeft;

		std::string targetPrefix = target.substr (0, getEndPositionH(seed));	// from read start til start seed (seed included)
		std::string queryPrefix  = query.substr  (0, getEndPositionV(seed));	// from read start til start seed (seed included)

		std::reverse (targetPrefix.begin(), targetPrefix.end());
		std::reverse (queryPrefix.begin(), queryPrefix.end());

		extLeft = LoganOneDirection (seed, targetPrefix, queryPrefix, scoringScheme, scoreDropOff);

		setBeginPositionH(seed, getEndPositionH(seed) - getEndPositionH(seed1));
		setBeginPositionV(seed, getEndPositionV(seed) - getEndPositionV(seed1));

		return extLeft;
	}
	else if (direction == LOGAN_EXTEND_RIGHT)
	{
		SeedL seed2 = seed;
		std::pair<int, int> extRight;

		std::string targetSuffix = target.substr (getBeginPositionH(seed), target.length()); 	// from end seed until the end (seed included)
		std::string querySuffix  = query.substr  (getBeginPositionV(seed), query.length());		// from end seed until the end (seed included)

		extRight = LoganOneDirection (seed, targetSuffix, querySuffix, scoringScheme, scoreDropOff);

		setEndPositionH(seed, getBeginPositionH(seed) + getEndPositionH(seed2));
		setEndPositionV(seed, getBeginPositionV(seed) + getEndPositionV(seed2));

		return extRight;
	}
	else
	{
		SeedL seed1 = seed;
		SeedL seed2 = seed;

		std::pair<int, int> extLeft;
		std::pair<int, int> extRight;

		std::string targetPrefix = target.substr (0, getEndPositionH(seed));	// from read start til start seed (seed not included)
		std::string queryPrefix  = query.substr  (0, getEndPositionV(seed));	// from read start til start seed (seed not included)

		std::reverse (targetPrefix.begin(), targetPrefix.end());
		std::reverse (queryPrefix.begin(), queryPrefix.end());

		std::cout << targetPrefix.size() << std::endl;
		std::cout << queryPrefix.size()  << std::endl;
		extLeft = LoganOneDirection (seed1, targetPrefix, queryPrefix, scoringScheme, scoreDropOff);

		std::cout << std::endl;
		std::cout << std::endl;

		std::string targetSuffix = target.substr (getEndPositionH(seed), target.length()); 	// from end seed until the end (seed included)
		std::string querySuffix  = query.substr  (getEndPositionV(seed), query.length());	// from end seed until the end (seed included)

		std::cout << targetSuffix.size() << std::endl;
		std::cout << querySuffix.size()  << std::endl;
		extRight = LoganOneDirection (seed2, targetSuffix, querySuffix, scoringScheme, scoreDropOff);

		std::cout << std::endl;
		std::cout << std::endl;

		setBeginPositionH(seed, getEndPositionH(seed) - getEndPositionH(seed1));
		setBeginPositionV(seed, getEndPositionV(seed) - getEndPositionV(seed1));
		setEndPositionH  (seed, getEndPositionH(seed) + getEndPositionH(seed2));
		setEndPositionV  (seed, getEndPositionV(seed) + getEndPositionV(seed2));

		int diffCol  = getEndPositionV(seed) - getBeginPositionV(seed);
		int diffRow  = getEndPositionH(seed) - getBeginPositionH(seed);
		int minLeft  = std::min(getBeginPositionV(seed), getBeginPositionH(seed));
		int minRight = std::min(query.length() - getEndPositionV(seed), target.length() - getEndPositionH(seed));
		int ov = minLeft+minRight+(diffCol+diffRow)/2;

		//std::cout << "ov " << ov << std::endl;

		return extLeft + extRight;
	}
}
#endif
