//==================================================================
// Title:  LOGAN: X-Drop Adaptive Banded Alignment
// Author: G. Guidi, E. Younis
// Date:   30 April 2019
//==================================================================

#include <cstdint>

#ifndef SIMD_UTILSA_H
#define SIMD_UTILSA_H

//======================================================================================
// GLOBAL FUNCTION DECLARATION
//======================================================================================

#ifdef __AVX2__ 	// Compile flag: -mavx2
#define VECTORWIDTH  (32)
#define LOGICALWIDTH (VECTORWIDTH - 1)
#define vector_t    __m256i
#define add_func    _mm256_adds_epi8   // saturated arithmetic
#define sub_func    _mm256_subs_epi8   // saturated arithmetic
#define max_func    _mm256_max_epi8    // max
#define set1_func   _mm256_set1_epi8   // set1 operation
#define blend_func  _mm256_blendv_epi8 // blending operation
#define cmpeq_func  _mm256_cmpeq_epi8  // compare equality operation
#elif __SSE4_2__ 	// Compile flag: -msse4.2
#define VECTORWIDTH  (16)
#define LOGICALWIDTH (VECTORWIDTH - 1)
#define vector_t    __m128i
#define add_func    _mm_adds_epi8 	// saturated arithmetic
#define max_func    _mm_max_epi8  	// max
#define set1_func   _mm_set1_epi8   // set1 operation
#define blend_func  _mm_blendv_epi8 // blending operation
#define cmpeq_func  _mm_cmpeq_epi8  // compare equality operation
#endif

//======================================================================================
// GLOBAL VARIABLE DEFINITION
//======================================================================================

#define NINF  		(std::numeric_limits<int8_t>::min())
#define myRIGHT 	(0)
#define myDOWN  	(1)
#define INTLIMIT	((std::numeric_limits<int8_t>::max()) - 25)
#define MIDDLE 		(LOGICALWIDTH / 2)

//======================================================================================
// SIMD UTILS
//======================================================================================

typedef int8_t element_t;

typedef union {
	vector_t  simd;
	element_t elem[VECTORWIDTH];
} vector_union_t;

void
print_vector_c(vector_t a) {

	vector_union_t tmp;
	tmp.simd = a;

	printf("{");
	for (int i = 0; i < VECTORWIDTH-1; ++i)
		printf("%c,", tmp.elem[i]);
	printf("%c}\n", tmp.elem[VECTORWIDTH-1]);
}

void
print_vector_d(vector_t a) {

	vector_union_t tmp;
	tmp.simd = a;

	printf("{");
	for (int i = 0; i < VECTORWIDTH-1; ++i)
		//printf("%d,", tmp.elem[i]);
		std::cout << (int)tmp.elem[i] << ",";
	printf("%d}\n", tmp.elem[VECTORWIDTH-1]);
}

enum ExtDirectionL
{
	LOGAN_EXTEND_NONE  = 0,
	LOGAN_EXTEND_LEFT  = 1,
	LOGAN_EXTEND_RIGHT = 2,
	LOGAN_EXTEND_BOTH  = 3
};

template <typename T,typename U>
std::pair<T,U> operator+(const std::pair<T,U> & l,const std::pair<T,U> & r) {
	return {l.first+r.first,l.second+r.second};
}

// TODO: REWRITE leftShift/RightShift for int8_t
#ifdef __AVX2__
inline vector_union_t
leftShift (const vector_union_t& a) { // this work for avx2

	vector_union_t b;
	// https://stackoverflow.com/questions/25248766/emulating-shifts-on-32-bytes-with-avxhttps://stackoverflow.com/questions/25248766/emulating-shifts-on-32-bytes-with-avx
	b.simd = _mm256_alignr_epi8(_mm256_permute2x128_si256(a.simd, a.simd, _MM_SHUFFLE(2, 0, 0, 1)), a.simd, 1);
	b.elem[VECTORWIDTH - 1] = NINF;
	return b;
}
inline vector_union_t
rightShift (const vector_union_t& a) { // this work for avx2

	vector_union_t b;
	// https://stackoverflow.com/questions/25248766/emulating-shifts-on-32-bytes-with-avx
	b.simd = _mm256_alignr_epi8(a.simd, _mm256_permute2x128_si256(a.simd, a.simd, _MM_SHUFFLE(0, 0, 2, 0)), 16 - 1);
	b.elem[0] = NINF;
	return b;
}
#elif __SSE4_2__
inline vector_union_t
leftShift (const vector_union_t& a) { // this work for avx2

	vector_union_t b;
	// https://stackoverflow.com/questions/25248766/emulating-shifts-on-32-bytes-with-avxhttps://stackoverflow.com/questions/25248766/emulating-shifts-on-32-bytes-with-avx
	b.simd = _mm256_alignr_epi8(_mm256_permute2x128_si256(a.simd, a.simd, _MM_SHUFFLE(2, 0, 0, 1)), a.simd, 2);
	b.elem[VECTORWIDTH - 1] = NINF;
	return b;
}
inline vector_union_t
rightShift (const vector_union_t& a) { // this work for avx2

	vector_union_t b;
	// https://stackoverflow.com/questions/25248766/emulating-shifts-on-32-bytes-with-avx
	b.simd = _mm256_alignr_epi8(a.simd, _mm256_permute2x128_si256(a.simd, a.simd, _MM_SHUFFLE(0, 0, 2, 0)), 16 - 2);
	b.elem[0] = NINF;
	return b;
}
#endif

static inline void
moveRight (vector_union_t& antiDiag1, vector_union_t& antiDiag2,
	vector_union_t& antiDiag3, int& hoffset, int& voffset, vector_union_t& vqueryh,
		vector_union_t& vqueryv, const int8_t queryh[], const int8_t queryv[])
{
	// (a) shift to the left on query horizontal
	vqueryh = leftShift (vqueryh);
	vqueryh.elem[LOGICALWIDTH-1] = queryh[hoffset++];
	// (b) shift left on updated vector 1 (this places the right-aligned vector 2 as a left-aligned vector 1)
	antiDiag1.simd = antiDiag2.simd;
	antiDiag1 = leftShift (antiDiag1);
	antiDiag2.simd = antiDiag3.simd;
}

static inline void
moveDown (vector_union_t& antiDiag1, vector_union_t& antiDiag2,
	vector_union_t& antiDiag3, int& hoffset, int& voffset, vector_union_t& vqueryh,
		vector_union_t& vqueryv, const int8_t queryh[], const int8_t queryv[])
{
	//(a) shift to the right on query vertical
	vqueryv = rightShift (vqueryv);
	vqueryv.elem[0] = queryv[voffset++];
	//(b) shift to the right on updated vector 2 (this places the left-aligned vector 3 as a right-aligned vector 2)
	antiDiag1.simd = antiDiag2.simd;
	antiDiag2.simd = antiDiag3.simd;
	antiDiag2 = rightShift (antiDiag2);
}

#endif
