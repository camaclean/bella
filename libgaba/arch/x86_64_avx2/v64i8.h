
/**
 * @file v64i8.h
 *
 * @brief struct and _Generic based vector class implementation
 */
#ifndef _V64I8_H_INCLUDED
#define _V64I8_H_INCLUDED

/* include header for intel / amd sse2 instruction sets */
#include <x86intrin.h>

/* 8bit 64cell */
typedef struct v64i8_s {
	__m256i v1;
	__m256i v2;
} v64i8_t;

/* expanders (without argument) */
#define _e_x_v64i8_1(u)
#define _e_x_v64i8_2(u)

/* expanders (without immediate) */
#define _e_v_v64i8_1(a)				(a).v1
#define _e_v_v64i8_2(a)				(a).v2
#define _e_vv_v64i8_1(a, b)			(a).v1, (b).v1
#define _e_vv_v64i8_2(a, b)			(a).v2, (b).v2
#define _e_vvv_v64i8_1(a, b, c)		(a).v1, (b).v1, (c).v1
#define _e_vvv_v64i8_2(a, b, c)		(a).v2, (b).v2, (c).v2

/* expanders with immediate */
#define _e_i_v64i8_1(imm)			(imm)
#define _e_i_v64i8_2(imm)			(imm)
#define _e_vi_v64i8_1(a, imm)		(a).v1, (imm)
#define _e_vi_v64i8_2(a, imm)		(a).v2, (imm)
#define _e_vvi_v64i8_1(a, b, imm)	(a).v1, (b).v1, (imm)
#define _e_vvi_v64i8_2(a, b, imm)	(a).v2, (b).v2, (imm)

/* address calculation macros */
#define _addr_v64i8_1(imm)			( (__m256i *)(imm) )
#define _addr_v64i8_2(imm)			( (__m256i *)(imm) + 1 )
#define _pv_v64i8(ptr)				( _addr_v64i8_1(ptr) )

/* expanders with pointers */
#define _e_p_v64i8_1(ptr)			_addr_v64i8_1(ptr)
#define _e_p_v64i8_2(ptr)			_addr_v64i8_2(ptr)
#define _e_pv_v64i8_1(ptr, a)		_addr_v64i8_1(ptr), (a).v1
#define _e_pv_v64i8_2(ptr, a)		_addr_v64i8_2(ptr), (a).v2

/* expand intrinsic name */
#define _i_v64i8(intrin) 			_mm256_##intrin##_epi8
#define _i_v64i8x(intrin)			_mm256_##intrin##_si256

/* apply */
#define _a_v64i8(intrin, expander, ...) ( \
	(v64i8_t) { \
		_i_v64i8(intrin)(expander##_v64i8_1(__VA_ARGS__)), \
		_i_v64i8(intrin)(expander##_v64i8_2(__VA_ARGS__)) \
	} \
)
#define _a_v64i8x(intrin, expander, ...) ( \
	(v64i8_t) { \
		_i_v64i8x(intrin)(expander##_v64i8_1(__VA_ARGS__)), \
		_i_v64i8x(intrin)(expander##_v64i8_2(__VA_ARGS__)) \
	} \
)
#define _a_v64i8xv(intrin, expander, ...) { \
	_i_v64i8x(intrin)(expander##_v64i8_1(__VA_ARGS__)); \
	_i_v64i8x(intrin)(expander##_v64i8_2(__VA_ARGS__)); \
}

/* load and store */
#define _load_v64i8(...)	_a_v64i8x(load, _e_p, __VA_ARGS__)
#define _loadu_v64i8(...)	_a_v64i8x(loadu, _e_p, __VA_ARGS__)
#define _store_v64i8(...)	_a_v64i8xv(store, _e_pv, __VA_ARGS__)
#define _storeu_v64i8(...)	_a_v64i8xv(storeu, _e_pv, __VA_ARGS__)

/* broadcast */
#define _set_v64i8(...)		_a_v64i8(set1, _e_i, __VA_ARGS__)
#define _zero_v64i8()		_a_v64i8x(setzero, _e_x, _unused)

/* swap (reverse) */
#define _swap_idx_v64i8() ( \
	_mm256_broadcastsi128_si256(_mm_set_epi8( \
		0, 1, 2, 3, 4, 5, 6, 7, \
		8, 9, 10, 11, 12, 13, 14, 15)) \
)
#define _swap_v64i8(a) ( \
	(v64i8_t) { \
		_mm256_permute2x128_si256( \
			_mm256_shuffle_epi8((a).v2, _swap_idx_v64i8()), \
			_mm256_shuffle_epi8((a).v2, _swap_idx_v64i8()), \
			0x01), \
		_mm256_permute2x128_si256( \
			_mm256_shuffle_epi8((a).v1, _swap_idx_v64i8()), \
			_mm256_shuffle_epi8((a).v1, _swap_idx_v64i8()), \
			0x01) \
	} \
)

/* logics */
#define _not_v64i8(...)		_a_v64i8x(not, _e_v, __VA_ARGS__)
#define _and_v64i8(...)		_a_v64i8x(and, _e_vv, __VA_ARGS__)
#define _or_v64i8(...)		_a_v64i8x(or, _e_vv, __VA_ARGS__)
#define _xor_v64i8(...)		_a_v64i8x(xor, _e_vv, __VA_ARGS__)
#define _andn_v64i8(...)	_a_v64i8x(andnot, _e_vv, __VA_ARGS__)

/* arithmetics */
#define _add_v64i8(...)		_a_v64i8(add, _e_vv, __VA_ARGS__)
#define _sub_v64i8(...)		_a_v64i8(sub, _e_vv, __VA_ARGS__)
#define _adds_v64i8(...)	_a_v64i8(adds, _e_vv, __VA_ARGS__)
#define _subs_v64i8(...)	_a_v64i8(subs, _e_vv, __VA_ARGS__)
#define _max_v64i8(...)		_a_v64i8(max, _e_vv, __VA_ARGS__)
#define _min_v64i8(...)		_a_v64i8(min, _e_vv, __VA_ARGS__)

/* shuffle */
#define _shuf_v64i8(...)	_a_v64i8(shuffle, _e_vv, __VA_ARGS__)

/* blend */
#define _sel_v64i8(...)		_a_v64i8(blendv, _e_vvv, __VA_ARGS__)

/* compare */
#define _eq_v64i8(...)		_a_v64i8(cmpeq, _e_vv, __VA_ARGS__)
#define _gt_v64i8(...)		_a_v64i8(cmpgt, _e_vv, __VA_ARGS__)

/* insert and extract */
#define _ins_v64i8(a, val, imm) { \
	if((imm) < sizeof(__m256i)) { \
		(a).v1 = _i_v64i8(insert)((a).v1, (val), (imm)); \
	} else { \
		(a).v2 = _i_v64i8(insert)((a).v2, (val), (imm) - sizeof(__m256i)); \
	} \
}
#define _ext_v64i8(a, imm) ( \
	(int8_t)(((imm) < sizeof(__m256i)) ? ( \
		_i_v64i8(extract)((a).v1, (imm)) \
	) : ( \
		_i_v64i8(extract)((a).v2, (imm) - sizeof(__m256i)) \
	)) \
)

/* shift */
#define _bsl_v64i8(a, imm) ( \
	(v64i8_t) { \
		_mm256_alignr_epi8( \
			(a).v1, \
			_mm256_permute2x128_si256((a).v1, (a).v2, 0x08), \
			15), \
		_mm256_alignr_epi8( \
			(a).v2, \
			_mm256_permute2x128_si256((a).v1, (a).v2, 0x21), \
			15) \
	} \
)
#define _bsr_v64i8(a, imm) ( \
	(v64i8_t) { \
		_mm256_alignr_epi8( \
			_mm256_permute2x128_si256((a).v1, (a).v2, 0x21), \
			(a).v1, \
			1), \
		_mm256_alignr_epi8( \
			_mm256_castsi128_si256(_mm256_extracti128_si256((a).v2, 1)), \
			(a).v2, \
			1) \
	} \
)
#define _shl_v64i8(a, imm) ( \
	(v64i8_t) { \
		_mm256_slli_epi32((a).v1, (imm)), \
		_mm256_slli_epi32((a).v2, (imm)) \
	} \
)
#define _shr_v64i8(a, imm) ( \
	(v64i8_t) { \
		_mm256_srli_epi32((a).v1, (imm)), \
		_mm256_srli_epi32((a).v2, (imm)) \
	} \
)
#define _sal_v64i8(a, imm) ( \
	(v64i8_t) { \
		_mm256_slai_epi32((a).v1, (imm)), \
		_mm256_slai_epi32((a).v2, (imm)) \
	} \
)
#define _sar_v64i8(a, imm) ( \
	(v64i8_t) { \
		_mm256_srai_epi32((a).v1, (imm)), \
		_mm256_srai_epi32((a).v2, (imm)), \
	} \
)

/* mask */
#define _mask_v64i8(a) ( \
	(v64_mask_t) { \
		.m1 = _i_v64i8(movemask)((a).v1), \
		.m2 = _i_v64i8(movemask)((a).v2) \
	} \
)

/* convert */
#define _cvt_v64i16_v64i8(a) ( \
	(v64i8_t) { \
		_mm256_packs_epi16( \
			_mm256_permute2x128_si256((a).v1, (a).v2, 0x20), \
			_mm256_permute2x128_si256((a).v1, (a).v2, 0x31)), \
		_mm256_packs_epi16( \
			_mm256_permute2x128_si256((a).v3, (a).v4, 0x20), \
			_mm256_permute2x128_si256((a).v3, (a).v4, 0x31)) \
	} \
)

/* debug print */
#ifdef DEBUG
#define _print_v64i8(a) { \
	int8_t _tmp[64]; v64i8_t _tmpv = (a); _storeu_v64i8(_tmp, _tmpv); \
	debug("(v64i8_t) %s(%d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, " \
				 "%d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, " \
				 "%d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, " \
				 "%d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d)", \
		#a, \
		_tmp[32 + 31], \
		_tmp[32 + 30], \
		_tmp[32 + 29], \
		_tmp[32 + 28], \
		_tmp[32 + 27], \
		_tmp[32 + 26], \
		_tmp[32 + 25], \
		_tmp[32 + 24], \
		_tmp[32 + 23], \
		_tmp[32 + 22], \
		_tmp[32 + 21], \
		_tmp[32 + 20], \
		_tmp[32 + 19], \
		_tmp[32 + 18], \
		_tmp[32 + 17], \
		_tmp[32 + 16], \
		_tmp[32 + 15], \
		_tmp[32 + 14], \
		_tmp[32 + 13], \
		_tmp[32 + 12], \
		_tmp[32 + 11], \
		_tmp[32 + 10], \
		_tmp[32 + 9], \
		_tmp[32 + 8], \
		_tmp[32 + 7], \
		_tmp[32 + 6], \
		_tmp[32 + 5], \
		_tmp[32 + 4], \
		_tmp[32 + 3], \
		_tmp[32 + 2], \
		_tmp[32 + 1], \
		_tmp[32 + 0], \
		_tmp[31], \
		_tmp[30], \
		_tmp[29], \
		_tmp[28], \
		_tmp[27], \
		_tmp[26], \
		_tmp[25], \
		_tmp[24], \
		_tmp[23], \
		_tmp[22], \
		_tmp[21], \
		_tmp[20], \
		_tmp[19], \
		_tmp[18], \
		_tmp[17], \
		_tmp[16], \
		_tmp[15], \
		_tmp[14], \
		_tmp[13], \
		_tmp[12], \
		_tmp[11], \
		_tmp[10], \
		_tmp[9], \
		_tmp[8], \
		_tmp[7], \
		_tmp[6], \
		_tmp[5], \
		_tmp[4], \
		_tmp[3], \
		_tmp[2], \
		_tmp[1], \
		_tmp[0]); \
}
#else
#define _print_v64i8(x)		;
#endif

#endif /* _V64I8_H_INCLUDED */
/**
 * end of v64i8.h
 */
