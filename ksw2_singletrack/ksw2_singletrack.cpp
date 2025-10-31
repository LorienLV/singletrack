#include "ksw2_singletrack.hpp"
//STL includes
#include <string>
#include <algorithm> 
#include <cmath>
#define unused(val) (void)val;

static inline int32_t get_difference_value_aff2p(int8_t* mat, int32_t i, int32_t j, int n_col, const int *off, int32_t gapo, int32_t gape, int32_t gapo2, int32_t gape2) {
	const int32_t gap_transition_point = std::ceil(
		(static_cast<double>(gapo2 - gapo) / static_cast<double>(gape - gape2)) - 1.0
    );
	const int32_t gap_transition = gap_transition_point * (gape - gape2) - (gapo2 - gapo) - gape2;
	if (i < 0) {
		if 		(j == 0) 					return gapo+gape;
		else if (j  < gap_transition_point)	return gape;
		else if (j == gap_transition_point)	return -gap_transition;
		else								return gape2;
	}
	else if (j < 0) {
		if 		(i == 0) 					return gapo+gape;
		else if (i  < gap_transition_point)	return gape;
		else if (i == gap_transition_point)	return -gap_transition;
		else								return gape2;
	}
	else {
		const int64_t r 	= i + j;
		const int64_t index = static_cast<int64_t>(r * n_col + i - off[r]);
		return static_cast<int32_t>(mat[index]);
	}
}

static inline int32_t get_difference_value_aff(int8_t* mat, int32_t i, int32_t j, int n_col, const int *off, int32_t gapo, int32_t gape) {
	if (i < 0) {
		return (j == 0) ? (gapo+gape) : gape;
	}
	else if (j < 0) {
		return (i == 0) ? (gapo+gape) : gape;
	}
	else {
		const int64_t r 	= i + j;
		const int64_t index = static_cast<int64_t>(r * n_col + i - off[r]);
		return static_cast<int32_t>(mat[index]) + gapo + gape;
	}
}

void ksw2_singletrack_backtrace_affine2p(void* km, int qlen, const uint8_t *query, int tlen, const uint8_t *target, const int8_t *mat,
				   int8_t q, int8_t e, int8_t q2, int8_t e2, const int n_col, const int *off, int8_t *p, int8_t *p2, ksw_extz_t *ez) 
{
	const int32_t match   = static_cast<int32_t>(mat[0]);
    const int32_t mismatch= static_cast<int32_t>(mat[1]);
    const int32_t gapo    = static_cast<int32_t>(-q);
    const int32_t gape    = static_cast<int32_t>(-e);
	const int32_t gapo2   = static_cast<int32_t>(-q2);
	const int32_t gape2   = static_cast<int32_t>(-e2);
	bool in_mmatrix = true; 
    int n_cigar 	= 0;
    int m_cigar 	= ez->m_cigar;
    uint32_t* cigar = ez->cigar;
	int32_t i = tlen-1, j = qlen-1, l = 0;
	int32_t diff_del = 0, diff_ins = 0;
	while (i >= 0 && j >= 0) {
		if (in_mmatrix) {
			int32_t h    = get_difference_value_aff2p(p,  i, j,   n_col, off, gapo, gape, gapo2, gape2);
			int32_t v    = get_difference_value_aff2p(p2, i, j-1, n_col, off, gapo, gape, gapo2, gape2);
			int32_t cost = (target[i] == query[j]) ? match : mismatch;
					cost = (target[i] == 4 || query[j] == 4) ? gape2 : cost;
			int32_t sum  = h + v;
			if ((i >= 0) && (j >= 0) && (sum == cost)) {
				cigar = ksw_push_cigar(km, &n_cigar, &m_cigar, cigar, KSW_CIGAR_MATCH, 1);
				j--;
				i--;
			}
			else {
				in_mmatrix	= false; 
				l			= 0;
				diff_ins	= get_difference_value_aff2p(p,  i, j, n_col, off, gapo, gape, gapo2, gape2);
				diff_del	= get_difference_value_aff2p(p2, i, j, n_col, off, gapo, gape, gapo2, gape2);
			}
		}
		else {
			l++;
			int32_t acc_af = gapo  + l * gape;
			int32_t acc_af2= gapo2 + l * gape2;
			int32_t j_ins  = j - l;
			int32_t i_del  = i - l;
			if ((i_del >= 0) && ((diff_del == acc_af) || (diff_del == acc_af2))) {
				cigar = ksw_push_cigar(km, &n_cigar, &m_cigar, cigar, KSW_CIGAR_DEL, l);
				i = i_del;
				in_mmatrix = true;
			}
			else if ((j_ins >= 0) && ((diff_ins == acc_af) || (diff_ins == acc_af2))) {
				cigar = ksw_push_cigar(km, &n_cigar, &m_cigar, cigar, KSW_CIGAR_INS, l);
				j = j_ins; 
				in_mmatrix = true; 
			}
			else {
				if (j_ins >= 0) {
					int32_t H = get_difference_value_aff2p(p, i, j_ins, n_col, off, gapo, gape, gapo2, gape2);
					diff_ins += H;
				}
				if (i_del >= 0) {
					int32_t V = get_difference_value_aff2p(p2, i_del, j, n_col, off, gapo, gape, gapo2, gape2);
					diff_del += V;
				}
			}
		}
	}
    if (i >= 0) cigar = ksw_push_cigar(km, &n_cigar, &m_cigar, cigar, KSW_CIGAR_DEL, i+1);
    if (j >= 0) cigar = ksw_push_cigar(km, &n_cigar, &m_cigar, cigar, KSW_CIGAR_INS, j+1);
	for(i = 0; i < n_cigar>>1; ++i) {
		uint32_t tmp = cigar[i];
		cigar[i] = cigar[n_cigar-1-i];
		cigar[n_cigar-1-i] = tmp;
	}
	ez->cigar = cigar;
    ez->m_cigar = m_cigar;
    ez->n_cigar = n_cigar;
}

void ksw2_singletrack_backtrace_affine(void* km, int qlen, const uint8_t *query, int tlen, const uint8_t *target, const int8_t *mat,
				   int8_t q, int8_t e, const int n_col, const int *off, int8_t *p, int8_t *p2, ksw_extz_t *ez) 
{
	const int32_t match   = static_cast<int32_t>(mat[0]);
    const int32_t mismatch= static_cast<int32_t>(mat[1]);
    const int32_t gapo    = static_cast<int32_t>(-q);
    const int32_t gape    = static_cast<int32_t>(-e);
	bool in_mmatrix = true; 
    int n_cigar 	= 0;
    int m_cigar 	= ez->m_cigar;
    uint32_t* cigar = ez->cigar;
	int32_t i = tlen-1, j = qlen-1, l = 0;
	int32_t diff_del = 0, diff_ins = 0;
	while (i >= 0 && j >= 0) {
		if (in_mmatrix) {
			int32_t h    = get_difference_value_aff(p,  i, j,   n_col, off, gapo, gape);
			int32_t v    = get_difference_value_aff(p2, i, j-1, n_col, off, gapo, gape);
			int32_t cost = (target[i] == query[j]) ? match : mismatch;
					cost = (target[i] == 4 || query[j] == 4) ? gape : cost;
			int32_t sum  = h + v;
			if ((i >= 0) && (j >= 0) && (sum == cost)) {
				cigar = ksw_push_cigar(km, &n_cigar, &m_cigar, cigar, KSW_CIGAR_MATCH, 1);
				j--;
				i--;
			}
			else {
				in_mmatrix	= false; 
				l			= 0;
				diff_ins	= get_difference_value_aff(p,  i, j, n_col, off, gapo, gape);
				diff_del	= get_difference_value_aff(p2, i, j, n_col, off, gapo, gape);
			}
		}
		else {
			l++;
			int32_t acc_af = gapo + l * gape;
			int32_t j_ins  = j - l;
			int32_t i_del  = i - l;
			if ((i_del >= 0) && (diff_del == acc_af)) {
                cigar = ksw_push_cigar(km, &n_cigar, &m_cigar, cigar, KSW_CIGAR_DEL, l);
				i = i_del;
				in_mmatrix = true;
			}
			else if ((j_ins >= 0) && (diff_ins == acc_af)) {
                cigar = ksw_push_cigar(km, &n_cigar, &m_cigar, cigar, KSW_CIGAR_INS, l);
				j = j_ins;
				in_mmatrix = true;
			}
			else {
				if (j_ins >= 0) {
					int32_t H = get_difference_value_aff(p, i, j_ins, n_col, off, gapo, gape);
					diff_ins += H;
				}
				if (i_del >= 0) {
					int32_t V = get_difference_value_aff(p2, i_del, j, n_col, off, gapo, gape);
					diff_del += V;
				}
			}
		}
	}
    if (i >= 0) cigar = ksw_push_cigar(km, &n_cigar, &m_cigar, cigar, KSW_CIGAR_DEL, i+1);
    if (j >= 0) cigar = ksw_push_cigar(km, &n_cigar, &m_cigar, cigar, KSW_CIGAR_INS, j+1);
	for(i = 0; i < n_cigar>>1; ++i) {
		uint32_t tmp = cigar[i];
		cigar[i] = cigar[n_cigar-1-i];
		cigar[n_cigar-1-i] = tmp;
	}
	ez->cigar = cigar;
    ez->m_cigar = m_cigar;
    ez->n_cigar = n_cigar;
}

#ifdef __SSE2__
#ifdef USE_SIMDE
#include <simde/x86/sse2.h>
#else
#include <emmintrin.h>
#endif

#ifdef KSW_SSE2_ONLY
#undef __SSE4_1__
#endif

#ifdef __SSE4_1__
#ifdef USE_SIMDE
#include <simde/x86/sse4.1.h>
#else
#include <smmintrin.h>
#endif
#endif

void ksw_extd2_singletrack_sse(void *km, int qlen, const uint8_t *query, int tlen, const uint8_t *target, int8_t m, const int8_t *mat,
				   int8_t q, int8_t e, int8_t q2, int8_t e2, int w, int zdrop, int end_bonus, int flag, ksw_extz_t *ez)
{
#define __dp_code_block1 \
	z = _mm_load_si128(&s[t]); \
	xt1 = _mm_load_si128(&x[t]);                     /* xt1 <- x[r-1][t..t+15] */ \
	tmp = _mm_srli_si128(xt1, 15);                   /* tmp <- x[r-1][t+15] */ \
	xt1 = _mm_or_si128(_mm_slli_si128(xt1, 1), x1_); /* xt1 <- x[r-1][t-1..t+14] */ \
	x1_ = tmp; \
	vt1 = _mm_load_si128(&v[t]);                     /* vt1 <- v[r-1][t..t+15] */ \
	tmp = _mm_srli_si128(vt1, 15);                   /* tmp <- v[r-1][t+15] */ \
	vt1 = _mm_or_si128(_mm_slli_si128(vt1, 1), v1_); /* vt1 <- v[r-1][t-1..t+14] */ \
	v1_ = tmp; \
	a = _mm_add_epi8(xt1, vt1);                      /* a <- x[r-1][t-1..t+14] + v[r-1][t-1..t+14] */ \
	ut = _mm_load_si128(&u[t]);                      /* ut <- u[t..t+15] */ \
	b = _mm_add_epi8(_mm_load_si128(&y[t]), ut);     /* b <- y[r-1][t..t+15] + u[r-1][t..t+15] */ \
	x2t1= _mm_load_si128(&x2[t]); \
	tmp = _mm_srli_si128(x2t1, 15); \
	x2t1= _mm_or_si128(_mm_slli_si128(x2t1, 1), x21_); \
	x21_= tmp; \
	a2= _mm_add_epi8(x2t1, vt1); \
	b2= _mm_add_epi8(_mm_load_si128(&y2[t]), ut);

#define __dp_code_block2 \
	tmp2 = _mm_sub_epi8(z, vt1);\
	tmp3 = _mm_sub_epi8(z, ut);\
	_mm_store_si128(&u[t], tmp2);    /* u[r][t..t+15] <- z - v[r-1][t-1..t+14] */ \
	_mm_store_si128(&v[t], tmp3);     /* v[r][t..t+15] <- z - u[r-1][t..t+15] */ \
	tmp = _mm_sub_epi8(z, q_); \
	a = _mm_sub_epi8(a, tmp); \
	b = _mm_sub_epi8(b, tmp); \
	tmp = _mm_sub_epi8(z, q2_); \
	a2= _mm_sub_epi8(a2, tmp); \
	b2= _mm_sub_epi8(b2, tmp);

	unused(end_bonus);
	int r, t, qe = q + e, n_col_, *off = 0, *off_end = 0, tlen_, qlen_, last_st, last_en, wl, wr, max_sc, min_sc, long_thres, long_diff;
	int with_cigar = !(flag&KSW_EZ_SCORE_ONLY), approx_max = !!(flag&KSW_EZ_APPROX_MAX);
	int32_t *H = 0, H0 = 0, last_H0_t = 0;
	uint8_t *qr, *sf, *mem, *mem2 = 0;
	__m128i q_, q2_, qe_, qe2_, zero_, sc_mch_, sc_mis_, m1_, sc_N_;
	__m128i *u, *v, *x, *y, *x2, *y2, *s, *p = 0, *p2 = 0;

	ksw_reset_extz(ez);
	if (m <= 1 || qlen <= 0 || tlen <= 0) return;

	if (q2 + e2 < q + e) t = q, q = q2, q2 = t, t = e, e = e2, e2 = t; // make sure q+e no larger than q2+e2

	zero_   = _mm_set1_epi8(0);
	q_      = _mm_set1_epi8(q);
	q2_     = _mm_set1_epi8(q2);
	qe_     = _mm_set1_epi8(q + e);
	qe2_    = _mm_set1_epi8(q2 + e2);
	sc_mch_ = _mm_set1_epi8(mat[0]);
	sc_mis_ = _mm_set1_epi8(mat[1]);
	sc_N_   = mat[m*m-1] == 0? _mm_set1_epi8(-e2) : _mm_set1_epi8(mat[m*m-1]);
	m1_     = _mm_set1_epi8(m - 1); // wildcard

	if (w < 0) w = tlen > qlen? tlen : qlen;
	wl = wr = w;
	tlen_ = (tlen + 15) / 16;
	n_col_ = qlen < tlen? qlen : tlen;
	n_col_ = ((n_col_ < w + 1? n_col_ : w + 1) + 15) / 16 + 1;
	qlen_ = (qlen + 15) / 16;
	for (t = 1, max_sc = mat[0], min_sc = mat[1]; t < m * m; ++t) {
		max_sc = max_sc > mat[t]? max_sc : mat[t];
		min_sc = min_sc < mat[t]? min_sc : mat[t];
	}
	if (-min_sc > 2 * (q + e)) return; // otherwise, we won't see any mismatches

	long_thres = e != e2? (q2 - q) / (e - e2) - 1 : 0;
	if (q2 + e2 + long_thres * e2 > q + e + long_thres * e)
		++long_thres;
	long_diff = long_thres * (e - e2) - (q2 - q) - e2;

	mem = (uint8_t*)kcalloc(km, tlen_ * 8 + qlen_ + 1, 16);
	u = (__m128i*)(((size_t)mem + 15) >> 4 << 4); // 16-byte aligned
	v = u + tlen_, x = v + tlen_, y = x + tlen_, x2 = y + tlen_, y2 = x2 + tlen_;
	s = y2 + tlen_, sf = (uint8_t*)(s + tlen_), qr = sf + tlen_ * 16;
	memset(u,  -q  - e,  tlen_ * 16);
	memset(v,  -q  - e,  tlen_ * 16);
	memset(x,  -q  - e,  tlen_ * 16);
	memset(y,  -q  - e,  tlen_ * 16);
	memset(x2, -q2 - e2, tlen_ * 16);
	memset(y2, -q2 - e2, tlen_ * 16);
	if (!approx_max) {
		H = (int32_t*)kmalloc(km, tlen_ * 16 * 4);
		for (t = 0; t < tlen_ * 16; ++t) H[t] = KSW_NEG_INF;
	}
	if (with_cigar) {
		mem2 = (uint8_t*)kmalloc(km, ((size_t)(qlen + tlen - 1) * n_col_ + 1) * 16 * 2);
		p = (__m128i*)(((size_t)mem2 + 15) >> 4 << 4);
		p2= p + ((qlen + tlen - 1) * n_col_ + 1);
		off = (int*)kmalloc(km, (qlen + tlen - 1) * sizeof(int) * 2);
		off_end = off + qlen + tlen - 1;
	}

	for (t = 0; t < qlen; ++t) qr[t] = query[qlen - 1 - t];
	memcpy(sf, target, tlen);

	for (r = 0, last_st = last_en = -1; r < qlen + tlen - 1; ++r) {
		int st = 0, en = tlen - 1, st0, en0, st_, en_;
		int8_t x1, x21, v1;
		uint8_t *qrr = qr + (qlen - 1 - r);
		int8_t *u8 = (int8_t*)u, *v8 = (int8_t*)v, *x8 = (int8_t*)x, *x28 = (int8_t*)x2;
		__m128i x1_, x21_, v1_;
		// find the boundaries
		if (st < r - qlen + 1) st = r - qlen + 1;
		if (en > r) en = r;
		if (st < (r-wr+1)>>1) st = (r-wr+1)>>1; // take the ceil
		if (en > (r+wl)>>1) en = (r+wl)>>1; // take the floor
		if (st > en) {
			ez->zdropped = 1;
			break;
		}
		st0 = st, en0 = en;
		st = st / 16 * 16, en = (en + 16) / 16 * 16 - 1;
		// set boundary conditions
		if (st > 0) {
			if (st - 1 >= last_st && st - 1 <= last_en) {
				x1 = x8[st - 1], x21 = x28[st - 1], v1 = v8[st - 1]; // (r-1,s-1) calculated in the last round
			} else {
				x1 = -q - e, x21 = -q2 - e2;
				v1 = -q - e;
			}
		} else {
			x1 = -q - e, x21 = -q2 - e2;
			v1 = r == 0? -q - e : r < long_thres? -e : r == long_thres? long_diff : -e2;
		}
		if (en >= r) {
			((int8_t*)y)[r] = -q - e, ((int8_t*)y2)[r] = -q2 - e2;
			u8[r] = r == 0? -q - e : r < long_thres? -e : r == long_thres? long_diff : -e2;
		}
		// loop fission: set scores first
		if (!(flag & KSW_EZ_GENERIC_SC)) {
			for (t = st0; t <= en0; t += 16) {
				__m128i sq, st, tmp, mask;
				sq = _mm_loadu_si128((__m128i*)&sf[t]);
				st = _mm_loadu_si128((__m128i*)&qrr[t]);
				mask = _mm_or_si128(_mm_cmpeq_epi8(sq, m1_), _mm_cmpeq_epi8(st, m1_));
				tmp = _mm_cmpeq_epi8(sq, st);
#ifdef __SSE4_1__
				tmp = _mm_blendv_epi8(sc_mis_, sc_mch_, tmp);
				tmp = _mm_blendv_epi8(tmp,     sc_N_,   mask);
#else
				tmp = _mm_or_si128(_mm_andnot_si128(tmp,  sc_mis_), _mm_and_si128(tmp,  sc_mch_));
				tmp = _mm_or_si128(_mm_andnot_si128(mask, tmp),     _mm_and_si128(mask, sc_N_));
#endif
				_mm_storeu_si128((__m128i*)((int8_t*)s + t), tmp);
			}
		} else {
			for (t = st0; t <= en0; ++t)
				((uint8_t*)s)[t] = mat[sf[t] * m + qrr[t]];
		}
		// core loop
		x1_  = _mm_cvtsi32_si128((uint8_t)x1);
		x21_ = _mm_cvtsi32_si128((uint8_t)x21);
		v1_  = _mm_cvtsi32_si128((uint8_t)v1);
		st_ = st / 16, en_ = en / 16;
		assert(en_ - st_ + 1 <= n_col_);
		if (!with_cigar) { // score only
			for (t = st_; t <= en_; ++t) {
				__m128i z, a, b, a2, b2, xt1, x2t1, vt1, ut, tmp, tmp2, tmp3;
				__dp_code_block1;
#ifdef __SSE4_1__
				z = _mm_max_epi8(z, a);
				z = _mm_max_epi8(z, b);
				z = _mm_max_epi8(z, a2);
				z = _mm_max_epi8(z, b2);
				z = _mm_min_epi8(z, sc_mch_);
				__dp_code_block2; // save u[] and v[]; update a, b, a2 and b2
				_mm_store_si128(&x[t],  _mm_sub_epi8(_mm_max_epi8(a,  zero_), qe_));
				_mm_store_si128(&y[t],  _mm_sub_epi8(_mm_max_epi8(b,  zero_), qe_));
				_mm_store_si128(&x2[t], _mm_sub_epi8(_mm_max_epi8(a2, zero_), qe2_));
				_mm_store_si128(&y2[t], _mm_sub_epi8(_mm_max_epi8(b2, zero_), qe2_));
#else
				tmp = _mm_cmpgt_epi8(a,  z);
				z = _mm_or_si128(_mm_andnot_si128(tmp, z), _mm_and_si128(tmp, a));
				tmp = _mm_cmpgt_epi8(b,  z);
				z = _mm_or_si128(_mm_andnot_si128(tmp, z), _mm_and_si128(tmp, b));
				tmp = _mm_cmpgt_epi8(a2, z);
				z = _mm_or_si128(_mm_andnot_si128(tmp, z), _mm_and_si128(tmp, a2));
				tmp = _mm_cmpgt_epi8(b2, z);
				z = _mm_or_si128(_mm_andnot_si128(tmp, z), _mm_and_si128(tmp, b2));
				tmp = _mm_cmplt_epi8(sc_mch_, z);
				z = _mm_or_si128(_mm_and_si128(tmp, sc_mch_), _mm_andnot_si128(tmp, z));
				__dp_code_block2;
				tmp = _mm_cmpgt_epi8(a, zero_);
				_mm_store_si128(&x[t],  _mm_sub_epi8(_mm_and_si128(tmp, a),  qe_));
				tmp = _mm_cmpgt_epi8(b, zero_);
				_mm_store_si128(&y[t],  _mm_sub_epi8(_mm_and_si128(tmp, b),  qe_));
				tmp = _mm_cmpgt_epi8(a2, zero_);
				_mm_store_si128(&x2[t], _mm_sub_epi8(_mm_and_si128(tmp, a2), qe2_));
				tmp = _mm_cmpgt_epi8(b2, zero_);
				_mm_store_si128(&y2[t], _mm_sub_epi8(_mm_and_si128(tmp, b2), qe2_));
#endif
			} 
		} else if (!(flag&KSW_EZ_RIGHT)) { // gap left-alignment
			__m128i *pr = p + (size_t)r * n_col_ - st_;
			__m128i *pr2= p2+ (size_t)r * n_col_ - st_;
			off[r] = st, off_end[r] = en;
			for (t = st_; t <= en_; ++t) {
				__m128i z, a, b, a2, b2, xt1, x2t1, vt1, ut, tmp, tmp2, tmp3;
				__dp_code_block1;
#ifdef __SSE4_1__
				z = _mm_max_epi8(z, a);
				z = _mm_max_epi8(z, b);
				z = _mm_max_epi8(z, a2);
				z = _mm_max_epi8(z, b2);
				z = _mm_min_epi8(z, sc_mch_);
#else // we need to emulate SSE4.1 intrinsics _mm_max_epi8() and _mm_blendv_epi8()
				tmp = _mm_cmpgt_epi8(a,  z);
				z = _mm_or_si128(_mm_andnot_si128(tmp, z), _mm_and_si128(tmp, a));
				tmp = _mm_cmpgt_epi8(b,  z);
				z = _mm_or_si128(_mm_andnot_si128(tmp, z), _mm_and_si128(tmp, b));
				tmp = _mm_cmpgt_epi8(a2, z);
				z = _mm_or_si128(_mm_andnot_si128(tmp, z), _mm_and_si128(tmp, a2));
				tmp = _mm_cmpgt_epi8(b2, z);
				z = _mm_or_si128(_mm_andnot_si128(tmp, z), _mm_and_si128(tmp, b2));
				tmp = _mm_cmplt_epi8(sc_mch_, z);
				z = _mm_or_si128(_mm_and_si128(tmp, sc_mch_), _mm_andnot_si128(tmp, z));
#endif
				__dp_code_block2;
				tmp = _mm_cmpgt_epi8(a, zero_);
				_mm_store_si128(&x[t],  _mm_sub_epi8(_mm_and_si128(tmp, a),  qe_));
				tmp = _mm_cmpgt_epi8(b, zero_);
				_mm_store_si128(&y[t],  _mm_sub_epi8(_mm_and_si128(tmp, b),  qe_));
				tmp = _mm_cmpgt_epi8(a2, zero_);
				_mm_store_si128(&x2[t], _mm_sub_epi8(_mm_and_si128(tmp, a2), qe2_));
				tmp = _mm_cmpgt_epi8(b2, zero_);
				_mm_store_si128(&y2[t], _mm_sub_epi8(_mm_and_si128(tmp, b2), qe2_));
				_mm_store_si128(&pr[t], tmp3); 
				_mm_store_si128(&pr2[t], tmp2);
			}
		} else { // gap right-alignment
			__m128i *pr = p + (size_t)r * n_col_ - st_;
			__m128i *pr2= p2+ (size_t)r * n_col_ - st_;
			off[r] = st, off_end[r] = en;
			for (t = st_; t <= en_; ++t) {
				__m128i z, a, b, a2, b2, xt1, x2t1, vt1, ut, tmp, tmp2, tmp3;
				__dp_code_block1;
#ifdef __SSE4_1__
				z = _mm_max_epi8(z, a);
				z = _mm_max_epi8(z, b);
				z = _mm_max_epi8(z, a2);
				z = _mm_max_epi8(z, b2);
				z = _mm_min_epi8(z, sc_mch_);
#else // we need to emulate SSE4.1 intrinsics _mm_max_epi8() and _mm_blendv_epi8()
				tmp = _mm_cmpgt_epi8(z, a);
				z = _mm_or_si128(_mm_and_si128(tmp, z), _mm_andnot_si128(tmp, a));
				tmp = _mm_cmpgt_epi8(z, b);
				z = _mm_or_si128(_mm_and_si128(tmp, z), _mm_andnot_si128(tmp, b));
				tmp = _mm_cmpgt_epi8(z, a2);
				z = _mm_or_si128(_mm_and_si128(tmp, z), _mm_andnot_si128(tmp, a2));
				tmp = _mm_cmpgt_epi8(z, b2);
				z = _mm_or_si128(_mm_and_si128(tmp, z), _mm_andnot_si128(tmp, b2));
				tmp = _mm_cmplt_epi8(sc_mch_, z);
				z = _mm_or_si128(_mm_and_si128(tmp, sc_mch_), _mm_andnot_si128(tmp, z));
#endif
				__dp_code_block2;
				tmp = _mm_cmpgt_epi8(zero_, a);
				_mm_store_si128(&x[t],  _mm_sub_epi8(_mm_andnot_si128(tmp, a),  qe_));
				tmp = _mm_cmpgt_epi8(zero_, b);
				_mm_store_si128(&y[t],  _mm_sub_epi8(_mm_andnot_si128(tmp, b),  qe_));
				tmp = _mm_cmpgt_epi8(zero_, a2);
				_mm_store_si128(&x2[t], _mm_sub_epi8(_mm_andnot_si128(tmp, a2), qe2_));
				tmp = _mm_cmpgt_epi8(zero_, b2);
				_mm_store_si128(&y2[t], _mm_sub_epi8(_mm_andnot_si128(tmp, b2), qe2_));
				_mm_store_si128(&pr[t], tmp3); 
				_mm_store_si128(&pr2[t], tmp2);
			}
		}
		if (!approx_max) { // find the exact max with a 32-bit score array
			int32_t max_H, max_t;
			// compute H[], max_H and max_t
			if (r > 0) {
				int32_t HH[4], tt[4], en1 = st0 + (en0 - st0) / 4 * 4, i;
				__m128i max_H_, max_t_;
				max_H = H[en0] = en0 > 0? H[en0-1] + u8[en0] : H[en0] + v8[en0]; // special casing the last element
				max_t = en0;
				max_H_ = _mm_set1_epi32(max_H);
				max_t_ = _mm_set1_epi32(max_t);
				for (t = st0; t < en1; t += 4) { // this implements: H[t]+=v8[t]-qe; if(H[t]>max_H) max_H=H[t],max_t=t;
					__m128i H1, tmp, t_;
					H1 = _mm_loadu_si128((__m128i*)&H[t]);
					t_ = _mm_setr_epi32(v8[t], v8[t+1], v8[t+2], v8[t+3]);
					H1 = _mm_add_epi32(H1, t_);
					_mm_storeu_si128((__m128i*)&H[t], H1);
					t_ = _mm_set1_epi32(t);
					tmp = _mm_cmpgt_epi32(H1, max_H_);
#ifdef __SSE4_1__
					max_H_ = _mm_blendv_epi8(max_H_, H1, tmp);
					max_t_ = _mm_blendv_epi8(max_t_, t_, tmp);
#else
					max_H_ = _mm_or_si128(_mm_and_si128(tmp, H1), _mm_andnot_si128(tmp, max_H_));
					max_t_ = _mm_or_si128(_mm_and_si128(tmp, t_), _mm_andnot_si128(tmp, max_t_));
#endif
				}
				_mm_storeu_si128((__m128i*)HH, max_H_);
				_mm_storeu_si128((__m128i*)tt, max_t_);
				for (i = 0; i < 4; ++i)
					if (max_H < HH[i]) max_H = HH[i], max_t = tt[i] + i;
				for (; t < en0; ++t) { // for the rest of values that haven't been computed with SSE
					H[t] += (int32_t)v8[t];
					if (H[t] > max_H)
						max_H = H[t], max_t = t;
				}
			} else H[0] = v8[0] - qe, max_H = H[0], max_t = 0; // special casing r==0
			// update ez
			if (en0 == tlen - 1 && H[en0] > ez->mte)
				ez->mte = H[en0], ez->mte_q = r - en;
			if (r - st0 == qlen - 1 && H[st0] > ez->mqe)
				ez->mqe = H[st0], ez->mqe_t = st0;
			if (ksw_apply_zdrop(ez, 1, max_H, r, max_t, zdrop, e2)) break;
			if (r == qlen + tlen - 2 && en0 == tlen - 1)
				ez->score = H[tlen - 1];
		} else { // find approximate max; Z-drop might be inaccurate, too.
			if (r > 0) {
				if (last_H0_t >= st0 && last_H0_t <= en0 && last_H0_t + 1 >= st0 && last_H0_t + 1 <= en0) {
					int32_t d0 = v8[last_H0_t];
					int32_t d1 = u8[last_H0_t + 1];
					if (d0 > d1) H0 += d0;
					else H0 += d1, ++last_H0_t;
				} else if (last_H0_t >= st0 && last_H0_t <= en0) {
					H0 += v8[last_H0_t];
				} else {
					++last_H0_t, H0 += u8[last_H0_t];
				}
			} else H0 = v8[0] - qe, last_H0_t = 0;
			if ((flag & KSW_EZ_APPROX_DROP) && ksw_apply_zdrop(ez, 1, H0, r, last_H0_t, zdrop, e2)) break;
			if (r == qlen + tlen - 2 && en0 == tlen - 1)
				ez->score = H0;
		}
		last_st = st, last_en = en;
		//for (t = st0; t <= en0; ++t) printf("(%d,%d)\t(%d,%d,%d,%d)\t%d\n", r, t, ((int8_t*)u)[t], ((int8_t*)v)[t], ((int8_t*)x)[t], ((int8_t*)y)[t], H[t]); // for debugging
	}
	kfree(km, mem);
	if (!approx_max) kfree(km, H);
	if (with_cigar) { // backtrack
		ksw2_singletrack_backtrace_affine2p(km, qlen, query, tlen, target, mat, q, e, q2, e2, n_col_*16, off, reinterpret_cast<int8_t*>(p), reinterpret_cast<int8_t*>(p2), ez);
		kfree(km, mem2); kfree(km, off);
	}
    #undef __dp_code_block1
    #undef __dp_code_block2
}

void ksw_extz2_singletrack_sse(void *km, int qlen, const uint8_t *query, int tlen, const uint8_t *target, int8_t m, const int8_t *mat, 
					int8_t q, int8_t e, int w, int zdrop, int end_bonus, int flag, ksw_extz_t *ez)
{
#define __dp_code_block1 \
	z = _mm_add_epi8(_mm_load_si128(&s[t]), qe2_); \
	xt1 = _mm_load_si128(&x[t]);                     /* xt1 <- x[r-1][t..t+15] */ \
	tmp = _mm_srli_si128(xt1, 15);                   /* tmp <- x[r-1][t+15] */ \
	xt1 = _mm_or_si128(_mm_slli_si128(xt1, 1), x1_); /* xt1 <- x[r-1][t-1..t+14] */ \
	x1_ = tmp; \
	vt1 = _mm_load_si128(&v[t]);                     /* vt1 <- v[r-1][t..t+15] */ \
	tmp = _mm_srli_si128(vt1, 15);                   /* tmp <- v[r-1][t+15] */ \
	vt1 = _mm_or_si128(_mm_slli_si128(vt1, 1), v1_); /* vt1 <- v[r-1][t-1..t+14] */ \
	v1_ = tmp; \
	a = _mm_add_epi8(xt1, vt1);                      /* a <- x[r-1][t-1..t+14] + v[r-1][t-1..t+14] */ \
	ut = _mm_load_si128(&u[t]);                      /* ut <- u[t..t+15] */ \
	b = _mm_add_epi8(_mm_load_si128(&y[t]), ut);     /* b <- y[r-1][t..t+15] + u[r-1][t..t+15] */

#define __dp_code_block2 \
	z = _mm_max_epu8(z, b);                          /* z = max(z, b); this works because both are non-negative */ \
	z = _mm_min_epu8(z, max_sc_); \
	tmp2 = _mm_sub_epi8(z, vt1);\
	tmp3 = _mm_sub_epi8(z, ut);\
	_mm_store_si128(&u[t], tmp2);    /* u[r][t..t+15] <- z - v[r-1][t-1..t+14] */ \
	_mm_store_si128(&v[t], tmp3);     /* v[r][t..t+15] <- z - u[r-1][t..t+15] */ \
	z = _mm_sub_epi8(z, q_); \
	a = _mm_sub_epi8(a, z); \
	b = _mm_sub_epi8(b, z);

	unused(end_bonus);
	int r, t, qe = q + e, n_col_, *off = 0, *off_end = 0, tlen_, qlen_, last_st, last_en, wl, wr, max_sc, min_sc;
	int with_cigar = !(flag&KSW_EZ_SCORE_ONLY), approx_max = !!(flag&KSW_EZ_APPROX_MAX);
	int32_t *H = 0, H0 = 0, last_H0_t = 0;
	uint8_t *qr, *sf, *mem, *mem2 = 0;
	__m128i q_, qe2_, zero_, sc_mch_, sc_mis_, sc_N_, m1_, max_sc_;
	__m128i *u, *v, *x, *y, *s, *p = 0, *p2 = 0;

	ksw_reset_extz(ez);
	if (m <= 0 || qlen <= 0 || tlen <= 0) return;

	zero_   = _mm_set1_epi8(0);
	q_      = _mm_set1_epi8(q);
	qe2_    = _mm_set1_epi8((q + e) * 2);
	sc_mch_ = _mm_set1_epi8(mat[0]);
	sc_mis_ = _mm_set1_epi8(mat[1]);
	sc_N_   = mat[m*m-1] == 0? _mm_set1_epi8(-e) : _mm_set1_epi8(mat[m*m-1]);
	m1_     = _mm_set1_epi8(m - 1); // wildcard
	max_sc_ = _mm_set1_epi8(mat[0] + (q + e) * 2);

	if (w < 0) w = tlen > qlen? tlen : qlen;
	wl = wr = w;
	tlen_ = (tlen + 15) / 16;
	n_col_ = qlen < tlen? qlen : tlen;
	n_col_ = ((n_col_ < w + 1? n_col_ : w + 1) + 15) / 16 + 1;
	qlen_ = (qlen + 15) / 16;
	for (t = 1, max_sc = mat[0], min_sc = mat[1]; t < m * m; ++t) {
		max_sc = max_sc > mat[t]? max_sc : mat[t];
		min_sc = min_sc < mat[t]? min_sc : mat[t];
	}
	if (-min_sc > 2 * (q + e)) return; // otherwise, we won't see any mismatches

	mem = (uint8_t*)kcalloc(km, tlen_ * 6 + qlen_ + 1, 16);
	u = (__m128i*)(((size_t)mem + 15) >> 4 << 4); // 16-byte aligned
	v = u + tlen_, x = v + tlen_, y = x + tlen_, s = y + tlen_, sf = (uint8_t*)(s + tlen_), qr = sf + tlen_ * 16;
	if (!approx_max) {
		H = (int32_t*)kmalloc(km, tlen_ * 16 * 4);
		for (t = 0; t < tlen_ * 16; ++t) H[t] = KSW_NEG_INF;
	}
	if (with_cigar) {
		mem2 = (uint8_t*)kmalloc(km, ((size_t)(qlen + tlen - 1) * n_col_ + 1) * 16 * 2);
		p = (__m128i*)(((size_t)mem2 + 15) >> 4 << 4);
		p2= p + ((qlen + tlen - 1) * n_col_ + 1);
		off = (int*)kmalloc(km, (qlen + tlen - 1) * sizeof(int) * 2);
		off_end = off + qlen + tlen - 1;
	}

	for (t = 0; t < qlen; ++t) qr[t] = query[qlen - 1 - t];
	memcpy(sf, target, tlen);

	for (r = 0, last_st = last_en = -1; r < qlen + tlen - 1; ++r) {
		int st = 0, en = tlen - 1, st0, en0, st_, en_;
		int8_t x1, v1;
		uint8_t *qrr = qr + (qlen - 1 - r), *u8 = (uint8_t*)u, *v8 = (uint8_t*)v;
		__m128i x1_, v1_;
		// find the boundaries
		if (st < r - qlen + 1) st = r - qlen + 1;
		if (en > r) en = r;
		if (st < (r-wr+1)>>1) st = (r-wr+1)>>1; // take the ceil
		if (en > (r+wl)>>1) en = (r+wl)>>1; // take the floor
		if (st > en) {
			ez->zdropped = 1;
			break;
		}
		st0 = st, en0 = en;
		st = st / 16 * 16, en = (en + 16) / 16 * 16 - 1;
		// set boundary conditions
		if (st > 0) {
			if (st - 1 >= last_st && st - 1 <= last_en)
				x1 = ((uint8_t*)x)[st - 1], v1 = v8[st - 1]; // (r-1,s-1) calculated in the last round
			else x1 = v1 = 0; // not calculated; set to zeros
		} else x1 = 0, v1 = r? q : 0;
		if (en >= r) ((uint8_t*)y)[r] = 0, u8[r] = r? q : 0;
		// loop fission: set scores first
		if (!(flag & KSW_EZ_GENERIC_SC)) {
			for (t = st0; t <= en0; t += 16) {
				__m128i sq, st, tmp, mask;
				sq = _mm_loadu_si128((__m128i*)&sf[t]);
				st = _mm_loadu_si128((__m128i*)&qrr[t]);
				mask = _mm_or_si128(_mm_cmpeq_epi8(sq, m1_), _mm_cmpeq_epi8(st, m1_));
				tmp = _mm_cmpeq_epi8(sq, st);
#ifdef __SSE4_1__
				tmp = _mm_blendv_epi8(sc_mis_, sc_mch_, tmp);
				tmp = _mm_blendv_epi8(tmp,     sc_N_,   mask);
#else
				tmp = _mm_or_si128(_mm_andnot_si128(tmp,  sc_mis_), _mm_and_si128(tmp,  sc_mch_));
				tmp = _mm_or_si128(_mm_andnot_si128(mask, tmp),     _mm_and_si128(mask, sc_N_));
#endif
				_mm_storeu_si128((__m128i*)((uint8_t*)s + t), tmp);
			}
		} else {
			for (t = st0; t <= en0; ++t)
				((uint8_t*)s)[t] = mat[sf[t] * m + qrr[t]];
		}
		// core loop
		x1_ = _mm_cvtsi32_si128(x1);
		v1_ = _mm_cvtsi32_si128(v1);
		st_ = st / 16, en_ = en / 16;
		assert(en_ - st_ + 1 <= n_col_);
		if (!with_cigar) { // score only
			for (t = st_; t <= en_; ++t) {
				__m128i z, a, b, xt1, vt1, ut, tmp, tmp2, tmp3;
				__dp_code_block1;
#ifdef __SSE4_1__
				z = _mm_max_epi8(z, a);                          // z = z > a? z : a (signed)
#else // we need to emulate SSE4.1 intrinsics _mm_max_epi8()
				z = _mm_and_si128(z, _mm_cmpgt_epi8(z, zero_));  // z = z > 0? z : 0;
				z = _mm_max_epu8(z, a);                          // z = max(z, a); this works because both are non-negative
#endif
				__dp_code_block2;
#ifdef __SSE4_1__
				_mm_store_si128(&x[t], _mm_max_epi8(a, zero_));
				_mm_store_si128(&y[t], _mm_max_epi8(b, zero_));
#else
				tmp = _mm_cmpgt_epi8(a, zero_);
				_mm_store_si128(&x[t], _mm_and_si128(a, tmp));
				tmp = _mm_cmpgt_epi8(b, zero_);
				_mm_store_si128(&y[t], _mm_and_si128(b, tmp));
#endif
			}
		} else if (!(flag&KSW_EZ_RIGHT)) { // gap left-alignment
			__m128i *pr = p + (size_t)r * n_col_ - st_;
			__m128i *pr2= p2+ (size_t)r * n_col_ - st_;
			off[r] = st, off_end[r] = en;
			for (t = st_; t <= en_; ++t) {
				__m128i z, a, b, xt1, vt1, ut, tmp, tmp2, tmp3;
				__dp_code_block1;
#ifdef __SSE4_1__
				z = _mm_max_epi8(z, a);                          // z = z > a? z : a (signed)
				tmp = _mm_cmpgt_epi8(b, z);
#else // we need to emulate SSE4.1 intrinsics _mm_max_epi8() and _mm_blendv_epi8()
				z = _mm_and_si128(z, _mm_cmpgt_epi8(z, zero_));  // z = z > 0? z : 0;
				z = _mm_max_epu8(z, a);                          // z = max(z, a); this works because both are non-negative
				tmp = _mm_cmpgt_epi8(b, z);
#endif
				__dp_code_block2;
				tmp = _mm_cmpgt_epi8(a, zero_);
				_mm_store_si128(&x[t], _mm_and_si128(tmp, a));
				tmp = _mm_cmpgt_epi8(b, zero_);
				_mm_store_si128(&y[t], _mm_and_si128(tmp, b));
				_mm_store_si128(&pr[t], tmp3); 
				_mm_store_si128(&pr2[t], tmp2);
			}
		} else { // gap right-alignment
			__m128i *pr = p + (size_t)r * n_col_ - st_;
			__m128i *pr2= p2+ (size_t)r * n_col_ - st_;
			off[r] = st, off_end[r] = en;
			for (t = st_; t <= en_; ++t) {
				__m128i z, a, b, xt1, vt1, ut, tmp, tmp2, tmp3;
				__dp_code_block1;
#ifdef __SSE4_1__
				z = _mm_max_epi8(z, a);                          // z = z > a? z : a (signed)
				tmp = _mm_cmpgt_epi8(z, b);
#else // we need to emulate SSE4.1 intrinsics _mm_max_epi8() and _mm_blendv_epi8()
				z = _mm_and_si128(z, _mm_cmpgt_epi8(z, zero_));  // z = z > 0? z : 0;
				z = _mm_max_epu8(z, a);                          // z = max(z, a); this works because both are non-negative
				tmp = _mm_cmpgt_epi8(z, b);
#endif
				__dp_code_block2;
				tmp = _mm_cmpgt_epi8(zero_, a);
				_mm_store_si128(&x[t], _mm_andnot_si128(tmp, a));
				tmp = _mm_cmpgt_epi8(zero_, b);
				_mm_store_si128(&y[t], _mm_andnot_si128(tmp, b));
				_mm_store_si128(&pr[t], tmp3); 
				_mm_store_si128(&pr2[t], tmp2);			
			}
		}
		if (!approx_max) { // find the exact max with a 32-bit score array
			int32_t max_H, max_t;
			// compute H[], max_H and max_t
			if (r > 0) {
				int32_t HH[4], tt[4], en1 = st0 + (en0 - st0) / 4 * 4, i;
				__m128i max_H_, max_t_, qe_;
				max_H = H[en0] = en0 > 0? H[en0-1] + u8[en0] - qe : H[en0] + v8[en0] - qe; // special casing the last element
				max_t = en0;
				max_H_ = _mm_set1_epi32(max_H);
				max_t_ = _mm_set1_epi32(max_t);
				qe_    = _mm_set1_epi32(q + e);
				for (t = st0; t < en1; t += 4) { // this implements: H[t]+=v8[t]-qe; if(H[t]>max_H) max_H=H[t],max_t=t;
					__m128i H1, tmp, t_;
					H1 = _mm_loadu_si128((__m128i*)&H[t]);
					t_ = _mm_setr_epi32(v8[t], v8[t+1], v8[t+2], v8[t+3]);
					H1 = _mm_add_epi32(H1, t_);
					H1 = _mm_sub_epi32(H1, qe_);
					_mm_storeu_si128((__m128i*)&H[t], H1);
					t_ = _mm_set1_epi32(t);
					tmp = _mm_cmpgt_epi32(H1, max_H_);
#ifdef __SSE4_1__
					max_H_ = _mm_blendv_epi8(max_H_, H1, tmp);
					max_t_ = _mm_blendv_epi8(max_t_, t_, tmp);
#else
					max_H_ = _mm_or_si128(_mm_and_si128(tmp, H1), _mm_andnot_si128(tmp, max_H_));
					max_t_ = _mm_or_si128(_mm_and_si128(tmp, t_), _mm_andnot_si128(tmp, max_t_));
#endif
				}
				_mm_storeu_si128((__m128i*)HH, max_H_);
				_mm_storeu_si128((__m128i*)tt, max_t_);
				for (i = 0; i < 4; ++i)
					if (max_H < HH[i]) max_H = HH[i], max_t = tt[i] + i;
				for (; t < en0; ++t) { // for the rest of values that haven't been computed with SSE
					H[t] += (int32_t)v8[t] - qe;
					if (H[t] > max_H)
						max_H = H[t], max_t = t;
				}
			} else H[0] = v8[0] - qe - qe, max_H = H[0], max_t = 0; // special casing r==0
			// update ez
			if (en0 == tlen - 1 && H[en0] > ez->mte)
				ez->mte = H[en0], ez->mte_q = r - en;
			if (r - st0 == qlen - 1 && H[st0] > ez->mqe)
				ez->mqe = H[st0], ez->mqe_t = st0;
			if (ksw_apply_zdrop(ez, 1, max_H, r, max_t, zdrop, e)) break;
			if (r == qlen + tlen - 2 && en0 == tlen - 1)
				ez->score = H[tlen - 1];
		} else { // find approximate max; Z-drop might be inaccurate, too.
			if (r > 0) {
				if (last_H0_t >= st0 && last_H0_t <= en0 && last_H0_t + 1 >= st0 && last_H0_t + 1 <= en0) {
					int32_t d0 = v8[last_H0_t] - qe;
					int32_t d1 = u8[last_H0_t + 1] - qe;
					if (d0 > d1) H0 += d0;
					else H0 += d1, ++last_H0_t;
				} else if (last_H0_t >= st0 && last_H0_t <= en0) {
					H0 += v8[last_H0_t] - qe;
				} else {
					++last_H0_t, H0 += u8[last_H0_t] - qe;
				}
				if ((flag & KSW_EZ_APPROX_DROP) && ksw_apply_zdrop(ez, 1, H0, r, last_H0_t, zdrop, e)) break;
			} else H0 = v8[0] - qe - qe, last_H0_t = 0;
			if (r == qlen + tlen - 2 && en0 == tlen - 1)
				ez->score = H0;
		}
		last_st = st, last_en = en;
		//for (t = st0; t <= en0; ++t) printf("(%d,%d)\t(%d,%d,%d,%d)\t%d\n", r, t, ((int8_t*)u)[t], ((int8_t*)v)[t], ((int8_t*)x)[t], ((int8_t*)y)[t], H[t]); // for debugging
	}
	kfree(km, mem);
	if (!approx_max) kfree(km, H);
	if (with_cigar) { // backtrack
		ksw2_singletrack_backtrace_affine(km, qlen, query, tlen, target, mat, q, e, n_col_*16, off, reinterpret_cast<int8_t*>(p), reinterpret_cast<int8_t*>(p2), ez);
		kfree(km, mem2); kfree(km, off);
	}
    #undef __dp_code_block1
	#undef __dp_code_block2
}
#endif // __SSE2__
