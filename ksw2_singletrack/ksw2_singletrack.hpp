#include <string.h>
#include <stdio.h>
#include <assert.h>
#include "ksw2.h"

void ksw2_singletrack_backtrace_affine2p(void* km, int qlen, const uint8_t *query, int tlen, const uint8_t *target, const int8_t *mat,
				   int8_t q, int8_t e, int8_t q2, int8_t e2, const int n_col, const int *off, int8_t *p, int8_t *p2, ksw_extz_t *ez);

void ksw2_singletrack_backtrace_affine(void* km, int qlen, const uint8_t *query, int tlen, const uint8_t *target, const int8_t *mat,
				   int8_t q, int8_t e, const int n_col, const int *off, int8_t *p, int8_t *p2, ksw_extz_t *ez);

#ifdef __SSE2__
void ksw_extd2_singletrack_sse(void *km, int qlen, const uint8_t *query, int tlen, const uint8_t *target, int8_t m, const int8_t *mat,
				   int8_t q, int8_t e, int8_t q2, int8_t e2, int w, int zdrop, int end_bonus, int flag, ksw_extz_t *ez);

void ksw_extz2_singletrack_sse(void *km, int qlen, const uint8_t *query, int tlen, const uint8_t *target, int8_t m, const int8_t *mat, 
                    int8_t q, int8_t e, int w, int zdrop, int end_bonus, int flag, ksw_extz_t *ez);
#endif 
