#include "structs_data.h"

#include <immintrin.h>

#include <limits.h>

// ----------------------------------------------------------------------------------------------------------
// Conversions
// -----------

/*__m256i RNS to int64_t RNS conversion using store, more efficient than extract.

BEFORE :
	- base contains the RNS base used to represent op (even if only base->size matters here)
	- op __m256i array to convert

AFTER :
	- rop contains the same values as op

NEEDS :
	- rop allocated

ENSURES :
	- op UNCHANGED
	- base UNCHANGED
*/
void from_m256i_to_int64_t_rns(int64_t *rop, struct rns_base_t *base, __m256i *op);

/*int64_t RNS to __m256i RNS conversion using store, more efficient than extract.

BEFORE :
	- base contains the RNS base used to represent op (even if only base->size matters here)
	- op int64_y array to convert

AFTER :
	- rop contains the same values as op

NEEDS :
	- rop allocated

ENSURES :
	- op UNCHANGED
	- base UNCHANGED
*/
void from_int64_t_to_m256i_rns(__m256i *rop, struct rns_base_t *base, int64_t *op);

void print_RNS(struct rns_base_t *base, int64_t *a);
void print_m256i(struct rns_base_t *base, __m256i *a);
void print_alone_m256i(__m256i a);

void avx_init_rns(struct rns_base_t *base);
void avx_initialize_inverses_base_conversion(struct conv_base_t *conv_base);

__m256i avx_add_mod_cr(__m256i a, __m256i b, __m256i k);
void avx_add_rns_cr(__m256i *rop, struct rns_base_t *base, __m256i *pa, __m256i *pb);

void avx_sub_rns_cr(__m256i *rop, struct rns_base_t *base, __m256i *pa, __m256i *pb);

void avx_add_aux_2e(__m256i *rop_up, __m256i *rop_lo, __m256i a, __m256i b);
void avx_add_aux_3e(__m256i *rop_up, __m256i *rop_lo, __m256i a, __m256i b, __m256i c);
void avx_mul_aux(__m256i *rop_up, __m256i *rop_lo, __m256i a, __m256i b);

__m256i avx_mul_mod_cr(__m256i a, __m256i b, __m256i k);
void avx_mul_rns_cr(__m256i *rop, struct rns_base_t *base, __m256i *pa, __m256i *pb);

//void avx_init_mrs(__m256i **rop,struct conv_base_t *conv_base);
void avx_init_mrs(struct conv_base_t *conv_base);
void avx_base_conversion_cr(__m256i *rop, struct conv_base_t *conv_base, __m256i *op, int64_t *a);
void avx_mult_mod_rns_cr(__m256i *rop, __m256i *pa, __m256i *pab, __m256i *pb,
						 __m256i *pbb, struct mod_mul_t *mult, __m256i *tmp0, __m256i *tmp1, __m256i *tmp2, int64_t *a);

// Cox functions
void avx_base_conversion_cox(__m256i *rop, struct conv_base_t *conv_base, __m256i *op, int r, int q, int alpha);
void avx_mult_mod_rns_cr_cox(__m256i *rop, __m256i *pa, __m256i *pab, __m256i *pb,
							 __m256i *pbb, struct mod_mul_t *mult, __m256i *tmp0, __m256i *tmp1, __m256i *tmp2, int64_t *a);
