#include "structs_data.h"
#include <immintrin.h>

// RNS arithmetic functions
void add_rns(int64_t *rop, struct rns_base_t  *base, int64_t *pa, int64_t *pb);
void sub_rns(int64_t *rop, struct rns_base_t  *base, int64_t *pa, int64_t *pb);
void mul_rns(int64_t *rop, struct rns_base_t  *base, int64_t *pa, int64_t *pb);
// void mult_mod_rns(int64_t *rop, int64_t *pa, int64_t *pb, struct mod_mul_t *mult, int64_t *tmp[3]);
void mult_mod_rns(int64_t *rop, int64_t *pa, int64_t *pab, int64_t *pb, 
	int64_t *pbb, struct mod_mul_t *mult, int64_t *tmp[3]);

void init_rns(struct rns_base_t *base);
void from_int_to_rns(int64_t *rop, struct rns_base_t *base, mpz_t op);
void from_rns_to_int_crt(mpz_t rop, struct rns_base_t *base, int64_t *op);
void initialize_inverses_base_conversion(struct conv_base_t *conv_base);
void base_conversion(int64_t *rop, struct conv_base_t *conv_base, int64_t *op);

void base_conversion_cr(int64_t *rop, struct conv_base_t *conv_base, int64_t *op, int64_t *a);
void add_rns_cr(int64_t *rop, struct rns_base_t  *base, int64_t *pa, int64_t *pb);
void sub_rns_cr(int64_t *rop, struct rns_base_t  *base, int64_t *pa, int64_t *pb);
void mul_rns_cr(int64_t *rop, struct rns_base_t  *base, int64_t *pa, int64_t *pb);
// void mult_mod_rns_cr(int64_t *rop, int64_t *pa, int64_t *pb, struct mod_mul_t *mult, int64_t *tmp[4]);
void mult_mod_rns_cr(int64_t *rop, int64_t *pa, int64_t *pab, int64_t *pb, 
	int64_t *pbb, struct mod_mul_t *mult, int64_t *tmp[4]);

int64_t add_mod_cr(int64_t a, int64_t b, int k) ;
int64_t mul_mod_cr(int64_t a, int64_t b, int k) ; 

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//

void from_m256i_to_rns(int64_t *rop, struct rns_base_t *base, __m256i *op);
void from_rns_to_m256i(__m256i *rop, struct rns_base_t *base, int64_t *op);

void print_RNS(struct rns_base_t *base, int64_t *a);
void print_m256i(struct rns_base_t *base, __m256i *a);
void print_alone_m256i(__m256i a);

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
int avx_compute_k_cox(__m256i *op, struct rns_base_t *base, int r, int q, int alpha);
void avx_base_conversion_cox(__m256i *rop, struct conv_base_t *conv_base, __m256i *op, int r, int q, int alpha);
void avx_mult_mod_rns_cr_cox(__m256i *rop, __m256i *pa, __m256i *pab, __m256i *pb, 
	__m256i *pbb, struct mod_mul_t *mult, __m256i *tmp[4]);