#include <stdlib.h>
#include <stdio.h>
#include <gmp.h>

#include <math.h>

#include "rns.h"
#include "tests.c"

#include "rnsv.h"

int main(void) {

    // Initializing random
	gmp_randstate_t state;
	gmp_randinit_default (state);

    int64_t op1[NB_COEFF];
	__m256i avx_op1[NB_COEFF/4];


	mpz_t A;
	mpz_inits(A, NULL);
    int64_t tmp;

	// Base
	struct rns_base_t rns_a;
	rns_a.size = NB_COEFF;
	
	int64_t m_tmp[NB_COEFF] = {9223372036854775805,
9223372036854775801,
9223372036854775789,
9223372036854775783,
9223372036854775777,
9223372036854775769,
9223372036854775757,
9223372036854775747};
	rns_a.m = m_tmp;
	

	int k_tmp[NB_COEFF] = {3, 7, 19, 25, 31, 39, 51, 61};
	rns_a.k = k_tmp;

	init_rns(&rns_a);

	int64_t tmp_k[NB_COEFF];
    int j;
    for (j=0;j<NB_COEFF;j++){
    	tmp_k[j]=(int64_t)k_tmp[j];
    }

    __m256i avx_k[NB_COEFF/4];
    from_rns_to_m256i(avx_k,&rns_a,tmp_k);
    rns_a.avx_k = avx_k;

/*
	__m256i avx_inv_Mi[NB_COEFF/4];
	from_rns_to_m256i(avx_inv_Mi, &rns_a, rns_a.inv_Mi);
	rns_a.avx_inv_Mi = avx_inv_Mi;
*/
	
	mpz_t M;
	mpz_inits(M, NULL);
	mpz_set(M, rns_a.M); // Get M from the base


    mpz_urandomm(A, state, M); // random A
    from_int_to_rns(op1, &rns_a, A);

    printf("A in RNS : ");
    print_RNS(&rns_a, op1);
    printf("\n");

    from_rns_to_int_crt(tmp, &rns_a, op1);
    printf("A before : %Zd\n", tmp);




}