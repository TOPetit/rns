#include <stdlib.h>

#include <stdio.h>

#include <gmp.h>

#include <math.h>

#include "rns.h"

#include "tests.c"

#include "rnsv.h"

int main(void)
{

    // Timing
    unsigned long long before_cycles, after_cycles = ULLONG_MAX;

    // Initializing random
    gmp_randstate_t state;
    gmp_randinit_default(state);

    int64_t op1[NB_COEFF];
    int64_t op2[NB_COEFF];
    __m256i avx_op1[NB_COEFF / 4];
    __m256i avx_op2[NB_COEFF / 4];

    mpz_t A, B;
    mpz_inits(A, B, NULL);

    // Base
    struct rns_base_t rns_a;
    rns_a.size = NB_COEFF;

    int64_t m_tmp[NB_COEFF] = {
        9223372036854775805,
        9223372036854775801,
        9223372036854775789,
        9223372036854775783,
        9223372036854775777,
        9223372036854775769,
        9223372036854775757,
        9223372036854775747};
    rns_a.m = m_tmp;

    int k_tmp[NB_COEFF] = {
        3,
        7,
        19,
        25,
        31,
        39,
        51,
        61};
    rns_a.k = k_tmp;

    init_rns(&rns_a);

    int64_t tmp_k[NB_COEFF];
    int j;
    for (j = 0; j < NB_COEFF; j++)
    {
        tmp_k[j] = (int64_t)k_tmp[j];
    }

    __m256i avx_k[NB_COEFF / 4];
    from_int64_t_to_m256i_rns(avx_k, &rns_a, tmp_k);
    rns_a.avx_k = avx_k;

    mpz_t M;
    mpz_inits(M, NULL);
    mpz_set(M, rns_a.M); // Get M from the base

    mpz_set_str(A, "100106136745050507346481674824002435247796236765700289941745547607771451581114", 10);
    gmp_printf("A before : %Zd\n", A);

    from_int_to_rns(op1, &rns_a, A);
    print_RNS(&rns_a, op1);
    printf("\n");
    from_int64_t_to_m256i_rns(avx_op1, &rns_a, op1);
    from_m256i_to_int64_t_rns(op2, &rns_a, avx_op1);
    from_m256i_to_int64_t_rns(op2, &rns_a, avx_op1);
    from_m256i_to_int64_t_rns(op2, &rns_a, avx_op1);
    from_m256i_to_int64_t_rns(op2, &rns_a, avx_op1);
    from_m256i_to_int64_t_rns(op2, &rns_a, avx_op1);
    from_m256i_to_int64_t_rns(op2, &rns_a, avx_op1);
    from_m256i_to_int64_t_rns(op2, &rns_a, avx_op1);
    from_m256i_to_int64_t_rns(op2, &rns_a, avx_op1);
    before_cycles = cpucyclesStart();
    from_m256i_to_int64_t_rns(op2, &rns_a, avx_op1);
    after_cycles = cpucyclesStop();
    printf("\n Cycles without extract : %lld\n", after_cycles - before_cycles);
    from_rns_to_int_crt(B, &rns_a, op2);
    gmp_printf("already built functions.\nA after  : %Zd\n", B);

    printf("tests.\n");
    mpz_set_str(A, "100106136745050507346481674824002435247796236765700289941745547607771451581114", 10);

    from_int_to_rns(op1, &rns_a, A);
    from_int64_t_to_m256i_rns(avx_op1, &rns_a, op1);

    from_m256i_to_int64_t_rns(op2, &rns_a, avx_op1);
    from_m256i_to_int64_t_rns(op2, &rns_a, avx_op1);
    from_m256i_to_int64_t_rns(op2, &rns_a, avx_op1);
    from_m256i_to_int64_t_rns(op2, &rns_a, avx_op1);
    from_m256i_to_int64_t_rns(op2, &rns_a, avx_op1);
    from_m256i_to_int64_t_rns(op2, &rns_a, avx_op1);
    from_m256i_to_int64_t_rns(op2, &rns_a, avx_op1);
    from_m256i_to_int64_t_rns(op2, &rns_a, avx_op1);
    before_cycles = cpucyclesStart();
    from_m256i_to_int64_t_rns(op2, &rns_a, avx_op1);
    after_cycles = cpucyclesStop();
    printf("\n Cycles with extract : %lld\n", after_cycles - before_cycles);

    from_rns_to_int_crt(B, &rns_a, op2);
    gmp_printf("\nA after  : %Zd\n", B);
}