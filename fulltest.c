#include <stdlib.h>

#include <stdio.h>

#include <gmp.h>

#include <math.h>

#include <time.h>

#include "rns.h"

#include "tests.c"

#include "rnsv.h"

#define bool unsigned int

#define true 1

#define false 0

int main(void)
{

    /////////////////////////////
    // INIT
    /////////////////////////////

    // Init Random
    mpz_t rand_Num;
    unsigned long int i, seed;
    gmp_randstate_t r_state;

    seed = clock();

    gmp_randinit_default(r_state);
    gmp_randseed_ui(r_state, seed);

    // Test variables
    bool test;

    // Variables
    int64_t op1[NB_COEFF];
    int64_t op2[NB_COEFF];
    int64_t res[NB_COEFF];

    mpz_t A, B, C;
    mpz_inits(A, B, C, NULL);

    // Base
    struct rns_base_t rns_a;
    rns_a.size = NB_COEFF;

    int64_t m1[NB_COEFF] = {
        9223372036854775805,
        9223372036854775801,
        9223372036854775789,
        9223372036854775783,
        9223372036854775777,
        9223372036854775769,
        9223372036854775757,
        9223372036854775747};
    rns_a.m = m1;

    int k1[NB_COEFF] = {
        3,
        7,
        19,
        25,
        31,
        39,
        51,
        61};
    rns_a.k = k1;

    init_rns(&rns_a);

    int64_t tmp_k[NB_COEFF];

    for (int j = 0; j < NB_COEFF; j++)
    {
        tmp_k[j] = (int64_t)k1[j];
    }

    __m256i avx_k1[NB_COEFF / 4];
    from_rns_to_m256i(avx_k1, &rns_a, tmp_k);
    rns_a.avx_k = avx_k1;

    // Second Base
    struct rns_base_t rns_b;
    rns_b.size = NB_COEFF;

    int64_t m2[NB_COEFF] = {
        9223372036854775807,
        9223372036854775803,
        9223372036854775799,
        9223372036854775787,
        9223372036854775781,
        9223372036854775771,
        9223372036854775763,
        9223372036854775753};
    rns_b.m = m2;

    int k2[NB_COEFF] = {
        1,
        5,
        9,
        21,
        27,
        37,
        45,
        55};
    rns_b.k = k2;

    init_rns(&rns_b);

    for (int j = 0; j < NB_COEFF; j++)
    {
        tmp_k[j] = (int64_t)k2[j];
    }
    __m256i avx_k2[NB_COEFF / 4];
    from_rns_to_m256i(avx_k2, &rns_b, tmp_k);
    rns_b.avx_k = avx_k2;

    /*
    // Base
    struct rns_base_t rns_a;
    rns_a.size = NB_COEFF;

    int64_t m_tmp[NB_COEFF] = {
        32767,
        32763,
        32759,
        32747,
        32741,
        32731,
        32723,
        32717};
    rns_a.m = m_tmp;

    int k_tmp[NB_COEFF] = {
        1,
        5,
        9,
        21,
        27,
        37,
        45,
        51};
    rns_a.k = k_tmp;

    init_rns(&rns_a);
    
    struct rns_base_t rns_b;
    rns_b.size = NB_COEFF;

    int64_t m_tmp_bis[NB_COEFF] = {
        32765,
        32761,
        32749,
        32743,
        32737,
        32729,
        32719,
        32713};
    rns_b.m = m_tmp_bis;

    int k_tmp_bis[NB_COEFF] = {
        3,
        7,
        19,
        25,
        31,
        37,
        45,
        55};
    rns_b.k = k_tmp_bis;

    init_rns(&rns_b);
    */

    gmp_printf("M = %Zd\n", rns_a.M);

    /////////////////////////////
    // TEST CONVERSION INT -> RNS
    /////////////////////////////

    mpz_urandomm(A, r_state, rns_a.M);
    from_int_to_rns(op1, &rns_a, A);

    mpz_t R;
    mpz_inits(R, NULL);
    long unsigned int r;

    test = true;

    for (i = 0; i < rns_a.size; ++i)
    {
        r = mpz_fdiv_r_ui(R, A, rns_a.m[i]);
        //printf("%lu and %lu\n", r, op1[i]);
    }

    mpz_clear(R);

    printf("Conversion from int to RNS... ");
    if (test)
        printf("OK\n");
    else
        printf("ERROR\n");

    /////////////////////////////
    // TEST CONVERSION RNS -> INT
    /////////////////////////////
    from_rns_to_int_crt(B, &rns_a, op1);

    //gmp_printf("%Zd\n%Zd\n", A, B);

    printf("Conversion from RNS to int... ");
    if (mpz_cmp(A, B) == 0)
        printf("OK\n");
    else
        printf("ERROR\n");

    /////////////////////////////
    // TEST SEQUENTIAL ADDITION
    /////////////////////////////

    mpz_t D;
    mpz_inits(D, NULL);

    mpz_urandomm(A, r_state, rns_a.M);
    mpz_urandomm(B, r_state, rns_a.M);
    from_int_to_rns(op1, &rns_a, A);
    from_int_to_rns(op2, &rns_a, B);

    add_rns_cr(res, &rns_a, op1, op2);

    from_rns_to_int_crt(D, &rns_a, res);

    mpz_add(C, A, B);
    from_int_to_rns(op2, &rns_a, C);

    printf("Int64_t RNS addition... ");
    if (rns_equal(rns_a, res, op2))
        printf("OK\n");
    else
        printf("ERROR\n");

    /////////////////////////////
    // TEST SEQUENTIAL SUBSTRACTION
    /////////////////////////////

    mpz_urandomm(A, r_state, rns_a.M);
    mpz_urandomm(B, r_state, rns_a.M);
    from_int_to_rns(op1, &rns_a, A);
    from_int_to_rns(op2, &rns_a, B);

    sub_rns_cr(res, &rns_a, op1, op2);

    from_rns_to_int_crt(D, &rns_a, res);

    mpz_sub(C, A, B);
    from_int_to_rns(op2, &rns_a, C);

    printf("Int64_t RNS substraction... ");
    if (rns_equal(rns_a, res, op2))
        printf("OK\n");
    else
        printf("ERROR\n");

    /////////////////////////////
    // TEST SEQUENTIAL MULTIPLICATION
    /////////////////////////////

    mpz_urandomm(A, r_state, rns_a.M);
    mpz_urandomm(B, r_state, rns_a.M);
    from_int_to_rns(op1, &rns_a, A);
    from_int_to_rns(op2, &rns_a, B);

    mul_rns_cr(res, &rns_a, op1, op2);

    from_rns_to_int_crt(D, &rns_a, res);

    mpz_mul(C, A, B);
    from_int_to_rns(op2, &rns_a, C);

    printf("Int64_t RNS multipliation... ");
    if (rns_equal(rns_a, res, op2))
        printf("OK\n");
    else
        printf("ERROR\n");

    /////////////////////////////
    // TEST SEQUENTIAL BASE CONVERSION
    /////////////////////////////

    struct conv_base_t conv;
    conv.rns_a = &rns_a;
    conv.rns_b = &rns_b;
    initialize_inverses_base_conversion(&conv);

    mpz_urandomm(A, r_state, rns_a.M);
    from_int_to_rns(op1, &rns_a, A);

    int64_t a[NB_COEFF];
    base_conversion_cr(op2, &conv, op1, a);
    from_rns_to_int_crt(A, &rns_a, op1);
    from_rns_to_int_crt(B, &rns_b, op2);

    printf("Int64_t RNS base conversion... ");
    if (mpz_cmp(A, B) == 0)
        printf("OK\n");
    else
        printf("ERROR\n");

    /////////////////////////////
    // TEST CONVERSION RNS -> AVX-2
    /////////////////////////////

    __m256i avx_op1[NB_COEFF / 8];

    mpz_urandomm(A, r_state, rns_a.M);
    from_int_to_rns(op1, &rns_a, A);
    from_rns_to_m256i(avx_op1, &rns_a, op1);

    test = true;
    for (int i = 0; i < rns_a.size / 4; i++)
    {
        test = test && (_mm256_extract_epi64(avx_op1[i], 0) == op1[4 * i]);
        test = test && (_mm256_extract_epi64(avx_op1[i], 1) == op1[4 * i + 1]);
        test = test && (_mm256_extract_epi64(avx_op1[i], 2) == op1[4 * i + 2]);
        test = test && (_mm256_extract_epi64(avx_op1[i], 3) == op1[4 * i + 3]);
    }

    printf("Conversion from RNS to AVX-2... ");
    if (test)
        printf("OK\n");
    else
        printf("ERROR\n");

    /////////////////////////////
    // TEST CONVERSION AVX-2 -> RNS
    /////////////////////////////

    from_m256i_to_rns(op2, &rns_a, avx_op1);

    printf("Conversion from AVX-2 to RNS... ");
    if (rns_equal(rns_a, op1, op2))
        printf("OK\n");
    else
        printf("ERROR\n");

    /////////////////////////////
    // TEST PARALLEL ADDITION
    /////////////////////////////

    __m256i avx_op2[NB_COEFF / 4];
    __m256i avx_res[NB_COEFF / 4];

    //mpz_urandomm(A, r_state, rns_a.M);
    //mpz_urandomm(B, r_state, rns_a.M);
    from_int_to_rns(op1, &rns_a, A);
    from_int_to_rns(op2, &rns_a, B);
    from_rns_to_m256i(avx_op1, &rns_a, op1);
    from_rns_to_m256i(avx_op2, &rns_a, op2);

    avx_add_rns_cr(avx_res, &rns_a, avx_op1, avx_op2);
    add_rns_cr(res, &rns_a, op1, op2);

    from_m256i_to_rns(op1, &rns_a, avx_res);

    printf("AVX-2 RNS addition... ");
    if (rns_equal(rns_a, op1, res))
        printf("OK\n");
    else
        printf("ERROR\n");

    /////////////////////////////
    // TEST PARALLEL SUBSTRACTION
    /////////////////////////////

    //mpz_urandomm(A, r_state, rns_a.M);
    //mpz_urandomm(B, r_state, rns_a.M);
    from_int_to_rns(op1, &rns_a, A);
    from_int_to_rns(op2, &rns_a, B);
    from_rns_to_m256i(avx_op1, &rns_a, op1);
    from_rns_to_m256i(avx_op2, &rns_a, op2);

    avx_sub_rns_cr(avx_res, &rns_a, avx_op1, avx_op2);
    sub_rns_cr(res, &rns_a, op1, op2);

    from_m256i_to_rns(op1, &rns_a, avx_res);

    for (int i = 0; i < rns_a.size; i++)
    {
        printf("%ld\n%ld\n\n", op1[i], res[i]);
    }

    printf("AVX-2 RNS substraction... ");
    if (rns_equal(rns_a, op1, res))
        printf("OK\n");
    else
        printf("ERROR\n");

    /////////////////////////////
    // TEST PARALLEL MULTIPLICATION
    /////////////////////////////

    //mpz_urandomm(A, r_state, rns_a.M);
    //mpz_urandomm(B, r_state, rns_a.M);
    from_int_to_rns(op1, &rns_a, A);
    from_int_to_rns(op2, &rns_a, B);
    from_rns_to_m256i(avx_op1, &rns_a, op1);
    from_rns_to_m256i(avx_op2, &rns_a, op2);

    avx_mul_rns_cr(avx_res, &rns_a, avx_op1, avx_op2);
    mul_rns_cr(res, &rns_a, op1, op2);

    from_m256i_to_rns(op1, &rns_a, avx_res);

    printf("AVX-2 RNS multiplication... ");
    if (rns_equal(rns_a, op1, res))
        printf("OK\n");
    else
        printf("ERROR\n");

    return 0;
}
