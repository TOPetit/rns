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
    unsigned long int i, seed;
    gmp_randstate_t r_state;

    //seed = clock();

    gmp_randinit_default(r_state);
    //gmp_randseed_ui(r_state, seed);

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
    avx_init_rns(&rns_a);

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
    avx_init_rns(&rns_b);

    for (int j = 0; j < NB_COEFF; j++)
    {
        tmp_k[j] = (int64_t)k2[j];
    }
    __m256i avx_k2[NB_COEFF / 4];
    from_rns_to_m256i(avx_k2, &rns_b, tmp_k);
    rns_b.avx_k = avx_k2;

    /////////////////////////////
    // TEST SEQUENTIAL ADDITION
    /////////////////////////////

    mpz_t D;
    mpz_inits(D, NULL);

    

    /////////////////////////////
    // TEST SEQUENTIAL MULTIPLICATION
    /////////////////////////////

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
    // TEST SEQUENTIAL MOD MULT
    /////////////////////////////

    mpz_t inv_p_modM, inv_M_modMp, modul_p;
    mpz_inits(inv_p_modM, inv_M_modMp, modul_p, NULL);

    int64_t pa[NB_COEFF];
    int64_t pb[NB_COEFF];
    int64_t pab[NB_COEFF];
    int64_t pbb[NB_COEFF];
    int64_t pc[NB_COEFF];
    int64_t pp1[NB_COEFF];
    int64_t pp2[NB_COEFF];
    int64_t pp3[NB_COEFF];

    // Set custom values
    mpz_set_str(modul_p, "115792089021636622262124715160334756877804245386980633020041035952359812890593", 10);
    mpz_set_str(inv_p_modM, "-7210642370083763919688086698199040857322895088554003933210287226647459666846134833419938084604981461493089686639677942359747717700454441525223348684285", 10);
    mpz_set_str(inv_M_modMp, "2926906825829426928727294150364906856635623568440932569450673109926460590684432927230290255276608760237299661987870702836538185953568700154975953006659", 10);

    //Modular multiplication

    struct mod_mul_t mult;
    mpz_t tmp_gcd, t, tmp_inv;

    mpz_init(tmp_gcd);
    mpz_init(t);
    mpz_init(tmp_inv);
    from_int_to_rns(pp2, &rns_b, modul_p); // P mod Mb

    mpz_sub(tmp_inv, rns_a.M, modul_p);
    mpz_gcdext(tmp_gcd, inv_p_modM, t, tmp_inv, rns_a.M);
    from_int_to_rns(pp1, &rns_a, inv_p_modM); //(-P)^-1 mod Ma

    mpz_gcdext(tmp_gcd, inv_M_modMp, t, rns_a.M, rns_b.M);
    from_int_to_rns(pp3, &rns_b, inv_M_modMp); // Ma^{-1} mod Mb

    mult.inv_p_modMa = pp1;
    mult.p_modMb = pp2;
    mult.inv_Ma_modMb = pp3;
    mult.conv = &conv;

    /////////////////////////////
    // TEST CONVERSION RNS -> AVX-2
    /////////////////////////////

    
    /////////////////////////////
    // TEST PARALLEL ADDITION
    /////////////////////////////

    
    __m256i avx_op1[NB_COEFF / 4];
    __m256i avx_op2[NB_COEFF / 4];
    __m256i avx_res[NB_COEFF / 4];

    
    /////////////////////////////
    // TEST PARALLEL BASE CONVERSION
    /////////////////////////////

    avx_init_mrs(&conv);
    avx_initialize_inverses_base_conversion(&conv);

    from_int_to_rns(op1, &rns_a, A);
    from_rns_to_m256i(avx_op1, &rns_a, op1);

    avx_base_conversion_cr(avx_op2, &conv, avx_op1, a);
    from_m256i_to_rns(op1, &rns_a, avx_op2);
    from_int_to_rns(op2, &rns_b, A);

    printf("AVX-2 RNS base conversion... ");
    if (rns_equal(rns_b, op1, op2))
        printf("OK\n");
    else
        printf("ERROR\n");

    
    return 0;
}
