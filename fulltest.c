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

    printf("Conversion from int to rns... ");
    if (test)
        printf("OK\n");
    else
        printf("ERROR\n");

    /////////////////////////////
    // TEST CONVERSION RNS -> INT
    /////////////////////////////
    from_rns_to_int_crt(B, &rns_a, op1);

    //gmp_printf("%Zd\n%Zd\n", A, B);

    printf("Conversion from rns to int... ");
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

    gmp_printf("%Zd\n%Zd\n", C, D);

    printf("Sequential addition... ");
    if (mpz_cmp(C, D) == 0)
        printf("OK\n");
    else
        printf("ERROR\n");

    gmp_randclear(r_state);
    mpz_clear(A);
    mpz_clear(B);
    mpz_clear(C);

    return 0;
}
