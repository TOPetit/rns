#include <stdlib.h>

#include <stdio.h>

#include <gmp.h>

#include <math.h>

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

    mpz_t M;
    mpz_inits(M, NULL);
    mpz_set(M, rns_a.M); // Get M from the base

    printf("M = %Zd\n", M);

    /////////////////////////////
    // TEST CONVERSION
    /////////////////////////////

    mpz_urandomm(A, r_state, rns_a.M);
    from_int_to_rns(op1, &rns_a, A);

    test = true;

    for (i = 0; i < rns_a.size; ++i)
    {
        test = test && (mpz_get_ui(A) % rns_a.m[i] == op1[i]);
        printf("%Zd and %lu\n", M, op1[i]);
    }

    printf("Conversion from int to rns... ");
    if (test)
        printf("OK\n ");
    else
        printf("ERROR\n");

    gmp_randclear(r_state);
    mpz_clear(A);
    mpz_clear(B);
    mpz_clear(C);

    return 0;
}
