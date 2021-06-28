#include <stdlib.h>

#include <stdio.h>

#include <gmp.h>

#include <math.h>

#include "rns.h"

#include "tests.c"

int main(void)
{
    int64_t a, b, u, v, p, mp, r;
    int64_t g;
    int64_t result;
    int64_t r_dw;

    a = 2147483648;
    b = 214748364;
    // p = 18446744073709551599; // 2^64-17 thus -p mod 2^64 = 17
    p = 4294967279; // 2^32-17 = 2^32-2^4-1
    r = 4042322161; // (-p^{-1} mod 2^n)

    // Initializing random
    gmp_randstate_t state;
    gmp_randinit_default(state);
    // Timers
    unsigned long long timer, t1, t2;
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

    int64_t tmp_k[NB_COEFF];
    int j;
    for (j = 0; j < NB_COEFF; j++)
    {
        tmp_k[j] = (int64_t)k_tmp[j];
    }
}