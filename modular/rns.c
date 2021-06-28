#include "rns.h"

#include <stdlib.h>

#include <stdio.h>

///////////////////////////////
// RNS addition
///////////////////////////////
// rop : result
// base : RNS base
// pa : A
// pb : B
inline void add_rns(int64_t *rop, struct rns_base_t *base, int64_t *pa, int64_t *pb)
{
    int j;
    int128 tmp;

    for (j = 0; j < base->size; j++)
    {
        tmp = (int128)pa[j] + pb[j];
        rop[j] = (int64_t)(tmp % base->m[j]);
    }
}

///////////////////////////////
// RNS substraction
///////////////////////////////
// rop : result
// base : RNS base
// pa : A
// pb : B
inline void sub_rns(int64_t *rop, struct rns_base_t *base, int64_t *pa, int64_t *pb)
{
    int j;
    int128 tmp;

    for (j = 0; j < base->size; j++)
    {
        tmp = (int128)pa[j] - pb[j];
        rop[j] = (int64_t)(tmp % base->m[j]);
    }
}

///////////////////////////////
// RNS multiplication
///////////////////////////////
// rop : result
// base : RNS base
// pa : A
// pb : B
inline void mul_rns(int64_t *rop, struct rns_base_t *base, int64_t *pa, int64_t *pb, int64_t n)
{
    int j;
    int128 tmp, tmp1;

    for (j = 0; j < base->size; j++)
    {
        tmp = (int128)pa[j] * pb[j];
        tmp1 = (int128)tmp * 1 << (2 * n);
        rop[j] = (int64_t)(tmp % base->m[j]);
    }
}

static inline int64_t mul_mod_rns(int64_t *rop, struct rns_base_t *base, int64_t *pa, int64_t *pb, int64_t n)
{
    int64_t c;
    int64_t tmp1_w;
    int64_t tmp2_w;
    int W = 2 * n;
    int size = base->size;

    for (int i = 0; i < size; i++)
    {
        c = (int64_t)pa[i] * pb[i];
        tmp1_w = ((c * base->r[i]) >> W) + 1;
        tmp2_w = ((int64_t)tmp1_w * base->m[i]) >> W;
        if (tmp2_w == base->m[i])
            rop[i] = 0;
        else
            rop[i] = tmp2_w;
    }
}

///////////////////////////////
// GMP to RNS convertion
///////////////////////////////
// TODO : NB_COEFF to NB_RES
//~ Assumes allocation already done for "rop".
void from_int_to_rns(int64_t *rop, struct rns_base_t *base, mpz_t op)
{
    int i;
    mpz_t tmp_residue;

    //	printf("On rentre ds la conversion\n");

    mpz_init(tmp_residue);
    if (op->_mp_size == 0)
        return;

    //	printf("avant la boucle \n");

    for (i = 0; i < base->size; i++)
    {
        //		printf(" i %d\n",i);

        mpz_fdiv_r_ui(tmp_residue, op, base->m[i]);
        rop[i] = mpz_get_ui(tmp_residue);
    }

    //printf("On sorte de la conversion \n");
}

/////////////////////////////////
// Initialize some space for
// the constants of the RNS base
// Computes M
// Computes Mi and inv_Mi
/////////////////////////////////
// TODO : NB_COEFF to NB_RES
void init_rns(struct rns_base_t *base)
{
    mpz_t *Mi;
    mpz_t *inv_Mi;
    int i;
    mpz_t tmp_gcd, tmp_mi, tmp, t;
    int64_t *int_inv_Mi;
    int64_t up;

    Mi = (mpz_t *)malloc(base->size * sizeof(mpz_t));
    inv_Mi = (mpz_t *)malloc(base->size * sizeof(mpz_t));
    int_inv_Mi = (int64_t *)malloc(base->size * sizeof(int64_t));
    for (i = 0; i < base->size; i++)
    {
        mpz_init(Mi[i]);
        mpz_init(inv_Mi[i]);
    }
    mpz_init(base->M);
    mpz_init(tmp_gcd);
    mpz_init(tmp_mi);
    mpz_init(tmp);
    mpz_init(t);

    // Computes M
    mpz_add_ui(base->M, base->M, 1);
    for (i = 0; i < base->size; i++)
    {
        mpz_mul_ui(base->M, base->M, base->m[i]);
    }
    // Computes Mi and inv_Mi
    for (i = 0; i < base->size; i++)
    {
        mpz_fdiv_q_ui(Mi[i], base->M, base->m[i]);
        mpz_set_ui(tmp_mi, base->m[i]);
        mpz_gcdext(tmp_gcd, inv_Mi[i], t, Mi[i], tmp_mi);
    }
    base->Mi = Mi;
    base->inv_Mi = inv_Mi;
    // Converts inv_Mi in RNS ie just Inv_Mi mod m_i
    mpz_clears(tmp_gcd, tmp_mi, tmp, t, NULL);
}
