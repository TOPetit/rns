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

inline void add_rns_cr(int64_t *rop, struct rns_base_t *base, int64_t *pa, int64_t *pb)
{
	int j;
	// int128 tmp;
	// int128 up;
	//int64_t tmp, up, lo, mask = ((int64_t)1<<63) -1;  /////////////////////////////////////

	for (j = 0; j < base->size; j++)
	{
		rop[j] = add_mod_cr(pa[j], pb[j], base->k[j]);
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

inline void sub_rns_cr(int64_t *rop, struct rns_base_t *base, int64_t *pa, int64_t *pb)
{
	int j;
	int128 tmp;
	int128 up;

	for (j = 0; j < base->size; j++)
	{
		tmp = (int128)pa[j] - pb[j];
		up = tmp >> 63; //////////////////////////////////////////////////////

		rop[j] = (int64_t)(tmp + up * base->k[j]);
	}
}

///////////////////////////////
// RNS multiplication
///////////////////////////////
// rop : result
// base : RNS base
// pa : A
// pb : B
inline void mul_rns(int64_t *rop, struct rns_base_t *base, int64_t *pa, int64_t *pb)
{
	int j;
	int128 tmp;

	for (j = 0; j < base->size; j++)
	{
		tmp = (int128)pa[j] * pb[j];
		rop[j] = (int64_t)(tmp % base->m[j]);
	}
}

inline void mul_rns_cr(int64_t *rop, struct rns_base_t *base, int64_t *pa, int64_t *pb)
{
	int j;

	for (j = 0; j < base->size; j++)
	{
		rop[j] = mul_mod_cr(pa[j], pb[j], base->k[j]);
	}
}

///////////////////////////////
// RNS Modular multiplication
///////////////////////////////
// rop : result
// pa : A
// pb : B
// mult : constants
// tmp : temporary arrays for intermediate results
void mult_mod_rns(int64_t *rop, int64_t *pa, int64_t *pab, int64_t *pb,
				  int64_t *pbb, struct mod_mul_t *mult, int64_t *tmp[3])
{

	int i;

	// mpz_t t;      // For debugging
	// mpz_init(t);  // For debugging

	mul_rns(tmp[0], mult->conv->rns_a, pa, pb); //A*B

	// from_rns_to_int_crt(t, mult->conv->rns_a, tmp[0]);
	// gmp_printf("A*B: %Zd\n", t);

	mul_rns(tmp[1], mult->conv->rns_b, pab, pbb); //A*B in b&se2

	// base_conversion(tmp[1], mult->conv, tmp[0]); // A*B in base2 !! NOT GOOD

	// from_rns_to_int_crt(t, mult->conv->rns_b, tmp[1]);
	// gmp_printf("A*B in base2: %Zd\n", t);

	mul_rns(tmp[2], mult->conv->rns_a, tmp[0], mult->inv_p_modMa); //Q*{P-1}

	// from_rns_to_int_crt(t, mult->conv->rns_a, tmp[2]);
	// gmp_printf("Q: %Zd\n", t);

	base_conversion(tmp[0], mult->conv, tmp[2]); // Q in base2

	// from_rns_to_int_crt(t, mult->conv->rns_b, tmp[0]);
	// gmp_printf("Q in base2 %Zd\n", t);

	mul_rns(tmp[2], mult->conv->rns_b, tmp[0], mult->p_modMb); // Q*P base2

	// from_rns_to_int_crt(t, mult->conv->rns_b, tmp[2]);
	// gmp_printf("QP in base2 %Zd\n", t);

	add_rns(tmp[0], mult->conv->rns_b, tmp[1], tmp[2]); // A*B + Q*P in base 2

	// from_rns_to_int_crt(t, mult->conv->rns_b, tmp[0]);
	// gmp_printf("AB + QP in base2 %Zd\n", t);

	mul_rns(rop, mult->conv->rns_b, tmp[0], mult->inv_Ma_modMb); // Division by Ma

	// from_rns_to_int_crt(t, mult->conv->rns_b, tmp[0]);
	// gmp_printf("res in base2 %Zd\n", t);
}

void mult_mod_rns_cr(int64_t *rop, int64_t *pa, int64_t *pab, int64_t *pb,
					 int64_t *pbb, struct mod_mul_t *mult, int64_t *tmp[4])
{

	int i;

	mul_rns_cr(tmp[0], mult->conv->rns_a, pa, pb);					  //A*B
	mul_rns_cr(tmp[1], mult->conv->rns_b, pab, pbb);				  //A*B in base2
	mul_rns_cr(tmp[2], mult->conv->rns_a, tmp[0], mult->inv_p_modMa); //Q*{P-1}
	base_conversion_cr(tmp[0], mult->conv, tmp[2], tmp[3]);			  // Q in base2
	mul_rns_cr(tmp[2], mult->conv->rns_b, tmp[0], mult->p_modMb);	  // Q*P base2
	add_rns_cr(tmp[0], mult->conv->rns_b, tmp[1], tmp[2]);			  // A*B + Q*P in base 2
	mul_rns_cr(rop, mult->conv->rns_b, tmp[0], mult->inv_Ma_modMb);	  // Division by Ma
}

void mult_mod_rns_cr_cox(int64_t *rop, int64_t *pa, int64_t *pab, int64_t *pb,
						 int64_t *pbb, struct mod_mul_t *mult, int64_t *tmp[4])
{

	int i;

	mul_rns_cr(tmp[0], mult->conv->rns_a, pa, pb);					  //A*B
	mul_rns_cr(tmp[1], mult->conv->rns_b, pab, pbb);				  //A*B in base2
	mul_rns_cr(tmp[2], mult->conv->rns_a, tmp[0], mult->inv_p_modMa); //Q*{P-1}
	base_conversion_cox(tmp[0], mult->conv, tmp[2], 0, 0, 0);		  //Q in base 2
	mul_rns_cr(tmp[2], mult->conv->rns_b, tmp[0], mult->p_modMb);	  // Q*P base2
	add_rns_cr(tmp[0], mult->conv->rns_b, tmp[1], tmp[2]);			  // A*B + Q*P in base 2
	mul_rns_cr(rop, mult->conv->rns_b, tmp[0], mult->inv_Ma_modMb);	  // Division by Ma
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
	for (i = 0; i < base->size; i++)
	{
		int_inv_Mi[i] = mpz_get_si(inv_Mi[i]);

		//printf("lll %ld \n",int_inv_Mi[i]);

		//	mpz_tdiv_q_2exp(tmp, inv_Mi[i], 32);  //gives the up part
		//	gmp_printf("inv_M[%d] = %Zd\n", i, inv_Mi[i]);
		//	up = mpz_get_si(tmp);
		//	printf("uuu %ld \n",up);

		//	int_inv_Mi[i] += up<<32;
		//	printf("iii %ld \n",int_inv_Mi[i]);
	}
	base->int_inv_Mi = int_inv_Mi;
	mpz_clears(tmp_gcd, tmp_mi, tmp, t, NULL);
}

/////////////////////////////////
// Clear the space reserved for
// the RNS base constants
/////////////////////////////////
// TODO generates error ????
void clear_rns(struct rns_base_t *base)
{
	int i;

	mpz_clear(base->M);
	for (i = 0; i < base->size; i++)
	{
		mpz_clear(base->Mi[i]);
		mpz_clear(base->inv_Mi[i]);
	}
	free(base->Mi);
	free(base->inv_Mi);
}
/////////////////////////////////
// Convert RNS to GMP using
// the CRT
// Assumes "rop" already initialized
/////////////////////////////////
// TODO : NB_COEFF to NB_RES
void from_rns_to_int_crt(mpz_t rop, struct rns_base_t *base, int64_t *op)
{
	mpz_t tmp;
	int i;

	// Initializations
	mpz_init(tmp);
	mpz_set_ui(rop, 0);

	for (i = 0; i < base->size; i++)
	{
		mpz_mul_ui(tmp, base->inv_Mi[i], op[i]); //xi*(Mi**(-1))
		mpz_fdiv_r_ui(tmp, tmp, base->m[i]);	 //xi*(Mi**(-1)) mod mi
		mpz_mul(tmp, tmp, base->Mi[i]);
		mpz_add(rop, rop, tmp);
	}
	//mpz_fdiv_q(tmp, rop, base->M);
	//gmp_printf("k with gmp : %Zd\n", tmp);

	// Modulo M
	mpz_fdiv_r(rop, rop, base->M);

	mpz_clear(tmp);
}

/////////////////////////////////////////////////
// C function for extended Euclidean Algorithm
// gcd = ax + by
/////////////////////////////////////////////////
int64_t gcdExtended(int64_t a, int64_t b, int64_t *x, int64_t *y)
{
	// Base Case
	if (a == 0)
	{
		*x = 0;
		*y = 1;
		return b;
	}

	int64_t x1, y1; // To store results of recursive call
	int64_t gcd = gcdExtended(b % a, a, &x1, &y1);

	// Update x and y using results of recursive
	// call
	*x = y1 - (b / a) * x1;
	*y = x1;

	return gcd;
}

///////////////////////////////////////////////////////
// Initializes the constant inversese for
// the base conversion using MRS
// inv and mrs are supposed to be already initialized
///////////////////////////////////////////////////////
void initialize_inverses_base_conversion(struct conv_base_t *conv_base)
{
	//int64_t **inv, int64_t **mrs, int64_t *m1, int64_t *m2){
	int i, j;
	int64_t tmp, x;
	int128 tmp2;
	//	int64_t res;
	//	int128 prod;
	mpz_t tmpz, tmp_residue, tmp_divisor;

	int size = conv_base->rns_a->size;

	mpz_init(tmpz);
	mpz_init(tmp_residue);
	mpz_init(tmp_divisor);

	// Memory allocation for mrs and inverse
	conv_base->inva_to_b = (int64_t **)malloc(size * sizeof(int64_t *));
	conv_base->mrsa_to_b = (int64_t **)malloc(size * sizeof(int64_t *));
	conv_base->Mi_modPi = (int64_t **)malloc(size * sizeof(int64_t *));
	conv_base->invM_modPi = (int64_t *)malloc(size * sizeof(int64_t));

	for (i = 0; i < NB_COEFF; i++)
	{
		conv_base->inva_to_b[i] = (int64_t *)malloc(size * sizeof(int64_t));
		conv_base->mrsa_to_b[i] = (int64_t *)malloc(size * sizeof(int64_t));
		conv_base->Mi_modPi[i] = (int64_t *)malloc(size * sizeof(int64_t));
	}
	// Initialization of the arrays for mrs and inverse
	for (i = 0; i < size; i++)
	{
		for (j = 0; j < size; j++)
		{
			conv_base->inva_to_b[i][j] = 0;
			conv_base->mrsa_to_b[i][j] = 1;
		}
	}
	// Modular inverses : inv[i][j] = inv (m_i) mod m_j
	for (i = 0; i < size - 1; i++)
	{
		for (j = i + 1; j < size; j++)
		{
			gcdExtended(conv_base->rns_a->m[i], conv_base->rns_a->m[j], &x, &tmp);
			if (x > 0) // Sinon ca part en couille avec 128 et modulo
				conv_base->inva_to_b[i][j] = x;
			else
				conv_base->inva_to_b[i][j] = x + conv_base->rns_a->m[j];
			//			printf(" %ld  %ld\n", inv[i][j], tmp);

			//			prod = ((int128)m1[i]*inv[i][j]);
			//			res = prod%m2[j];
			//			printf("prod %ld  tmp %ld\n", res, tmp);

			//			printf("Inverses m%d = %ld modulo m'%d = %ld: gcd %ld\n",  i, m1[i], j, m2[j], g);
			//			printf(" %ld  %ld\n\n", x, tmp);
		}
	}
	// Modular value of the mrs base
	for (j = 0; j < size; j++)
		conv_base->mrsa_to_b[0][j] = conv_base->rns_a->m[0] % conv_base->rns_b->m[j];
	for (i = 1; i < size - 1; i++)
		for (j = 0; j < size; j++)
		{
			tmp2 = (int128)conv_base->mrsa_to_b[i - 1][j] * conv_base->rns_a->m[i];
			conv_base->mrsa_to_b[i][j] = tmp2 % conv_base->rns_b->m[j];
		}

	//	for(i=0; i<NB_COEFF-1; i++)
	//		for(j = 0; j<NB_COEFF; j++)
	//			printf("mrs[%d, %d] = %ld\n",i ,j, mrs[i][j]);
	// -M in base b
	mpz_neg(tmpz, conv_base->rns_a->M);

	//gmp_printf("M = %Zd \n-M=%Zd\n", conv_base->rns_a->M, tmpz);
	//printf("Taille base b %d\n", conv_base->rns_b->size);

	for (i = 0; i < conv_base->rns_b->size; i++)
	{
		//mpz_set_ui(tmp_divisor, conv_base->rns_b->m[i]);
		mpz_fdiv_r_ui(tmp_residue, tmpz, conv_base->rns_b->m[i]);

		//gmp_printf("p[%d]=%ld MmodPi=%Zd \n", i, conv_base->rns_b->m[i], tmp_residue);

		conv_base->invM_modPi[i] = mpz_get_ui(tmp_residue);
	}
	// Mi mod pj
	for (i = 0; i < conv_base->rns_a->size; i++)
	{
		//gmp_printf("M[%d]=%Zd \n", i, conv_base->rns_a->Mi[i]);

		for (j = 0; j < conv_base->rns_b->size; j++)
		{
			//mpz_set_ui(tmp_divisor, conv_base->rns_b->m[i]);
			mpz_fdiv_r_ui(tmp_residue, conv_base->rns_a->Mi[i], conv_base->rns_b->m[j]);
			conv_base->Mi_modPi[i][j] = mpz_get_ui(tmp_residue);
			//printf("   p[%d]=%ld   Mi_modPi[%d][%d]=%ld\n", j, conv_base->rns_b->m[j], i, j, conv_base->Mi_modPi[i][j]);
		}
	}
	mpz_clears(tmpz, tmp_residue, tmp_divisor, NULL);
}

///////////////////////////////////////////////////////
// Converts a RNS number from a base into an other
// using the MRS conversion
///////////////////////////////////////////////////////
void base_conversion(int64_t *rop, struct conv_base_t *conv_base, int64_t *op)
{
	int i, j;
	int64_t a[NB_COEFF * 10]; ///////////////////////////ARGH !! Le co??t de ce truc n'est pas ma??tris?? !
	int64_t tmp;
	int128 tmp2;
	int size = conv_base->rns_a->size;

	// Set target number to 0
	// for(j=0; j<NB_COEFF; j++)
	// 	rop[j]=0;
	for (j = 0; j < size; j++)
		a[j] = op[j];
	for (i = 0; i < size - 1; i++)
	{
		for (j = i + 1; j < size; j++)
		{
			tmp = a[j] - a[i];
			tmp2 = (int128)tmp * conv_base->inva_to_b[i][j];
			a[j] = tmp2 % conv_base->rns_a->m[j];
			if (a[j] < 0) // Sinon ca part en couille ????????????
				a[j] += conv_base->rns_a->m[j];
		}
	}

	// Residue of the MRS radix
	for (j = 0; j < size; j++)
		rop[j] = a[0] % conv_base->rns_b->m[j];
	for (j = 0; j < size; j++)
	{
		for (i = 1; i < size; i++)
		{
			tmp2 = (int128)a[i] * conv_base->mrsa_to_b[i - 1][j];
			tmp = tmp2 % conv_base->rns_b->m[j];
			rop[j] = ((int128)rop[j] + tmp) % conv_base->rns_b->m[j]; //Overflow risk

			if (rop[j] < 0) //Sinon, ca part en couille  ??????
				rop[j] += conv_base->rns_b->m[j];
		}
	}
}

///////////////////////////////////////////////////////
// Converts a RNS number from a base into an other
// using the MRS conversion. The RNS base uses Crandall
// numbers
///////////////////////////////////////////////////////
void base_conversion_cr(int64_t *rop, struct conv_base_t *conv_base, int64_t *op, int64_t *a)
{
	int i, j;
	//	int64_t a[NB_COEFF];  // En externe, car ??a prend du temps
	int64_t tmp;
	int128 tmp2;
	int128 tmp3;
	int64_t up, up2, lo, lo2;
	int64_t mask = ((int64_t)1 << 63) - 1; /////////////////////////////////////
	int size = conv_base->rns_a->size;

	// Set target number to 0
	// for(j=0; j<NB_COEFF; j++)
	// 	rop[j]=0;
	for (j = 0; j < size; j++) // Indispensable
		a[j] = op[j];
	for (i = 0; i < size - 1; i++)
	{
		for (j = i + 1; j < size; j++)
		{
			tmp = a[j] - a[i];
			a[j] = mul_mod_cr_t(tmp, conv_base->inva_to_b[i][j], conv_base->rns_a->k[j]);
			// if(a[j]<0)	// Sinon ca part en couille ??????? To be Checked
			// 	a[j]+=conv_base->rns_a->m[j];
		}
	}

	// Residue of the MRS radix
	for (j = 0; j < size; j++)
		rop[j] = a[0] % conv_base->rns_b->m[j];
	for (j = 0; j < size; j++)
	{
		for (i = 1; i < size; i++)
		{
			tmp = mul_mod_cr(a[i], conv_base->mrsa_to_b[i - 1][j], conv_base->rns_b->k[j]);
			rop[j] = add_mod_cr(rop[j], tmp, conv_base->rns_b->k[j]);

			// if (rop[j]<0){	//Sinon, ca part en couille ????????????? To be checked
			// 	rop[j] += conv_base->rns_b->m[j];
			// 	printf("conv");
			// }
		}
	}
}

//////////////////////////////////////////////
// Cox-rower method for conversion
//////////////////////////////////////////////
int compute_k_cox(int64_t *op, struct rns_base_t *base, int r, int q, int alpha)
{
	int i;
	int n, sigma, k = 0, k_i;
	int64_t mask, mask2, xhi, trunk; //We certainly could compute with 32bits words, Probably less

	//int32_t low_mask, up_mask;

	r = 63;							 //////////////////////////////////////
	q = 7;							 //////////////////////////////////
	alpha = ((int64_t)1 << (q - 1)); //////////////////////////////////////////

	mask = ((int64_t)1 << r) - ((int64_t)1 << (r - q));
	mask2 = ((int64_t)1 << q);
	n = base->size;

	//printf(" mask = %ld mask2 = %ld  n = %d \n",mask, mask2, n);

	//low_mask = (int32_t)mask;
	//up_mask= (int32_t)((uint64_t)mask>>32);

	//printf("upmask = %d  low mask = %d \n", up_mask, low_mask);

	sigma = alpha;

	//printf("alpha = %d\n",alpha);

	for (i = 0; i < n; i++)
	{
		//		xhi = op[i] * base->int_invM_i[i] % base->m[i];

		xhi = mul_mod_cr(op[i], base->int_inv_Mi[i], base->k[i]); // x_i*invM_i mod m_i
		trunk = xhi & mask;

		//trunk = ((xhi>>32) & up_mask)<<32;  // There are some unusefull operations. Works if q<r/2
		//printf("zzzz xhi = %ld, xhi_shift=%ld  xhi_s_masked=%ld trunk = %ld \n", xhi, xhi>>32, (xhi>>32) & up_mask, trunk);

		sigma += trunk >> (r - q);
		k_i = sigma & mask2;

		//	gmp_printf("inv_Mi[%d] = %Zd ", i, base->inv_Mi[i]);
		//	printf("op[%d]=%ld int_inv_Mi[%d]=%ld xhi[%d]= %ld  trunk[%d]=%ld  sigma[%d]=%d   k[%d]=%d \n ", i, op[i], i, base->int_inv_Mi[i], i, xhi, i, trunk, i, sigma, i, k_i);

		sigma -= k_i;
		k_i = k_i >> q;

		//printf("k[%d]=%d \n", i, k_i);

		k += k_i;
	}
	return k;
}

/////////////////////////////////////////////
// Base conversion using Cox-Rower metod
/////////////////////////////////////////////
void base_conversion_cox(int64_t *rop, struct conv_base_t *conv_base, int64_t *op, int r, int q, int alpha)
{
	int i, j;
	int n, k_i;
	int sigma;
	//unsigned char sigma;
	int64_t mask, mask2, xhi, trunk;   //We certainly could compute with 32bits words, Probably less
	int size = conv_base->rns_a->size; //Should be the size of secondary base
	int64_t tmp, tmp2, tmp3;

	r = 63;							 //////////////////////////////////////
	q = 7;							 //////////////////////////////////
	alpha = ((int64_t)1 << (q - 1)); //////////////////////////////////////////

	mask = ((int64_t)1 << r) - ((int64_t)1 << (r - q));
	mask2 = ((int64_t)1 << q);
	n = conv_base->rns_a->size;

	//printf(" mask = %ld mask2 = %ld  n = %d \n",mask, mask2, n);

	sigma = alpha;

	//printf("alpha = %d\n",alpha);

	// Initialize rop[]
	for (i = 0; i < n; i++)
		rop[i] = 0;
	for (i = 0; i < n; i++)
	{
		//		xhi = op[i] * base->int_invM_i[i] % base->m[i];

		xhi = mul_mod_cr(op[i], conv_base->rns_a->int_inv_Mi[i], conv_base->rns_a->k[i]); // x_i*invM_i mod m_i
		trunk = xhi & mask;

		//printf("zzzz xhi = %ld, xhi_shift=%ld  xhi_s_masked=%ld trunk = %ld \n", xhi, xhi>>32, (xhi>>32) & up_mask, trunk);

		sigma += trunk >> (r - q);
		k_i = sigma & mask2;

		//		gmp_printf("inv_Mi[%d] = %Zd ", i, base->inv_Mi[i]);
		//		printf("??????op[%d]=%ld int_inv_Mi[%d]=%ld xhi[%d]= %ld  trunk[%d]=%ld  sigma[%d]=%d   k[%d]=%d \n ", i, op[i], i, conv_base->rns_a->int_inv_Mi[i], i, xhi, i, trunk, i, sigma, i, k_i);

		sigma -= k_i;
		k_i = k_i >> q; // 0 or 1

		//printf("k[%d]=%d \n", i, k_i);

		for (j = 0; j < size; j++) // Computation of tmp2 has been simplified
		{
			tmp = mul_mod_cr(xhi, conv_base->Mi_modPi[i][j], conv_base->rns_b->k[j]);
			tmp2 = conv_base->invM_modPi[j] * k_i;
			tmp3 = add_mod_cr(tmp, tmp2, conv_base->rns_b->k[j]);
			//rop[j]+=tmp3;
			rop[j] = add_mod_cr(rop[j], tmp3, conv_base->rns_b->k[j]);
		}
	}
}

///////////////////////////////////////////////////
// Modular addition and multiplication using
// Crandall moduli
///////////////////////////////////////////////////
int64_t add_mod_cr(int64_t a, int64_t b, int k)
{
	int64_t tmp, up, lo, mask = ((int64_t)1 << 63) - 1; /////////////////////////////////////
	int64_t res = 0;									// Required, else computation on 32 bits with a sign extension to 64 bits
	uint64_t u_res, u_mod;

	u_mod = ((int64_t)1 << 63) - k; // only 32 shifts without int64_t extension
	tmp = a + b;
	up = (uint64_t)tmp >> 63; //////////////////////////////////////////////////////
							  /// Unsigned conversion to avoid sign extension
	lo = tmp & mask;
	res += lo + up * k;
	u_res = (uint64_t)res;
	//	if ( u_res  >  u_mod	){	// For the unfrequent case :  2^63 > res >=mod ???? To Be Checked
	//		res -=u_mod;
	//		printf("a");
	//	}
	return res;
}

// Crandal modular multiplication with test
int64_t mul_mod_cr_t(int64_t a, int64_t b, int k)
{
	int128 prod;
	uint128 tmp;
	uint128 tmp2;
	int64_t up;
	int64_t lo;
	int64_t up2, up3;
	int64_t lo2, lo3;
	int64_t mask = ((int64_t)1 << 63) - 1; /////////////////////////////////////
	int64_t res = 0;					   // Required, else computation on 32 bits with a sign extension to 64 bits
	uint64_t u_res, u_mod;

	u_mod = ((int64_t)1 << 63) - k; // only 32 shifts without int64_t extension
	prod = (int128)a * b;
	up = (int64_t)((uint128)prod >> 63); //////////////////////////////////
	lo = (int64_t)prod & mask;

	tmp = (int128)up * k;
	lo2 = (int64_t)tmp & mask;
	up2 = (int64_t)((uint128)tmp >> 63); ///////////////////////////////////
	tmp2 = (uint128)lo + lo2 + up2 * k;
	up3 = (uint64_t)(tmp2 >> 63); ///////////////////////////////////
	lo3 = (uint64_t)tmp2 & mask;

	res += lo3 + up3 * k;

	u_res = (uint64_t)res;

	// printf("i");
	if (u_res > u_mod)
	{ // For the unfrequent case :  2^63 > res >=mod ???? TO be Checked
		res -= u_mod;
		printf("m\n");
	}
	return res;
}

int64_t mul_mod_cr(int64_t a, int64_t b, int k)
{
	int128 prod;
	uint128 tmp;
	uint128 tmp2;
	int64_t up;
	int64_t lo;
	int64_t up2, up3;
	int64_t lo2, lo3;
	int64_t mask = ((int64_t)1 << 63) - 1; /////////////////////////////////////
	int64_t res = 0;					   // Required, else computation on 32 bits with a sign extension to 64 bits
	uint64_t u_res, u_mod;

	u_mod = ((int64_t)1 << 63) - k; // only 32 shifts without int64_t extension
	prod = (int128)a * b;
	up = (int64_t)((uint128)prod >> 63); //////////////////////////////////
	lo = (int64_t)prod & mask;

	tmp = (int128)up * k;
	lo2 = (int64_t)tmp & mask;
	up2 = (int64_t)((uint128)tmp >> 63); ///////////////////////////////////
	tmp2 = (uint128)lo + lo2 + up2 * k;
	up3 = (uint64_t)(tmp2 >> 63); ///////////////////////////////////
	lo3 = (uint64_t)tmp2 & mask;

	res += lo3 + up3 * k;

	u_res = (uint64_t)res;

	// printf("i");
	// if ( u_res  >  u_mod ){	// For the unfrequent case :  2^63 > res >=mod ???? TO be Checked
	// 	res -=u_mod;
	// 	printf("m\n");
	// }
	return res;
}
