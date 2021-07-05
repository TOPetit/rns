#include "rnsv.h"

#include <stdlib.h>

#include <stdio.h>

#include <unistd.h>

#include <stdint.h>

#include <string.h>

#include <time.h>

#include <gmp.h>

#include <immintrin.h>

#include <time.h>

///////////////////////////////
// RNS to __m256i convertion
///////////////////////////////
//~ Assumes allocation already done for "rop".
inline void from_m256i_to_rns_bis(int64_t *rop, struct rns_base_t *base, __m256i *op)
{
	int j;
	for (j = 0; j < (base->size) / 4; j += 1)
	{

		rop[4 * j] = _mm256_extract_epi64(op[j], 0);
		rop[4 * j + 1] = _mm256_extract_epi64(op[j], 1);
		rop[4 * j + 2] = _mm256_extract_epi64(op[j], 2);
		rop[4 * j + 3] = _mm256_extract_epi64(op[j], 3);
	}
}

///////////////////////////////
// RNS to __m256i convertion without extract (more efficient ~100 cycles less)
// Better to code it directly in function and removing for when timing is important
///////////////////////////////
//~ Assumes allocation already done for "rop".
inline void from_m256i_to_rns(int64_t *rop, struct rns_base_t *base, __m256i *op)
{
	for (int i = 0; i < NB_COEFF / 4; i++)
	{
		_mm256_storeu_si256((__m256i *)&rop[4 * i], op[i]);
	}
}

///////////////////////////////
// RNS to __m256i convertion
///////////////////////////////
//~ Assumes allocation already done for "rop".
inline void from_rns_to_m256i(__m256i *rop, struct rns_base_t *base, int64_t *op)
{
	int j;
	for (j = 0; j < (base->size) / 4; j += 1)
	{
		rop[j] = _mm256_set_epi64x(op[4 * j + 3], op[4 * j + 2], op[4 * j + 1], op[4 * j]);
	}
}

///////////////////////////////
// prints the 8 RNS compounds
///////////////////////////////
inline void print_RNS(struct rns_base_t *base, int64_t *a)
{
	int j;
	for (j = 0; j < base->size; j++)
	{
		printf("%ld ", a[j]);
	}
	printf("\n");
}

///////////////////////////////
// prints the 8 RNS compounds of the 2 vectors
///////////////////////////////
inline void print_m256i(struct rns_base_t *base, __m256i *a)
{
	int j;
	for (j = 0; j < (base->size) / 4; j++)
	{

		printf("%lld %lld %lld %lld ", _mm256_extract_epi64(a[j], 3),
			   _mm256_extract_epi64(a[j], 2),
			   _mm256_extract_epi64(a[j], 1),
			   _mm256_extract_epi64(a[j], 0));
	}
}

///////////////////////////////
// prints the 4 RNS compounds of only one vector
///////////////////////////////
inline void print_alone_m256i(__m256i a)
{
	printf("->  ");
	printf("%lld %lld %lld %lld\n", _mm256_extract_epi64(a, 3),
		   _mm256_extract_epi64(a, 2),
		   _mm256_extract_epi64(a, 1),
		   _mm256_extract_epi64(a, 0));
}

inline void avx_init_rns(struct rns_base_t *base)
{

	int n = base->size;

	__m256i *avx_inv_Mi = (__m256i *)malloc(n * sizeof(__m256i));

	for (int i = 0; i < n; i++)
	{
		avx_inv_Mi[i] = _mm256_set1_epi64x(base->int_inv_Mi[i]);
	}

	base->avx_inv_Mi = avx_inv_Mi;

	__m256i *tmp = (__m256i *)malloc(n * sizeof(__m256i) / 4);

	for (int i = 0; i < n / 4; i++)
	{
		tmp[i] = _mm256_set_epi64x(base->m[4 * i + 3], base->m[4 * i + 2], base->m[4 * i + 1], base->m[4 * i]);
	}

	base->avx_m = tmp;
}

//%%%%%%%%%%%%  ADD  %%%%%%%%%%%%%%%%%%%%%//

///////////////////////////////////////////////////
// Modular addition and multiplication using
// Crandall moduli
///////////////////////////////////////////////////
inline __m256i avx_add_mod_cr(__m256i a, __m256i b, __m256i k)
{

	__m256i tmp_mask = _mm256_slli_epi64(_mm256_set1_epi64x(1), 63);
	__m256i mask = _mm256_sub_epi64(tmp_mask, _mm256_set1_epi64x(1));

	__m256i tmp = _mm256_add_epi64(a, b);

	__m256i up = _mm256_srli_epi64(tmp, 63); // La retenue sortante

	__m256i lo = _mm256_and_si256(tmp, mask); // La partie basse de la somme

	__m256i tmp_res = _mm256_madd_epi16(up, k); // mul et add ? Il ne manque pas un terme ?

	__m256i res = _mm256_add_epi64(lo, tmp_res);

	//on verra après pour le if

	return res;
}

///////////////////////////////
// RNS addition
///////////////////////////////
// rop : result
// base : RNS base
// pa : A
// pb : B
inline void avx_add_rns_cr(__m256i *rop, struct rns_base_t *base, __m256i *pa, __m256i *pb)
{
	int j;

	for (j = 0; j < (base->size) / 4; j += 1)
	{
		rop[j] = avx_add_mod_cr(pa[j], pb[j], base->avx_k[j]);
		//attention, les coeffs sont à l'envers !
	}
}

//%%%%%%%%%%%%  SUB  %%%%%%%%%%%%%%%%%%%%%//

///////////////////////////////////////////////////
// Modular substraction and multiplication using
// Crandall moduli
///////////////////////////////////////////////////
__m256i avx_sub_mod_cr(__m256i a, __m256i b, __m256i k, __m256i m)
{

	__m256i tmp_mask = _mm256_slli_epi64(_mm256_set1_epi64x(1), 63);
	__m256i mask = _mm256_sub_epi64(tmp_mask, _mm256_set1_epi64x(1));

	__m256i tmp1 = _mm256_add_epi64(a, m);

	__m256i tmp = _mm256_sub_epi64(tmp1, b);

	__m256i up = _mm256_srli_epi64(tmp, 63); // La retenue sortante

	__m256i lo = _mm256_and_si256(tmp, mask); // La partie basse de la somme

	__m256i tmp_res = _mm256_madd_epi16(up, k); // mul et add ? Il ne manque pas un terme ?

	__m256i res = _mm256_add_epi64(lo, tmp_res);

	return res;
}

///////////////////////////////
// RNS substraction
///////////////////////////////
// rop : result
// base : RNS base
// pa : A
// pb : B
inline void avx_sub_rns_cr(__m256i *rop, struct rns_base_t *base, __m256i *pa, __m256i *pb)
{
	int j;

	for (j = 0; j < (base->size) / 4; j += 1)
	{

		rop[j] = avx_sub_mod_cr(pa[j], pb[j], base->avx_k[j], base->avx_m[j]);
	}
}

//%%%%%%%%%%%%  MUL  %%%%%%%%%%%%%%%%%%%%%//

///////////////////////////////
// auxiliar addition between two terms
///////////////////////////////
// result = 2^63 * rop_up + rop_lo
inline void avx_add_aux_2e(__m256i *rop_up, __m256i *rop_lo, __m256i a, __m256i b)
{
	//given two positive numbers on 63 bits

	__m256i tmp_mask3 = _mm256_slli_epi64(_mm256_set1_epi64x(1), 62);
	__m256i tmp_mask2 = _mm256_sub_epi64(tmp_mask3, _mm256_set1_epi64x(1));
	__m256i mask1 = _mm256_add_epi64(tmp_mask2,
									 _mm256_slli_epi64(_mm256_set1_epi64x(1), 62));
	__m256i mask2 = _mm256_slli_epi64(_mm256_set1_epi64x(1), 63);

	__m256i sum = _mm256_add_epi64(a, b);
	*rop_lo = _mm256_and_si256(sum, mask1);
	*rop_up = _mm256_srli_epi64(_mm256_and_si256(sum, mask2), 63);
}

///////////////////////////////
// auxiliar addition between three terms
///////////////////////////////
// result = 2^63 * rop_up + rop_lo
inline void avx_add_aux_3e(__m256i *rop_up, __m256i *rop_lo, __m256i a, __m256i b, __m256i c)
{

	__m256i up, lo, up2, lo2;
	avx_add_aux_2e(&up, &lo, a, b);
	avx_add_aux_2e(&up2, rop_lo, lo, c);
	*rop_up = _mm256_add_epi64(up, up2);
}

///////////////////////////////
// auxiliar multiplication between 2 63-bits numbers
///////////////////////////////
// result = 2^63 * rop_up + rop_lo
// had to divide the numbers into asymetrical parts, hence the tmp's
// having 62, 63 or 64 bits
inline void avx_mul_aux(__m256i *rop_up, __m256i *rop_lo, __m256i a, __m256i b)
{

	__m256i tmp_mask31 = _mm256_slli_epi64(_mm256_set1_epi64x(1), 31);
	__m256i tmp_mask32 = _mm256_slli_epi64(_mm256_set1_epi64x(1), 32);
	__m256i mask_lo31 = _mm256_sub_epi64(tmp_mask31, _mm256_set1_epi64x(1));
	__m256i mask_lo32 = _mm256_sub_epi64(tmp_mask32, _mm256_set1_epi64x(1));
	__m256i mask_up31 = _mm256_slli_epi64(mask_lo31, 32);
	__m256i mask_up32 = _mm256_slli_epi64(mask_lo32, 31);

	__m256i tmp_mask33 = _mm256_slli_epi64(_mm256_set1_epi64x(1), 33);
	__m256i mask_lo33 = _mm256_sub_epi64(tmp_mask33, _mm256_set1_epi64x(1));

	__m256i a_lo = _mm256_and_si256(a, mask_lo31);						  //31 bits
	__m256i a_up = _mm256_srli_epi64(_mm256_and_si256(a, mask_up32), 31); //32 bits
	__m256i b_lo = _mm256_and_si256(b, mask_lo32);						  // 32 bits
	__m256i b_up = _mm256_srli_epi64(_mm256_and_si256(b, mask_up31), 32); //31 bits

	//we could rework the way we apply the masks, so it's optimised
	//for example, we should initialize them in extern

	__m256i tmp1 = _mm256_mul_epu32(a_lo, b_lo); //63 bits
	__m256i tmp2 = _mm256_mul_epu32(a_lo, b_up); //62 bits
	__m256i tmp3 = _mm256_mul_epu32(a_up, b_lo); //64 bits
	__m256i tmp4 = _mm256_mul_epu32(a_up, b_up); // 63 bits

	__m256i ret;
	avx_add_aux_3e(&ret, rop_lo, tmp1,
				   _mm256_srli_epi64(_mm256_slli_epi64(tmp2, 33), 1),
				   _mm256_srli_epi64(_mm256_slli_epi64(tmp3, 32), 1));

	*rop_up = _mm256_add_epi64(_mm256_srli_epi64(tmp2, 31),
							   _mm256_add_epi64(_mm256_srli_epi64(tmp3, 32),
												_mm256_add_epi64(ret, tmp4)));
}

inline __m256i avx_mul_mod_cr(__m256i a, __m256i b, __m256i k)
{

	__m256i tmp_mask = _mm256_slli_epi64(_mm256_set1_epi64x(1), 63);
	__m256i mask = _mm256_sub_epi64(tmp_mask, _mm256_set1_epi64x(1));

	__m256i tmp_u_mod = _mm256_slli_epi64(_mm256_set1_epi64x(1), 63);
	__m256i u_mod = _mm256_sub_epi64(tmp_u_mod, k);

	__m256i up, lo, up2, lo2, up3, lo3;

	avx_mul_aux(&up, &lo, a, b);

	avx_mul_aux(&up2, &lo2, up, k);

	__m256i up2_times_k, ret1, ret2;
	avx_mul_aux(&ret1, &up2_times_k, up2, k);
	avx_add_aux_3e(&ret2, &lo3, lo, lo2, up2_times_k);
	up3 = _mm256_add_epi64(ret1, ret2);

	__m256i res = _mm256_add_epi64(lo3, _mm256_madd_epi16(up3, k));

	//on verra après pour le if

	return res;
}

///////////////////////////////
// RNS multiplication
///////////////////////////////
// rop : result
// base : RNS base
// pa : A
// pb : B
inline void avx_mul_rns_cr(__m256i *rop, struct rns_base_t *base, __m256i *pa, __m256i *pb)
{
	int j;

	for (j = 0; j < (base->size) / 4; j += 1)
	{
		rop[j] = avx_mul_mod_cr(pa[j], pb[j], base->avx_k[j]); //attention, les coeffs sont à l'envers !
	}
}

//%%%%%%%%%%%%  MOD MUL  %%%%%%%%%%%%%%%%%%%%%//

inline void avx_init_mrs(struct conv_base_t *conv_base)
{
	int i;
	int size = conv_base->rns_a->size;
	conv_base->avx_mrsa_to_b = (__m256i **)malloc(size * sizeof(__m256i *));
	//don't know why but sizeof has to be divided by 8 so it work (maybe a bit/byte problem)

	for (i = 0; i < size; i++)
	{
		conv_base->avx_mrsa_to_b[i] = (__m256i *)malloc(size * sizeof(__m256i) / 4);
	}
	for (i = 0; i < size; i++)
	{
		from_rns_to_m256i(conv_base->avx_mrsa_to_b[i], conv_base->rns_a, conv_base->mrsa_to_b[i]);
		//print_m256i(conv_base->rns_a, conv_base->avx_mrsa_to_b[i]);
		//printf("\n");
	}
} //call before calling function below

inline void avx_initialize_inverses_base_conversion(struct conv_base_t *conv_base)
{

	int size = conv_base->rns_a->size;

	__m256i **tmp_Arr;

	tmp_Arr = (__m256i **)_mm_malloc(size * sizeof(__m256i *), 32);

	for (int i = 0; i < size; i++)
	{
		tmp_Arr[i] = (__m256i *)_mm_malloc(size * sizeof(__m256i), 32);
	}

	for (int i = 0; i < size; i++)
	{
		for (int j = 0; j < size; j++)
		{
			tmp_Arr[i][j] = _mm256_set1_epi64x(1) * conv_base->Mi_modPi[i][j];
			//printf("%d %d %ld\n",i,j, conv_base->Mi_modPi[i][j]);
		}
	}
	conv_base->avx_Mi_modPi = tmp_Arr;
}

///////////////////////////////////////////////////////
// Converts a RNS number from a base into an other
// using the MRS conversion. The RNS base uses Crandall
// numbers
///////////////////////////////////////////////////////
inline void avx_base_conversion_cr(__m256i *rop, struct conv_base_t *conv_base, __m256i *op, int64_t *a)
{
	int i, j;
	//	int64_t a[NB_COEFF];  // En externe, car ça prend du temps
	int64_t tmp;
	__m256i avx_tmp;
	int128 tmp2;
	int128 tmp3;
	int64_t up, up2, lo, lo2;

	int size = conv_base->rns_a->size;

	from_m256i_to_rns(a, conv_base->rns_a, op);

	for (i = 0; i < size - 1; i++)
	{
		for (j = i + 1; j < size; j++)
		{
			tmp = a[j] - a[i];
			a[j] = mul_mod_cr(tmp, conv_base->inva_to_b[i][j], conv_base->rns_a->k[j]);
			if (a[j] < 0) // Sinon ca part en couille ?
				a[j] += conv_base->rns_a->m[j];
		}
	}

	// Residue of the MRS radix

	int64_t tmp_rop[NB_COEFF];

	for (j = 0; j < size; j++)
	{
		tmp_rop[j] = a[0] % conv_base->rns_b->m[j];
	}

	from_rns_to_m256i(rop, conv_base->rns_b, tmp_rop);

	for (j = 0; j < size / 4; j++)
	{
		for (i = 1; i < size; i++)
		{
			int64_t teee = a[i];
			//avx_tmp = avx_mul_mod_cr(_mm256_set1_epi64x(a[i]), conv_base->avx_mrsa_to_b[i - 1][j], conv_base->rns_b->avx_k[j]);
			//rop[j] = avx_add_mod_cr(rop[j], avx_tmp, conv_base->rns_b->avx_k[j]);

			if (&rop[j] < 0)
			{ //Sinon, ca part en couille ????????????? jamais atteint en pratique
				rop[j] += conv_base->rns_b->m[j];
				printf("conv");
			}
			//on verra pour le if après
		}
	}
}

//TODO : function to initialize avx_p_modMa and the others

///////////////////////////////
// RNS Modular multiplication
///////////////////////////////
// rop : result
// pa : A
// pb : B
// mult : constants
// tmp and a: temporary arrays for intermediate results
inline void avx_mult_mod_rns_cr(__m256i *rop, __m256i *pa, __m256i *pab, __m256i *pb,
								__m256i *pbb, struct mod_mul_t *mult, __m256i *tmp0, __m256i *tmp1, __m256i *tmp2, int64_t *a)
{

	avx_mul_rns_cr(tmp0, mult->conv->rns_a, pa, pb);   //A*B
	avx_mul_rns_cr(tmp1, mult->conv->rns_b, pab, pbb); //A*B in base2

	avx_mul_rns_cr(tmp2, mult->conv->rns_a, tmp0, mult->avx_inv_p_modMa); //Q*{P-1}

	avx_base_conversion_cr(tmp0, mult->conv, tmp2, a);					  // Q in base2
	avx_mul_rns_cr(tmp2, mult->conv->rns_b, tmp0, mult->avx_p_modMb);	  // Q*P base2
	avx_add_rns_cr(tmp0, mult->conv->rns_b, tmp1, tmp2);				  // A*B + Q*P in base 2
	avx_mul_rns_cr(rop, mult->conv->rns_b, tmp0, mult->avx_inv_Ma_modMb); // Division by Ma
}

// COX

inline void avx_base_conversion_cox(__m256i *rop, struct conv_base_t *conv_base, __m256i *op, int r, int q, int alpha)
{
	__m256i unit = _mm256_set1_epi64x(1);
	int i, j;
	int n;
	__m256i sigma, k_i;
	//unsigned char sigma;
	__m256i mask, mask2, xhi, trunk;   //We certainly could compute with 32bits words, Probably less
	int size = conv_base->rns_b->size; //Should be the size of secondary base
	__m256i tmp, tmp2, tmp3;

	r = 63; //////////////////////////////////////
	q = 7;	//////////////////////////////////

	mask = _mm256_slli_epi64(_mm256_set1_epi64x(1), r) - _mm256_slli_epi64(_mm256_set1_epi64x(1), r - q);
	mask2 = _mm256_slli_epi64(_mm256_set1_epi64x(1), q);
	n = conv_base->rns_a->size;

	//printf(" mask = %ld mask2 = %ld  n = %d \n",mask, mask2, n);

	sigma = _mm256_slli_epi64(_mm256_set1_epi64x(1), q - 1);

	//printf("alpha = %d\n",alpha);

	// Initialize rop[]

	for (i = 0; i < n; i++)
	{
		//		xhi = op[i] * base->int_invM_i[i] % base->m[i];

		xhi = avx_mul_mod_cr(op[i], conv_base->rns_a->avx_inv_Mi[i], conv_base->rns_a->avx_k[i]); // x_i*invM_i mod m_i
		trunk = xhi & mask;

		//printf("zzzz xhi = %ld, xhi_shift=%ld  xhi_s_masked=%ld trunk = %ld \n", xhi, xhi>>32, (xhi>>32) & up_mask, trunk);

		sigma += trunk >> (r - q);
		k_i = sigma & mask2;

		//		gmp_printf("inv_Mi[%d] = %Zd ", i, base->inv_Mi[i]);
		//		printf("£££op[%d]=%ld int_inv_Mi[%d]=%ld xhi[%d]= %ld  trunk[%d]=%ld  sigma[%d]=%d   k[%d]=%d \n ", i, op[i], i, conv_base->rns_a->int_inv_Mi[i], i, xhi, i, trunk, i, sigma, i, k_i);

		sigma -= k_i;
		k_i = k_i >> q; // 0 or 1

		//printf("k[%d]=%d \n", i, k_i);

		for (j = 0; j < size / 4; j++) // Computation of tmp2 has been simplified
		{
			// TO SHOW IT'S NOT WORTH IT
			//tmp = avx_mul_mod_cr(xhi, unit, conv_base->rns_b->avx_k[j]);
			tmp = avx_mul_mod_cr(xhi, conv_base->avx_Mi_modPi[i][j], conv_base->rns_b->avx_k[j]);
			tmp2 = conv_base->invM_modPi[j] * k_i;
			tmp3 = avx_add_mod_cr(tmp, tmp2, conv_base->rns_b->avx_k[j]);
			//rop[j]+=tmp3;
			rop[j] = avx_add_mod_cr(rop[j], tmp3, conv_base->rns_b->avx_k[j]);
		}
	}
}

inline void avx_mult_mod_rns_cr_cox(__m256i *rop, __m256i *pa, __m256i *pab, __m256i *pb,
									__m256i *pbb, struct mod_mul_t *mult, __m256i *tmp0, __m256i *tmp1, __m256i *tmp2, int64_t *a)
{

	int i;

	avx_mul_rns_cr(tmp0, mult->conv->rns_a, pa, pb);					  //A*B
	avx_mul_rns_cr(tmp1, mult->conv->rns_b, pab, pbb);					  //A*B in base2
	avx_mul_rns_cr(tmp2, mult->conv->rns_a, tmp0, mult->avx_inv_p_modMa); //Q*{P-1}
	avx_base_conversion_cr(tmp0, mult->conv, tmp2, a);					  //Q in base 2
	avx_mul_rns_cr(tmp2, mult->conv->rns_b, tmp0, mult->avx_p_modMb);	  // Q*P base2
	avx_add_rns_cr(tmp0, mult->conv->rns_b, tmp1, tmp2);				  // A*B + Q*P in base 2
	avx_mul_rns_cr(rop, mult->conv->rns_b, tmp0, mult->avx_inv_Ma_modMb); // Division by Ma
}