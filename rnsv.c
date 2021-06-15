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
void from_m256i_to_rns(int64_t *rop, struct rns_base_t *base, __m256i *op){
	int j;
	for (j=0; j<(base->size)/4; j+=1)
	{

		rop[4*j]=_mm256_extract_epi64(op[j],3);
		rop[4*j+1]=_mm256_extract_epi64(op[j],2);
		rop[4*j+2]=_mm256_extract_epi64(op[j],1);
		rop[4*j+3]=_mm256_extract_epi64(op[j],0);
	}

}

///////////////////////////////
// RNS to __m256i convertion
///////////////////////////////
//~ Assumes allocation already done for "rop".
void from_rns_to_m256i(__m256i *rop, struct rns_base_t *base, int64_t *op){
	int j;
	for (j=0; j<(base->size)/4; j+=1)
	{
		rop[j] = _mm256_set_epi64x(op[4*j],op[4*j+1],op[4*j+2],op[4*j+3]);
	}

}


///////////////////////////////
// prints the 8 RNS compounds
///////////////////////////////
void print_RNS(struct rns_base_t *base, int64_t *a){
	int j;
	for (j=0; j<base->size; j++)
	{
		printf("%ld ",a[j]);
	}
	printf("\n");
}

///////////////////////////////
// prints the 8 RNS compounds of the 2 vectors
///////////////////////////////
inline void print_m256i(struct rns_base_t *base, __m256i *a){
	int j;
	for (j=0; j<(base->size)/4; j++)
	{

		printf("%lld %lld %lld %lld\n",_mm256_extract_epi64(a[j],3),
									_mm256_extract_epi64(a[j],2),
									_mm256_extract_epi64(a[j],1),
									_mm256_extract_epi64(a[j],0));
	}
}

///////////////////////////////
// prints the 4 RNS compounds of only one vector
///////////////////////////////
inline void print_alone_m256i(__m256i a){
	printf("->  ");
	printf("%lld %lld %lld %lld\n",_mm256_extract_epi64(a,3),
									_mm256_extract_epi64(a,2),
									_mm256_extract_epi64(a,1),
									_mm256_extract_epi64(a,0));

}


//%%%%%%%%%%%%  ADD  %%%%%%%%%%%%%%%%%%%%%//


///////////////////////////////////////////////////
// Modular addition and multiplication using 
// Crandall moduli
///////////////////////////////////////////////////
__m256i avx_add_mod_cr(__m256i a, __m256i b, __m256i k){

	__m256i ini_res=_mm256_set1_epi64x(0);    // Juste 0. Est-ce vraiment utile ?

	__m256i tmp_mask= _mm256_slli_epi64(_mm256_set1_epi64x(1),63);
	__m256i mask =_mm256_sub_epi64 (tmp_mask,_mm256_set1_epi64x(1));

	__m256i tmp_u_mod = _mm256_slli_epi64(_mm256_set1_epi64x(1),63);  // Pas utilisé
	
	__m256i u_mod =_mm256_sub_epi64 (tmp_u_mod,k);  // pas utilisé
	
	__m256i tmp = _mm256_add_epi64(a,b);

	

	__m256i up = _mm256_srli_epi64 (tmp,63);  // La retenue sortante

	__m256i lo = _mm256_and_si256(tmp,mask);  // La partie basse de la somme

	__m256i tmp_res =  _mm256_madd_epi16 (up,k);   // mul et add ? Il ne manque pas un terme ?

	__m256i tmp_res2 = _mm256_add_epi64(lo,tmp_res);

	__m256i res = _mm256_add_epi64(ini_res,tmp_res2);  // 0 + tmp_res2 ?????

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
inline void avx_add_rns_cr(__m256i *rop, struct rns_base_t *base, __m256i *pa, __m256i *pb){
	int j;
	
	for (j=0; j<(base->size)/4; j+=1)
	{
		rop[j] = avx_add_mod_cr(pa[j], pb[j], base->avx_k[j]);
		//attention, les coeffs sont à l'envers !

	}
}


//%%%%%%%%%%%%  SUB  %%%%%%%%%%%%%%%%%%%%%//

///////////////////////////////
// RNS substraction
///////////////////////////////
// rop : result
// base : RNS base 
// pa : A
// pb : B
inline void avx_sub_rns_cr(__m256i *rop, struct rns_base_t *base, __m256i *pa, __m256i *pb){
	int j;

	for (j=0; j<(base->size)/4; j+=1)
	{

		__m256i tmp = _mm256_sub_epi64(pa[j],pb[j]);

		__m256i up = _mm256_srli_epi64 (tmp,63);

		__m256i tmp_rop =  _mm256_madd_epi16 (up,base->avx_k[j]);

		rop[j] = _mm256_add_epi64(tmp,tmp_rop);

		

	}
}


//%%%%%%%%%%%%  MUL  %%%%%%%%%%%%%%%%%%%%%//

///////////////////////////////
// auxiliar addition between two terms
///////////////////////////////
// result = 2^63 * rop_up + rop_lo
void avx_add_aux_2e(__m256i *rop_up, __m256i *rop_lo, __m256i a, __m256i b){
	//given two positive numbers on 63 bits

	__m256i tmp_mask3= _mm256_slli_epi64(_mm256_set1_epi64x(1),62);
	__m256i tmp_mask2 =_mm256_sub_epi64 (tmp_mask3,_mm256_set1_epi64x(1));
	__m256i mask1 = _mm256_add_epi64(tmp_mask2,
				_mm256_slli_epi64(_mm256_set1_epi64x(1),62));
	__m256i mask2 = _mm256_slli_epi64(_mm256_set1_epi64x(1),63);

	__m256i sum = _mm256_add_epi64(a,b);
	*rop_lo = _mm256_and_si256(sum,mask1);
	*rop_up = _mm256_srli_epi64(_mm256_and_si256(sum,mask2),63);

}

///////////////////////////////
// auxiliar addition between three terms
///////////////////////////////
// result = 2^63 * rop_up + rop_lo
void avx_add_aux_3e(__m256i *rop_up, __m256i *rop_lo, __m256i a, __m256i b, __m256i c){

	__m256i up,lo,up2,lo2;
	avx_add_aux_2e(&up,&lo,a,b);
	avx_add_aux_2e(&up2,rop_lo,lo,c);
	*rop_up = _mm256_add_epi64(up,up2);

	}

///////////////////////////////
// auxiliar multiplication between 2 63-bits numbers
///////////////////////////////
// result = 2^63 * rop_up + rop_lo
// had to divide the numbers into asymetrical parts, hence the tmp's
// having 62, 63 or 64 bits
void avx_mul_aux(__m256i *rop_up, __m256i *rop_lo, __m256i a, __m256i b){

	__m256i tmp_mask31 = _mm256_slli_epi64(_mm256_set1_epi64x(1),31);
	__m256i tmp_mask32 = _mm256_slli_epi64(_mm256_set1_epi64x(1),32);
	__m256i mask_lo31 = _mm256_sub_epi64 (tmp_mask31,_mm256_set1_epi64x(1));
	__m256i mask_lo32 = _mm256_sub_epi64 (tmp_mask32,_mm256_set1_epi64x(1));
	__m256i mask_up31 = _mm256_slli_epi64(mask_lo31,32);
	__m256i mask_up32 = _mm256_slli_epi64(mask_lo32,31);

	__m256i tmp_mask33 = _mm256_slli_epi64(_mm256_set1_epi64x(1),33);
	__m256i mask_lo33 = _mm256_sub_epi64 (tmp_mask33,_mm256_set1_epi64x(1));

	__m256i a_lo = _mm256_and_si256(a,mask_lo31); //31 bits
	__m256i a_up = _mm256_srli_epi64(_mm256_and_si256(a,mask_up32),31); //32 bits
	__m256i b_lo = _mm256_and_si256(b,mask_lo32); // 32 bits
	__m256i b_up = _mm256_srli_epi64(_mm256_and_si256(b,mask_up31),32); //31 bits

	//we could rework the way we apply the masks, so it's optimised
	//for example, we should initialize them in extern

	__m256i tmp1 = _mm256_mul_epu32(a_lo,b_lo); //63 bits 
	__m256i tmp2 = _mm256_mul_epu32(a_lo,b_up); //62 bits
	__m256i tmp3 = _mm256_mul_epu32(a_up,b_lo); //64 bits
	__m256i tmp4 = _mm256_mul_epu32(a_up,b_up); // 63 bits

	__m256i ret;
	avx_add_aux_3e(&ret,rop_lo,tmp1,
		_mm256_srli_epi64(_mm256_slli_epi64(tmp2,33),1),
		_mm256_srli_epi64(_mm256_slli_epi64(tmp3,32),1));

	*rop_up = _mm256_add_epi64(_mm256_srli_epi64(tmp2,31),
		_mm256_add_epi64(_mm256_srli_epi64(tmp3,32), 
			_mm256_add_epi64(ret,tmp4)));

	}


__m256i avx_mul_mod_cr(__m256i a, __m256i b, __m256i k){

	__m256i tmp_mask= _mm256_slli_epi64(_mm256_set1_epi64x(1),63);
	__m256i mask =_mm256_sub_epi64 (tmp_mask,_mm256_set1_epi64x(1));

	__m256i tmp_u_mod = _mm256_slli_epi64(_mm256_set1_epi64x(1),63);	
	__m256i u_mod =_mm256_sub_epi64 (tmp_u_mod,k);

	__m256i up, lo, up2, lo2, up3, lo3;

	avx_mul_aux(&up,&lo,a,b);

	avx_mul_aux(&up2,&lo2,up,k);

	__m256i up2_times_k, ret1, ret2;
	avx_mul_aux(&ret1,&up2_times_k,up2,k);
	avx_add_aux_3e(&ret2,&lo3,lo,lo2,up2_times_k);
	up3 = _mm256_add_epi64(ret1,ret2);

	__m256i res = _mm256_add_epi64(lo3,_mm256_madd_epi16(up3,k));

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
inline void avx_mul_rns_cr(__m256i *rop, struct rns_base_t *base, __m256i *pa, __m256i *pb){
	int j;
	
	for (j=0; j<(base->size)/4; j+=1)
	{
		rop[j] = avx_mul_mod_cr(pa[j], pb[j], base->avx_k[j]); //attention, les coeffs sont à l'envers !

	}
}

//%%%%%%%%%%%%  MOD MUL  %%%%%%%%%%%%%%%%%%%%%//



void avx_init_mrs(struct conv_base_t *conv_base){
	int i;
	int size = conv_base->rns_a->size;
	conv_base->avx_mrsa_to_b = (__m256i**) malloc(size*sizeof(__m256i*)/4);
	//don't know why but sizeof has to be divided by 8 so it work (maybe a bit/byte problem)

	for(i=0; i<size; i++){
		conv_base->avx_mrsa_to_b[i] = (__m256i*) malloc(size*sizeof(__m256i)/16);
		
	}
	for(i=0; i<size; i++){
		from_rns_to_m256i(conv_base->avx_mrsa_to_b[i],
			conv_base->rns_a, conv_base->mrsa_to_b[i]);
	}
}//call before calling function below


///////////////////////////////////////////////////////
// Converts a RNS number from a base into an other
// using the MRS conversion. The RNS base uses Crandall 
// numbers
///////////////////////////////////////////////////////
void avx_base_conversion_cr(__m256i *rop, struct conv_base_t *conv_base, __m256i *op, int64_t *a){
	int i, j;
//	int64_t a[NB_COEFF];  // En externe, car ça prend du temps 
	int64_t tmp;
	__m256i avx_tmp;
	int128 tmp2;
	int128 tmp3;
	int64_t up, up2, lo, lo2;

	int size = conv_base->rns_a->size;

	from_m256i_to_rns(a,conv_base->rns_a,op);

	for(i=0; i<size-1; i++) 
	{
		for(j = i+1; j<size; j++)
		{
			tmp = a[j]-a[i];
			a[j] = mul_mod_cr(tmp, conv_base->inva_to_b[i][j],  conv_base->rns_a->k[j]);
			if(a[j]<0)	// Sinon ca part en couille ?
				a[j]+=conv_base->rns_a->m[j];
		}

	}
	
	// Residue of the MRS radix

	int64_t tmp_rop[NB_COEFF];

	for(j=0; j<size; j++)
	{
		tmp_rop[j]=a[0] % conv_base->rns_b->m[j];	
	}

	from_rns_to_m256i(rop,conv_base->rns_b,tmp_rop);

	for(j=0; j<size/4; j++)
	{
		for(i=1; i<size ; i++)
		{
			avx_tmp = (_mm256_set1_epi64x(a[i]),
				conv_base->avx_mrsa_to_b[i-1][j], conv_base->rns_b->avx_k[j]);
			rop[j]= avx_add_mod_cr(rop[j], avx_tmp, conv_base->rns_b->avx_k[j]);

			if (&rop[j]<0){	//Sinon, ca part en couille ????????????? jamais atteint en pratique
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
void avx_mult_mod_rns_cr(__m256i *rop, __m256i *pa, __m256i *pab, __m256i *pb, 
	__m256i *pbb, struct mod_mul_t *mult, __m256i *tmp0, __m256i *tmp1, __m256i *tmp2, int64_t *a){


	avx_mul_rns_cr(tmp0, mult->conv->rns_a, pa, pb); //A*B
	avx_mul_rns_cr(tmp1, mult->conv->rns_b, pab, pbb); //A*B in base2
	
	avx_mul_rns_cr(tmp2, mult->conv->rns_a, tmp0, mult->avx_inv_p_modMa); //Q*{P-1}

	avx_base_conversion_cr(tmp0, mult->conv, tmp2, a); // Q in base2
	avx_mul_rns_cr(tmp2 , mult->conv->rns_b, tmp0, mult->avx_p_modMb); // Q*P base2
	avx_add_rns_cr(tmp0 , mult->conv->rns_b, tmp1, tmp2); // A*B + Q*P in base 2
	avx_mul_rns_cr(rop, mult->conv->rns_b, tmp0, mult->avx_inv_Ma_modMb); // Division by Ma

}


// COX

int avx_compute_k_cox(__m256i *op, struct rns_base_t *base, int r, int q, int alpha) {
	return 0;
}

void avx_base_conversion_cox(__m256i *rop, struct conv_base_t *conv_base, __m256i *op, int r, int q, int alpha) {

}

void avx_mult_mod_rns_cr_cox(__m256i *rop, __m256i *pa, __m256i *pab, __m256i *pb, 
	__m256i *pbb, struct mod_mul_t *mult, __m256i *tmp[4]) {

	}