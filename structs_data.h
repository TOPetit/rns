#include <gmp.h>
#include <stdint.h>
#include <immintrin.h>

#ifndef STRUCTS_DATA

#define STRUCTS_DATA 


#define K_VAL 32	 	
#define RED_KE 62  		 
#define POLY_DEG 7		 // the degree of the polynomials
#define NB_COEFF 8		 // number of coeffs for every polynomial
#define COEFF_SIZE 33	 // coeff size of polynomials in the amns
//~ #define MAT_W 3  		 // matrix weight
//~ #define AMNS_C 2		 // E(X) = X^NB_COEFF - AMNS_C


typedef __int128 int128;
typedef unsigned __int128 uint128;

struct rns_base_t
{
	unsigned int size;  // RNS base size
	int64_t *m;         // moduli
	int *k;				// m=2^XX -k for crandal base
	__m256i *avx_k;     // k vectored
	int64_t *p;         // modulus p expressed in the RNS base
	int64_t *inv_p;     // p^{-1} mod M in the RNS base
	mpz_t *Mi;          // Mi for the CRT conversion
	mpz_t *inv_Mi;      // Mi^{-1} mod mi for the CRT conversion
	int64_t *int_inv_Mi; // Mi^{-1} mod mi for the Cox conversion
	mpz_t M;            // M for the CRT conversion
	__m256i *avx_inv_Mi; // inv_Mi vectored
};

struct conv_base_t //Constants for the RNSa -> RNSb conversion
{
	struct rns_base_t *rns_a;  // RNSa
	struct rns_base_t *rns_b;  // RNSb
	int64_t **inva_to_b;       // modular inverses of RNSa modulo RNSb
	int64_t **mrsa_to_b;       // MRS Radices modulo RNSb
	__m256i **avx_mrsa_to_b; // mrsa_to_b vectored

	int64_t *invM_modPi;       // -M^{-1}mod pi for Cox conversion /////////////////////
	int64_t **Mi_modPi;        // Mi mod pj  for Cox conversion    /////////////////////
	__m256i *avx_invM_modPi;
	__m256i **avx_Mi_modPi;
};

struct mod_mul_t // Constants for modular multiplication from RNSa to RNSb
{
	int64_t *inv_p_modMa ; //(-P)^-1 mod M
	__m256i *avx_inv_p_modMa ;
	int64_t *p_modMb;      // P mod M2
	__m256i *avx_p_modMb ;
	int64_t *inv_Ma_modMb;  // M^-1 mod M2
	__m256i *avx_inv_Ma_modMb; ;
	struct conv_base_t *conv; //Constants for the RNSa -> RNSb conversion
};
#endif






