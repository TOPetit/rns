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
inline void add_rns(int64_t *rop, struct rns_base_t *base, int64_t *pa, int64_t *pb){
	int j;
	int128 tmp;
	
	for (j=0; j<base->size; j++)
	{
		tmp = (int128)pa[j] + pb[j];
		rop[j] = (int64_t)(tmp % base->m[j]);
	}
}

inline void add_rns_cr(int64_t *rop, struct rns_base_t *base, int64_t *pa, int64_t *pb){
	int j;
	// int128 tmp;
	// int128 up;
	int64_t tmp, up, lo, mask = ((int64_t)1<<63) -1;  /////////////////////////////////////
	
	for (j=0; j<base->size; j++)
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
inline void sub_rns(int64_t *rop, struct rns_base_t *base, int64_t *pa, int64_t *pb){
	int j;
	int128 tmp;
	
	for (j=0; j<base->size; j++)
	{
		tmp = (int128)pa[j] - pb[j];
		rop[j] = (int64_t)(tmp % base->m[j]);
	}
}

inline void sub_rns_cr(int64_t *rop, struct rns_base_t *base, int64_t *pa, int64_t *pb){
	int j;
	int128 tmp;
	int128 up;

	for (j=0; j<base->size; j++)
	{
		tmp = (int128)pa[j] - pb[j];
		up = tmp >> 63 ; //////////////////////////////////////////////////////

		rop[j] = (int64_t)(tmp + up*base->k[j]);
	}
}

///////////////////////////////
// RNS multiplication
///////////////////////////////
// rop : result
// base : RNS base 
// pa : A
// pb : B
inline void mul_rns(int64_t *rop, struct rns_base_t *base, int64_t *pa, int64_t *pb){
	int j;
	int128 tmp;
	
	for (j=0; j<base->size; j++)
	{
		tmp = (int128)pa[j] * pb[j];
		rop[j] = (int64_t)(tmp % base->m[j]);
	}
}


inline void mul_rns_cr(int64_t *rop, struct rns_base_t *base, int64_t *pa, int64_t *pb){
	int j;

	for (j=0; j<base->size; j++)
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
	int64_t *pbb, struct mod_mul_t *mult, int64_t *tmp[3]){

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

	mul_rns(tmp[2] , mult->conv->rns_b, tmp[0], mult->p_modMb); // Q*P base2

// from_rns_to_int_crt(t, mult->conv->rns_b, tmp[2]);
// gmp_printf("QP in base2 %Zd\n", t);

	add_rns(tmp[0] , mult->conv->rns_b, tmp[1], tmp[2]); // A*B + Q*P in base 2

// from_rns_to_int_crt(t, mult->conv->rns_b, tmp[0]);
// gmp_printf("AB + QP in base2 %Zd\n", t);

	mul_rns(rop, mult->conv->rns_b, tmp[0], mult->inv_Ma_modMb); // Division by Ma

// from_rns_to_int_crt(t, mult->conv->rns_b, tmp[0]);
// gmp_printf("res in base2 %Zd\n", t);

}

void mult_mod_rns_cr(int64_t *rop, int64_t *pa, int64_t *pab, int64_t *pb, 
	int64_t *pbb, struct mod_mul_t *mult, int64_t *tmp[4]){

	int i;

	mul_rns_cr(tmp[0], mult->conv->rns_a, pa, pb); //A*B
	mul_rns_cr(tmp[1], mult->conv->rns_b, pab, pbb); //A*B in base2

	// base_conversion_cr(tmp[1], mult->conv, tmp[0], tmp[3]); // A*B in base2
	
	mul_rns_cr(tmp[2], mult->conv->rns_a, tmp[0], mult->inv_p_modMa); //Q*{P-1}

	base_conversion_cr(tmp[0], mult->conv, tmp[2], tmp[3]); // Q in base2
	mul_rns_cr(tmp[2] , mult->conv->rns_b, tmp[0], mult->p_modMb); // Q*P base2
	add_rns_cr(tmp[0] , mult->conv->rns_b, tmp[1], tmp[2]); // A*B + Q*P in base 2
	mul_rns_cr(rop, mult->conv->rns_b, tmp[0], mult->inv_Ma_modMb); // Division by Ma

}

///////////////////////////////
// GMP to RNS convertion
///////////////////////////////
// TODO : NB_COEFF to NB_RES
//~ Assumes allocation already done for "rop".
void from_int_to_rns(int64_t *rop, struct rns_base_t *base, mpz_t op){
	int i;
	mpz_t tmp_residue;

//	printf("On rentre ds la conversion\n");
	
	mpz_init(tmp_residue);
	if(op->_mp_size == 0)
		return;

//	printf("avant la boucle \n");
	
	for(i=0; i<base->size; i++){
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

	Mi = (mpz_t*)malloc(base->size*sizeof(mpz_t));
	inv_Mi = (mpz_t*)malloc(base->size*sizeof(mpz_t));
	for(i=0; i<base->size;i++)
	{
		mpz_init(Mi[i]);
		mpz_init(inv_Mi[i]);
	}	
	mpz_init(base->M);
	mpz_init(tmp_gcd);
	mpz_init(tmp_mi);
	mpz_init(t);

	// Computes M
	mpz_add_ui(base->M, base->M, 1);
	for(i=0; i<base->size; i++){         
		mpz_mul_ui(base->M, base->M, base->m[i]);
	}
    // Computes Mi and inv_Mi
	for(i=0; i<base->size; i++){         
		mpz_fdiv_q_ui(Mi[i], base->M, base->m[i]);
		mpz_set_ui(tmp_mi, base->m[i]);
		mpz_gcdext(tmp_gcd, inv_Mi[i], t, Mi[i], tmp_mi);  
	}
	base->Mi = Mi;
	base->inv_Mi = inv_Mi;
	mpz_clears(tmp_gcd, tmp_mi, t, NULL);
	
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
	for(i=0; i<base->size;i++)
	{
		mpz_clear(base->Mi[i]);
		mpz_clear(base->inv_Mi[i]);
	}
	free(base->Mi); free(base->inv_Mi);
}
/////////////////////////////////
// Convert RNS to GMP using 
// the CRT
// Assumes "rop" already initialized
/////////////////////////////////
// TODO : NB_COEFF to NB_RES	
void from_rns_to_int_crt(mpz_t rop, struct rns_base_t *base, int64_t *op){
	mpz_t tmp; 
	int i;
	
	// Initializations
	mpz_init(tmp);
	mpz_set_ui(rop,0);
	
	for(i=0; i<base->size; i++){ 
	    mpz_mul_ui(tmp, base->inv_Mi[i], op[i]);   //xi*(Mi**(-1))
		mpz_fdiv_r_ui(tmp, tmp, base->m[i]);       //xi*(Mi**(-1)) mod mi
		mpz_mul(tmp, tmp, base->Mi[i]);
		mpz_add(rop, rop, tmp);		
	}
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
    int64_t gcd = gcdExtended(b%a, a, &x1, &y1); 
  
    // Update x and y using results of recursive 
    // call 
    *x = y1 - (b/a) * x1; 
    *y = x1; 
  
    return gcd; 
} 

///////////////////////////////////////////////////////
// Initializes the constant inversese for 
// the base conversion using MRS
// inv and mrs are supposed to be already initialized
///////////////////////////////////////////////////////
void initialize_inverses_base_conversion(struct conv_base_t *conv_base){
//int64_t **inv, int64_t **mrs, int64_t *m1, int64_t *m2){
	int i, j;
	int64_t tmp, x;
	int128 tmp2;
//	int64_t res;
//	int128 prod;
	
	int size = conv_base->rns_a->size;

	// Memory allocation for mrs and inverse
	conv_base->inva_to_b = (int64_t **) malloc(size*sizeof(int64_t*));
	conv_base->mrsa_to_b = (int64_t **) malloc(size*sizeof(int64_t*));		
	for(i=0; i<NB_COEFF; i++)
	{
		conv_base->inva_to_b[i] = (int64_t *) malloc(size*sizeof(int64_t));
		conv_base->mrsa_to_b[i] = (int64_t *) malloc(size*sizeof(int64_t));	
	}
    // Initialization of the arrays for mrs and inverse
	for(i=0; i<size; i++)  
	{
		for(j=0; j<size; j++)
		{
			conv_base->inva_to_b[i][j]=0;
			conv_base->mrsa_to_b[i][j]=1;
		}
	}
	// Modular inverses : inv[i][j] = inv (m_i) mod m_j
	for(i=0; i<size-1; i++)
	{
		for(j = i+1; j<size; j++)
		{
			gcdExtended(conv_base->rns_a->m[i], conv_base->rns_a->m[j], &x, &tmp);	
			if (x >0)		// Sinon ca part en couille avec 128 et modulo
				conv_base->inva_to_b[i][j]=x; 
			else
				conv_base->inva_to_b[i][j]=x + conv_base->rns_a->m[j];
//			printf(" %ld  %ld\n", inv[i][j], tmp);
			
//			prod = ((int128)m1[i]*inv[i][j]);
//			res = prod%m2[j];
//			printf("prod %ld  tmp %ld\n", res, tmp);
			
//			printf("Inverses m%d = %ld modulo m'%d = %ld: gcd %ld\n",  i, m1[i], j, m2[j], g);
//			printf(" %ld  %ld\n\n", x, tmp);
		}
	}
	// Modular value of the mrs base
	for(j = 0; j<size; j++)
		conv_base->mrsa_to_b[0][j] = conv_base->rns_a->m[0] % conv_base->rns_b->m[j];
	for(i=1; i<size-1; i++)
		for(j = 0; j<size; j++)
		{	tmp2 = (int128)conv_base->mrsa_to_b[i-1][j] * conv_base->rns_a->m[i];
			conv_base->mrsa_to_b[i][j] = tmp2 % conv_base->rns_b->m[j];
		}
		
//	for(i=0; i<NB_COEFF-1; i++)
//		for(j = 0; j<NB_COEFF; j++)
//			printf("mrs[%d, %d] = %ld\n",i ,j, mrs[i][j]); 
}

///////////////////////////////////////////////////////
// Converts a RNS number from a base into an other
// using the MRS conversion
///////////////////////////////////////////////////////
void base_conversion(int64_t *rop, struct conv_base_t *conv_base, int64_t *op){	
	int i, j;
	int64_t a[NB_COEFF*10]; ///////////////////////////ARGH !! Le coût de ce truc n'est pas maîtrisé !
	int64_t tmp;
	int128 tmp2;
	int size = conv_base->rns_a->size;

	// Set target number to 0
	// for(j=0; j<NB_COEFF; j++)
	// 	rop[j]=0;	
	for(j=0; j<size; j++) 
		a[j]=op[j];
	for(i=0; i<size-1; i++) 
	{
		for(j = i+1; j<size; j++)
		{
			tmp = a[j]-a[i];
			tmp2 = (int128)tmp*conv_base->inva_to_b[i][j];
			a[j] = tmp2 % conv_base->rns_a->m[j];
			if(a[j]<0)	// Sinon ca part en couille ????????????
				a[j]+=conv_base->rns_a->m[j];
		}

	}
	
	// Residue of the MRS radix
	for(j=0; j<size; j++)
		rop[j]=a[0] % conv_base->rns_b->m[j];	
	for(j=0; j<size; j++)
	{
		for(i=1; i<size ; i++)
		{
			tmp2 = (int128)a[i]*conv_base->mrsa_to_b[i-1][j];
			tmp = tmp2 % conv_base->rns_b->m[j];
			rop[j]=((int128)rop[j] + tmp) % conv_base->rns_b->m[j];  //Overflow risk

			if (rop[j]<0)	//Sinon, ca part en couille  ??????
				rop[j] += conv_base->rns_b->m[j];
		}

	} 
	
}

///////////////////////////////////////////////////////
// Converts a RNS number from a base into an other
// using the MRS conversion. The RNS base uses Crandall 
// numbers
///////////////////////////////////////////////////////
void base_conversion_cr(int64_t *rop, struct conv_base_t *conv_base, int64_t *op, int64_t *a){	
	int i, j;
//	int64_t a[NB_COEFF];  // En externe, car ça prend du temps 
	int64_t tmp;
	int128 tmp2;
	int128 tmp3;
	int64_t up, up2, lo, lo2;
	int64_t mask = ((int64_t)1<<63) -1;  /////////////////////////////////////
	int size = conv_base->rns_a->size;

	
	// Set target number to 0
	// for(j=0; j<NB_COEFF; j++)
	// 	rop[j]=0;	
	for(j=0; j<size; j++)   // Indispensable
		a[j]=op[j];
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
	for(j=0; j<size; j++)
		rop[j]=a[0] % conv_base->rns_b->m[j];	
	for(j=0; j<size; j++)
	{
		for(i=1; i<size ; i++)
		{
			tmp = mul_mod_cr(a[i], conv_base->mrsa_to_b[i-1][j], conv_base->rns_b->k[j]);
			rop[j]= add_mod_cr(rop[j], tmp, conv_base->rns_b->k[j]);

			if (rop[j]<0){	//Sinon, ca part en couille ?????????????
				rop[j] += conv_base->rns_b->m[j];
				printf("conv");
			}
		}
	} 
}

///////////////////////////////////////////////////////
// Converts a RNS number from a base into an other
// using the Shenoy - Kumaresan approximate conversion
// IN DEVELOPMENT
///////////////////////////////////////////////////////

void base_conversion_sk(int64_t *rop, struct conv_base_t *conv_base, int64_t *op, int64_t *a){	
	int i, j;
//	int64_t a[NB_COEFF];  // En externe, car ça prend du temps 
	int64_t tmp;
	int128 tmp2;
	int128 tmp3;
	int64_t up, up2, lo, lo2;
	int64_t mask = ((int64_t)1<<63) -1;  /////////////////////////////////////
	int size = conv_base->rns_a->size;


	int64_t xi[NB_COEFF];
	int64_t sigma=0;		// Sigma_0
	int64_t t1, t2, t, q;
	mul_rns(xi, conv_base->rns_a, op, conv_base->mrsa_to_b[0]);  // Juste pour le LOL
	for(i=0; j<size; i++){
		t = xi[i]/conv_base->rns_a->m[i];
		sigma += t;
		q=sigma;
		sigma -=q;
		for(j=0; j<size; j++){
			t1 = mul_mod_cr(xi[i],conv_base->mrsa_to_b[1][j], conv_base->rns_b->k[i]); // Juste pour le LOL
			rop[j] = add_mod_cr(rop[j], t1, conv_base->rns_b->k[i]);
			t1= mul_mod_cr(q,conv_base->mrsa_to_b[1][j], conv_base->rns_b->k[i]); // Juste pour le LOL
			rop[j] = add_mod_cr(rop[j], t1, conv_base->rns_b->k[i]);
		}
	}


// 	for(j=0; j<size; j++)
// 	// Set target number to 0
// 	// for(j=0; j<NB_COEFF; j++)
// 	// 	rop[j]=0;	
// 	for(j=0; j<size; j++)   // Indispensable
// 		a[j]=op[j];
// 	for(i=0; i<size-1; i++) 
// 	{
// 		for(j = i+1; j<size; j++)
// 		{
// 			tmp = a[j]-a[i];
// 			a[j] = mul_mod_cr(tmp, conv_base->inva_to_b[i][j],  conv_base->rns_a->k[j]);
// 			if(a[j]<0)	// Sinon ca part en couille ?
// 				a[j]+=conv_base->rns_a->m[j];
// 		}

// 	}
	
// 	// Residue of the MRS radix
// 	for(j=0; j<size; j++)
// 		rop[j]=a[0] % conv_base->rns_b->m[j];	
// 	for(j=0; j<size; j++)
// 	{
// 		for(i=1; i<size ; i++)
// 		{
// 			tmp = mul_mod_cr(a[i], conv_base->mrsa_to_b[i-1][j], conv_base->rns_b->k[j]);
// 			rop[j]= add_mod_cr(rop[j], tmp, conv_base->rns_b->k[j]);

// 			if (rop[j]<0){	//Sinon, ca part en couille ?????????????
// 				rop[j] += conv_base->rns_b->m[j];
// 				printf("conv");
// 			}
// 		}
// 	} 
}

///////////////////////////////////////////////////
// Modular addition and multiplication using 
// Crandall moduli
///////////////////////////////////////////////////
int64_t add_mod_cr(int64_t a, int64_t b, int k)   
{
	int64_t tmp, up, lo, mask = ((int64_t)1<<63) -1; /////////////////////////////////////
	int64_t res=0;  // Required, else computation on 32 bits with a sign extension to 64 bits
	uint64_t u_res, u_mod;

	u_mod = ((int64_t)1<<63)-k;  // only 32 shifts without int64_t extension
	tmp = a + b;
	up = (uint64_t)tmp >> 63 ; //////////////////////////////////////////////////////
								/// Unsigned conversion to avoid sign extension
	lo = tmp & mask;
	res += lo + up*k;  
	u_res = (uint64_t) res;
//	if ( u_res  >  u_mod	){	// For the unfrequent case :  2^63 > res >=mod ????
//		res -=u_mod;
//		printf("a");
//	}
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
	int64_t mask = ((int64_t)1<<63) -1;  /////////////////////////////////////
	int64_t res=0;  // Required, else computation on 32 bits with a sign extension to 64 bits
	uint64_t u_res, u_mod;

	u_mod = ((int64_t)1<<63)-k;  // only 32 shifts without int64_t extension
	prod = (int128)a * b;
	up = (int64_t) ((uint128)prod >> 63); //////////////////////////////////
	lo = (int64_t)prod & mask;

	tmp = (int128)up * k;
	lo2 = (int64_t)tmp & mask;
	up2 = (int64_t)((uint128)tmp >> 63); ///////////////////////////////////
	tmp2 = (uint128)lo + lo2 + up2*k;
	up3 = (uint64_t)(tmp2 >> 63); ///////////////////////////////////
	lo3 = (uint64_t)tmp2 & mask;

	res += lo3 + up3*k;

	u_res = (uint64_t) res;
// printf("i");
	if ( u_res  >  u_mod ){	// For the unfrequent case :  2^63 > res >=mod ????
		res -=u_mod;
		// printf("m");
	}
	return res;
}