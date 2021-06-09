#include <stdlib.h>
#include <stdio.h>
//#include <stdint.h>
#include <gmp.h>

//#include "structs_data.h"
//#include "print_functs.c"
//#include "red_functs.c"
//#include "add_mult_poly.c"
//#include "conv_functs.c"


//#include "amns_check.c"
//#include "tests.c"

#include "rns.h"
#include "tests.c"

//~ used in poly prod and internal reduction
extern int128 tmp_prod_result[NB_COEFF];

//~ used in the internal reductions
extern char tmp_sign_coeff_poly[NB_COEFF];

//~ used in the the little internal reduction
extern int64_t tmp_red_int_poly_low[NB_COEFF];
extern int64_t tmp_red_int_poly_high[NB_COEFF];

//~ used in the total internal reduction
extern int64_t tmp_red_coeff_poly_high[NB_COEFF];
extern int128  tmp_red_coeff_poly_low[NB_COEFF];

extern uint64_t red_int_mask;


//~ The module p
mpz_t modul_p;

 
int main(void){
	
	mpz_t A, B, C, E, F, inv_p_modM, inv_M_modMp;
	mpz_t Arns, Brns;
	mpz_inits (A, B, C, E, F, inv_p_modM, inv_M_modMp, NULL);
	mpz_inits (Arns, Brns, NULL);
	
	int64_t pa[NB_COEFF];
	int64_t pb[NB_COEFF];
	int64_t pab[NB_COEFF];
	int64_t pbb[NB_COEFF];
	int64_t pc[NB_COEFF];
	int64_t pp1[NB_COEFF];
	int64_t pp2[NB_COEFF];
	int64_t pp3[NB_COEFF];
	
	int i;



	// mpz_set_str (A, "115792089021636622621247151603347568778042451706330200410359523598128907", 10);
	// // mpz_set_str (A, "1", 10);
	// // mpz_set_str (A, "-115792089021636622621247151603347568778042451706330200410359523598128907", 10);
	// // mpz_set_str (B, "1", 10);
	// mpz_set_str (B, "451258902198277337595903956877804245177045777827200410351952359812890591", 10);

	mpz_set_str (A, "100106136745050507346481674824002435247796236765700289941745547607771451581114", 10); 
	mpz_set_str (B, "14037175231330152997227542512191303093263562666109397624121151367256237560831", 10);

    mpz_set_str (modul_p, "115792089021636622262124715160334756877804245386980633020041035952359812890593", 10);

	mpz_mul (E, A, B);
	mpz_mod (E, E, modul_p);
	gmp_printf("A   : %Zd\n", A);
	gmp_printf("B   : %Zd\n", B);
	gmp_printf("P   : %Zd\n", modul_p);
	
	mpz_set_str (inv_p_modM, "-7210642370083763919688086698199040857322895088554003933210287226647459666846134833419938084604981461493089686639677942359747717700454441525223348684285", 10);
		
	mpz_set_str (inv_M_modMp, "2926906825829426928727294150364906856635623568440932569450673109926460590684432927230290255276608760237299661987870702836538185953568700154975953006659", 10);
	

	////////////////////////////////////////////////////////////////////////
	//Test RNS
    ////////////////////////////////////////////////////////////////////////

	printf("RNS results : \n");

    // Initialization

    struct rns_base_t rns_a;
	int64_t base[NB_COEFF]={ 9223372036854775783, 9223372036854775643, 
				 9223372036854775549, 9223372036854775507, 
				 9223372036854775433, 9223372036854775421, 
				 9223372036854775417, 9223372036854775399};
    rns_a.size = NB_COEFF;
	rns_a.m = base;
    init_rns(&rns_a);

    // Conversion and display

	gmp_printf("B to be converted   : %Zd\n", B);
	from_int_to_rns(pa, &rns_a, B);
	mpz_set_ui(C,0);
	from_rns_to_int_crt(C, &rns_a, pa);
	gmp_printf("B converted  : %Zd\n", C);
	
    // Base extension

	struct rns_base_t rns_b;
	int64_t base2[NB_COEFF]={ 9223372036854775351, 9223372036854775337, 
				 9223372036854775291, 9223372036854775279, 
				 9223372036854775259, 9223372036854775181, 
				 9223372036854775159, 9223372036854775139};
    rns_b.size = NB_COEFF;
	rns_b.m = base2;
    init_rns(&rns_b);
	
	struct conv_base_t conv;
	conv.rns_a = &rns_a;
	conv.rns_b = &rns_b;
	initialize_inverses_base_conversion(&conv);			
	base_conversion(pb, &conv, pa);		
	from_rns_to_int_crt(C, &rns_b, pb);
	gmp_printf("B ext: %Zd\n", C);

    // Base extension crandall

    int64_t base_bis[NB_COEFF]={9223372036854775805,
9223372036854775801,
9223372036854775789,
9223372036854775783,
9223372036854775777,
9223372036854775769,
9223372036854775757,
9223372036854775747};

	int64_t ttt[NB_COEFF];
	int k_bis[NB_COEFF] = {3, 7, 19, 25, 31, 39, 51, 61};
    rns_a.size = NB_COEFF;
	rns_a.m = base_bis;
	rns_a.k = k_bis;
    init_rns(&rns_a);
	int64_t base2_bis[NB_COEFF]={ 9223372036854775807,
9223372036854775803,
9223372036854775799,
9223372036854775787,
9223372036854775781,
9223372036854775771,
9223372036854775763,
9223372036854775753};
	int k2_bis[NB_COEFF] = {1, 5, 9, 21, 27, 37, 45, 55};
    rns_b.size = NB_COEFF;
	rns_b.m = base2_bis;
	rns_b.k = k2_bis;
    init_rns(&rns_b);


	conv.rns_a = &rns_a;
	conv.rns_b = &rns_b;
	initialize_inverses_base_conversion(&conv);		

	gmp_printf("\nRNS Crandall results : \nB to be converted   : %Zd\n", B);
	from_int_to_rns(pa, &rns_a, B);
	from_rns_to_int_crt(C, &rns_a, pa);
	gmp_printf("B converted  : %Zd\n", C);

	base_conversion(pb, &conv, pa);		
	from_rns_to_int_crt(C, &rns_b, pb);
	gmp_printf("B ext: %Zd\n", C);

	base_conversion_cr(pb, &conv, pa, ttt);		
	from_rns_to_int_crt(C, &rns_b, pb);
	gmp_printf("B ext crandall: %Zd\n", C);


	int64_t ta, tb, tc, tk, plop;
	uint128 mod, tmpo;
	uint128 ta_l, tb_l;

	for(i=0; i< (1<<20); i++)
	{
		ta=rand();   // Rand() produces an integer in [0, 2^30[
		ta_l = ta<<32;
		ta_l+= rand();
		tb=rand();
		tb_l = tb<<32;
		tb_l+= rand();
	
		tk=rand()%128 + 1;

		mod = ((uint128)1<<63)-tk;
		tmpo = ta_l+tb_l;
		tc = (int64_t)(tmpo % mod);
		plop = add_mod_cr(ta_l, tb_l, tk);  
		if(tc!=plop)
			printf("%lu + %lu mod %lu = %lu  mais %lu  %lu\n", ta, tb, mod, tc, plop, tk);

	}

	printf(" test of the modular addition Done\n");


	for(i=0; i< (1<<20); i++)
	{
		ta=rand();   // Rand() produces an integer in [0, 2^30[
		ta_l = ta<<32;
		ta_l+= rand();
		tb=rand();
		tb_l = tb<<32;
		tb_l+= rand();
	
		tk=rand()%64 + 1;

		mod = ((uint128)1<<63)-tk;
		tmpo = ta_l * tb_l;
		tc = (int64_t)(tmpo % mod);
		plop = mul_mod_cr(ta_l, tb_l, tk);  
		if(tc!=plop)
			printf("%lu * %lu mod %lu = %lu  mais %lu  %lu\n", ta_l, tb_l, mod, tc, plop, tk);

	}

	printf(" test of the modular multiplication Done\n");

	gmp_randstate_t st;
	gmp_randinit_default (st);
 	mpz_t AA, BB;
 	mpz_inits (AA, BB, NULL);
	for(i=0; i< (1<<15); i++)
	{
		mpz_urandomm (AA, st, modul_p);  //Randomly generates A < P
		// mpz_urandomm (BB, st, modul_p);  //Randomly generates B < P
		from_int_to_rns(pa, &rns_a, AA);

		base_conversion_cr(pb, &conv, pa, ttt);		
		from_rns_to_int_crt(BB, &rns_b, pb);

 		if(mpz_cmp(AA,BB))
 			gmp_printf("B ext crandall: %Zd\n", AA);

	}
	printf("\n test of the base extension done\n");


	//Modular multiplication

    struct mod_mul_t mult;
	mpz_t tmp_gcd, t, tmp_inv; 

	mpz_init(tmp_gcd);
	mpz_init(t);
	mpz_init(tmp_inv);
	from_int_to_rns(pa, &rns_a, A);
	from_int_to_rns(pb, &rns_a, B);
	from_int_to_rns(pab, &rns_b, A);
	from_int_to_rns(pbb, &rns_b, B);
	from_int_to_rns(pp2, &rns_b, modul_p);      // P mod Mb

	mpz_sub (tmp_inv, rns_a.M, modul_p);
	mpz_gcdext(tmp_gcd, inv_p_modM, t, tmp_inv, rns_a.M);  
	from_int_to_rns(pp1, &rns_a, inv_p_modM);   //(-P)^-1 mod Ma

	mpz_gcdext(tmp_gcd, inv_M_modMp, t, rns_a.M, rns_b.M);  
	from_int_to_rns(pp3, &rns_b, inv_M_modMp);  // Ma^{-1} mod Mb

	mult.inv_p_modMa = pp1;
	mult.p_modMb = pp2; 
	mult.inv_Ma_modMb = pp3;
	mult.conv = &conv;

	from_rns_to_int_crt(Arns, &rns_a, pa);
	from_rns_to_int_crt(Brns, &rns_a, pb);
	from_rns_to_int_crt(AA, &rns_b, pp3);  // M^-1

	// gmp_printf("Arns %Zd\n", Arns);
	// gmp_printf("Brns %Zd\n", Brns);
	gmp_printf("inv_M_modMp RNS %Zd\n", AA);
	// gmp_printf("-P mod M %Zd\n", tmp_inv);
	// from_rns_to_int_crt(AA, &rns_a, pp1);  // P^-1
	// gmp_printf("inv_-p_modMa RNS %Zd\n", AA);
	gmp_printf("M %Zd\n", rns_a.M);
	// gmp_printf("M2 %Zd\n", rns_b.M);

    int64_t *tmp[4];  // RNS modular multiplication intermediate results
     				  // One more for the base convertion
    tmp[0]=(int64_t *)malloc(NB_COEFF*sizeof(int64_t));
    tmp[1]=(int64_t *)malloc(NB_COEFF*sizeof(int64_t));
    tmp[2]=(int64_t *)malloc(NB_COEFF*sizeof(int64_t));
    tmp[3]=(int64_t *)malloc(NB_COEFF*sizeof(int64_t));
 

 printf("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n");

    mult_mod_rns(pc, pa, pab, pb, pbb, &mult, tmp);
	from_rns_to_int_crt(F, &rns_b, pc);
	gmp_printf("F  %Zd\n", F);

    mult_mod_rns_cr(pc, pa, pab, pb, pbb, &mult, tmp);
	from_rns_to_int_crt(F, &rns_b, pc);
	gmp_printf("F Crandall %Zd\n", F);

	// // Verif 
	mpz_t MM; 
	mpz_init(MM);

	mpz_gcdext(tmp_gcd, MM, t, rns_a.M, modul_p);  
	mpz_mul (E, A, B);
	mpz_mul (E, E, MM);
	mpz_mod (E, E, modul_p);
	gmp_printf("F GMP %Zd\n", E);

 	mpz_t CC, DD;
 	mpz_inits (CC, DD, NULL);

	for(i=0; i< (1<<15); i++)
	{
		mpz_urandomm (AA, st, modul_p);  //Randomly generates A < P
		mpz_urandomm (BB, st, modul_p);  //Randomly generates B < P
		from_int_to_rns(pa, &rns_a, AA);
		from_int_to_rns(pb, &rns_a, BB);
		from_int_to_rns(pab, &rns_b, AA);
		from_int_to_rns(pbb, &rns_b, BB);

    	mult_mod_rns_cr(pc, pa, pab, pb, pbb, &mult, tmp);
		from_rns_to_int_crt(CC, &rns_b, pc);

		mpz_mul (E, AA, BB);
		mpz_mul (E, E, MM);
		mpz_mod (DD, E, modul_p);

 		if(mpz_cmp(CC,DD))
 			gmp_printf("RNS Mod mult: %Zd %Zd = %Zd  vs  %Zd \n", AA, BB, CC, DD);

	}
	printf("\n test of the RNS modular multiplication done\n");


	////////////////////////////////////////////////////////////////////////
	// Timing tests
	////////////////////////////////////////////////////////////////////////
	unsigned long long timer , meanTimer =0, t1, t2;
	unsigned long long timer2 , meanTimer2 =0;	
	unsigned long long timer3 , meanTimer3 =0;	
	unsigned long long nb1, nb2, nb_ins;
	gmp_randstate_t state;
	gmp_randinit_default (state);
 

    //////////////////////
	// RNS  timing tests
    //////////////////////
	
	meanTimer =0;
	meanTimer2 =0;
	// Heating caches
	for(int i=0;i<NTEST;i++)
	{
		// appel de la fonction a mesurer a mettre ici
		// juste pour chauffer les caches
		mult_mod_rns(pc, pa, pab, pb, pbb, &mult, tmp);

	}
  
	// timing
	for(int i=0;i<NSAMPLES;i++)
	{

     // initialiser un nouveau jeu de donnees a tester
		mpz_urandomm (A, state, modul_p);  //Randomly generates A < P
		mpz_urandomm (B, state, modul_p);  //Randomly generates B < P
		from_int_to_rns(pa, &rns_a, A);
		from_int_to_rns(pb, &rns_a, B);
		from_int_to_rns(pab, &rns_b, A);
		from_int_to_rns(pb, &rns_b, B);
		timer = (unsigned long long int)0x1<<63;
		timer2 = (unsigned long long int)0x1<<63;
		for(int j=0;j<NTEST;j++)
		{
			t1 = cpucyclesStart();

			// t1 = rdpmc_actual_cycles();
			// nb1 = rdpmc_instructions();

			// appeler la fonction ici avec toujours le meme jeu de donnees
			mult_mod_rns(pc, pa, pab, pb, pbb, &mult, tmp);

			// nb2= rdpmc_instructions();
			// t2 = rdpmc_actual_cycles();
			t2 = cpucyclesStop();
			//if(timer2>nb2-nb1) timer2 = nb2-nb1;
			if(timer>t2-t1) timer = t2-t1;
		}
		
		meanTimer += timer;
		//meanTimer2 += timer2;

		//printf("%lld \n", meanTimer3);
	}

 //   printf("\n RNS Modular multiplication :  %lld instruction %lld CPU cycles %f IPC \n", meanTimer2/NSAMPLES, meanTimer/NSAMPLES, (double)meanTimer2/meanTimer);
    printf("\n RNS Modular multiplication : %lld CPU cycles  \n", meanTimer/NSAMPLES);

 

    ///////////////////////////////////////////////////////////////////
    // Crandall RNS base
    ///////////////////////////////////////////////////////////////////

       // Initialization

    struct rns_base_t rns_c;
	int64_t base3[NB_COEFF]={ 9223372036854775807,9223372036854775803,
		9223372036854775799,9223372036854775795,9223372036854775787,
		9223372036854775783,9223372036854775771,9223372036854775763};
	int k3[NB_COEFF] = {1, 5, 9, 13, 21, 25, 37, 45};
    rns_c.size = NB_COEFF;
	rns_c.m = base3;
	rns_c.k = k3;
    init_rns(&rns_c);

    struct rns_base_t rns_d;
	int64_t base4[NB_COEFF]={ 9223372036854775805,9223372036854775801,
		9223372036854775797,9223372036854775793,9223372036854775789,
		9223372036854775781,9223372036854775777,9223372036854775769};
	int k4[NB_COEFF] = { 3, 7, 11, 15, 19, 27, 31, 39};
    rns_d.size = NB_COEFF;
	rns_d.m = base4;
	rns_d.k = k4;
    init_rns(&rns_d);

    // New constants for these new bases
	mpz_set_str (inv_p_modM, "21693078482973340509982853143798430575319512745533392255515620251574074455759948281164430134922732925648364830157422883443108685399968728543675726190767", 10);
		
	mpz_set_str (inv_M_modMp, "129815111504256952841543163587079826784843206503137094084041764328371834023366811606378043065985956427717294679212436277317010435733086197442429103533", 10);
	struct conv_base_t conv2;

	conv2.rns_a = &rns_c;
	conv2.rns_b = &rns_d;
	initialize_inverses_base_conversion(&conv2);	

	from_int_to_rns(pa, &rns_c, A);
	from_int_to_rns(pb, &rns_c, B);
	from_int_to_rns(pab, &rns_d, A);
	from_int_to_rns(pbb, &rns_d, B);
	from_int_to_rns(pp1, &rns_c, inv_p_modM);   //(-P)^-1 mod Ma
	from_int_to_rns(pp2, &rns_d, modul_p);      // P mod Mb
	from_int_to_rns(pp3, &rns_d, inv_M_modMp);  // Ma^{-1} mod Mb

	mult.inv_p_modMa = pp1;
	mult.p_modMb = pp2;
	mult.inv_Ma_modMb = pp3;
	mult.conv = &conv2;

    meanTimer =0;
	meanTimer2 =0;
	// Heating caches
	for(int i=0;i<NTEST;i++)
	{
		// appel de la fonction a mesurer a mettre ici
		// juste pour chauffer les caches
		mult_mod_rns_cr(pc, pa, pab, pb, pbb, &mult, tmp);

	}
  
	// timing
	for(int i=0;i<NSAMPLES;i++)
	{

     // initialiser un nouveau jeu de donnees a tester
		mpz_urandomm (A, state, modul_p);  //Randomly generates A < P
		mpz_urandomm (B, state, modul_p);  //Randomly generates B < P
		from_int_to_rns(pa, &rns_c, A);
		from_int_to_rns(pb, &rns_c, B);
		from_int_to_rns(pab, &rns_d, A);
		from_int_to_rns(pbb, &rns_d, B);
		timer = (unsigned long long int)0x1<<63;
		timer2 = (unsigned long long int)0x1<<63;
		for(int j=0;j<NTEST;j++)
		{
			t1 = cpucyclesStart();

			// t1 = rdpmc_actual_cycles();
			// nb1 = rdpmc_instructions();
			// appeler la fonction ici avec toujours le meme jeu de donnees
			mult_mod_rns_cr(pc, pa, pab, pb, pbb, &mult, tmp);

			// nb2= rdpmc_instructions();
			// t2 = rdpmc_actual_cycles();
			t2 = cpucyclesStop();
			//if(timer2>nb2-nb1) timer2 = nb2-nb1;
			if(timer>t2-t1) timer = t2-t1;
		}
		
		meanTimer += timer;
		//meanTimer2 += timer2;

		//printf("%lld \n", meanTimer3);
	}

 //   printf("\n RNS Crandall Modular multiplication :  %lld instruction %lld CPU cycles %f IPC \n", meanTimer2/NSAMPLES, meanTimer/NSAMPLES, (double)meanTimer2/meanTimer);
    printf("\n RNS Crandall Modular multiplication : %lld CPU cycles  \n", meanTimer/NSAMPLES);



	mpz_clears (A, B, C, E, F, inv_p_modM, modul_p, inv_M_modMp, NULL);

	printf("mpz_clear OK \n");

	free(tmp[0]);
	free(tmp[1]);
	free(tmp[2]);


	return 0;



}







