#include <stdlib.h>

#include <stdio.h>

#include <gmp.h>

#include <math.h>

#include "rns.h"

#include "tests.c"

#include "rnsv.h"

int main(void)
{

	mpz_t A, B, C;
	mpz_inits(A, B, C, NULL);

	mpz_set_str(A, "100106136745050507346481674824002435247796236765700289941745547607771451581114", 10);
	mpz_set_str(B, "14037175231330152997227542512191303093263562666109397624121151367256237560831", 10);

	struct rns_base_t rns_a;

	int64_t base_bis[NB_COEFF] = {
		9223372036854775805,
		9223372036854775801,
		9223372036854775789,
		9223372036854775783,
		9223372036854775777,
		9223372036854775769,
		9223372036854775757,
		9223372036854775747};

	//int64_t ttt[NB_COEFF];
	int k_bis[NB_COEFF] = {
		3,
		7,
		19,
		25,
		31,
		39,
		51,
		61};
	rns_a.size = NB_COEFF;
	rns_a.m = base_bis;
	rns_a.k = k_bis;
	init_rns(&rns_a);

	int64_t tmp_k[NB_COEFF];
	int j;
	for (j = 0; j < NB_COEFF; j++)
	{
		tmp_k[j] = (int64_t)k_bis[j];
	}

	__m256i avx_k[NB_COEFF / 4];
	from_rns_to_m256i(avx_k, &rns_a, tmp_k);
	rns_a.avx_k = avx_k;

	struct rns_base_t rns_b;

	int64_t base2_bis[NB_COEFF] = {
		9223372036854775807,
		9223372036854775803,
		9223372036854775799,
		9223372036854775787,
		9223372036854775781,
		9223372036854775771,
		9223372036854775763,
		9223372036854775753};
	int k2_bis[NB_COEFF] = {
		1,
		5,
		9,
		21,
		27,
		37,
		45,
		55};
	rns_b.size = NB_COEFF;
	rns_b.m = base2_bis;
	rns_b.k = k2_bis;
	init_rns(&rns_b);

	int64_t tmp_k2[NB_COEFF];
	for (j = 0; j < NB_COEFF; j++)
	{
		tmp_k2[j] = (int64_t)k2_bis[j];
	}
	__m256i avx_k2[NB_COEFF / 4];
	from_rns_to_m256i(avx_k2, &rns_b, tmp_k2);
	rns_b.avx_k = avx_k2;

	printf("TEST ADD\n");

	int64_t rop1[NB_COEFF];
	int64_t rop2[NB_COEFF];

	mpz_t testint1, testint2, C2, C3;
	mpz_inits(testint1, testint2, C2, C3, NULL);

	//mpz_set_str(testint1,"150",10);
	//mpz_set_str(testint2,"125",10);

	mpz_set(testint1, A);
	mpz_set(testint2, B);

	gmp_printf("addition des deux entiers %Zd et %Zd\n", testint1, testint2);

	//////////////////////////////////////////////////

	gmp_printf("RNS non vectorisée : \n");

	from_int_to_rns(rop1, &rns_a, testint1);
	from_int_to_rns(rop2, &rns_a, testint2);

	//print_RNS(&rns_a,rop2);

	int64_t result[NB_COEFF];

	add_rns_cr(result, &rns_a, rop1, rop2);

	from_rns_to_int_crt(C2, &rns_a, result);
	gmp_printf("résultat : %Zd\n", C2);

	//////////////////////////////////////////////////

	gmp_printf("RNS vectorisée : \n");

	from_int_to_rns(rop1, &rns_a, testint1);
	from_int_to_rns(rop2, &rns_a, testint2);

	__m256i rop1_aux[NB_COEFF / 4];
	__m256i rop2_aux[NB_COEFF / 4];

	from_rns_to_m256i(rop1_aux, &rns_a, rop1);
	from_rns_to_m256i(rop2_aux, &rns_a, rop2);

	//print_m256i(&rns_a,rop1_aux);
	//print_m256i(&rns_a,rop2_aux);

	int64_t result_avx[NB_COEFF];
	__m256i result_avx_aux[NB_COEFF / 4];

	avx_add_rns_cr(result_avx_aux, &rns_a, rop1_aux, rop2_aux);

	//print_m256i(&rns_a,result_avx_aux);

	from_m256i_to_int64_t_rns(result_avx, &rns_a, result_avx_aux);

	//print_RNS(&rns_a,result_avx);

	from_rns_to_int_crt(C3, &rns_a, result_avx);
	gmp_printf("résultat : %Zd\n", C3);

	///%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%///

	printf("\nTEST SUB\n");

	int64_t rop3[NB_COEFF];
	int64_t rop4[NB_COEFF];

	mpz_t testint3, testint4, S2, S3;
	mpz_inits(testint3, testint4, S2, S3, NULL);

	//mpz_set_str(testint3,"150",10);
	//mpz_set_str(testint4,"125",10); //doesn't work

	mpz_set(testint3, A);
	mpz_set(testint4, B); //problème de segmentation à régler

	gmp_printf("soustraction entre deux entiers : %Zd - %Zd\n", testint3, testint4);

	//////////////////////////////////////////////////

	gmp_printf("RNS non vectorisée : \n");

	from_int_to_rns(rop3, &rns_a, testint3);
	from_int_to_rns(rop4, &rns_a, testint4);

	//print_RNS(&rns_a,rop3);

	int64_t result_sub[NB_COEFF];

	sub_rns_cr(result_sub, &rns_a, rop3, rop4);

	//print_RNS(&rns_a,result_sub);

	from_rns_to_int_crt(S2, &rns_a, result_sub);
	gmp_printf("résultat : %Zd\n", S2);

	//////////////////////////////////////////////////

	gmp_printf("RNS vectorisée : \n");

	from_int_to_rns(rop3, &rns_a, testint3);
	from_int_to_rns(rop4, &rns_a, testint4);

	__m256i rop3_aux[NB_COEFF / 4];
	__m256i rop4_aux[NB_COEFF / 4];

	from_rns_to_m256i(rop3_aux, &rns_a, rop3);
	from_rns_to_m256i(rop4_aux, &rns_a, rop4);

	//print_m256i(&rns_a,rop3_aux);
	//print_m256i(&rns_a,rop4_aux);

	int64_t result_sub_avx[NB_COEFF];
	__m256i result_sub_avx_aux[NB_COEFF / 4];

	avx_sub_rns_cr(result_sub_avx_aux, &rns_a, rop3_aux, rop4_aux);

	//print_m256i(&rns_a,result_sub_avx_aux);

	from_m256i_to_int64_t_rns(result_sub_avx, &rns_a, result_sub_avx_aux);

	from_rns_to_int_crt(S3, &rns_a, result_sub_avx);
	gmp_printf("résultat : %Zd\n", S3);

	///%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%///

	printf("\nTEST MUL\n");

	printf("Test ADD AUX : \n");

	int64_t value = ((int64_t)1 << 62) + 10;
	printf("Avec %ld + %ld \n", value, value);
	__m256i i = _mm256_set1_epi64x(value);
	__m256i up, lo;

	avx_add_aux_2e(&up, &lo, i, i);
	print_alone_m256i(up);
	print_alone_m256i(lo);

	int64_t value2 = ((int64_t)1 << 62) + 10;
	printf("Avec %ld + %ld + %ld \n", value2, value2, value2);
	__m256i i2 = _mm256_set1_epi64x(value2);
	__m256i up2, lo2;

	avx_add_aux_3e(&up2, &lo2, i2, i2, i2);
	print_alone_m256i(up2);
	print_alone_m256i(lo2);

	printf("Test MUL AUX : \n");

	int64_t value3 = ((int64_t)1 << 60) + 1000;

	printf("Avec %ld * %ld \n", value3, value3);
	__m256i i3 = _mm256_set1_epi64x(value3);
	__m256i up3, lo3;

	avx_mul_aux(&up3, &lo3, i3, i3);
	print_alone_m256i(up3);
	print_alone_m256i(lo3);

	int64_t value_a = 6332183518700779750;
	int64_t value_b = 5590294625489548320;

	printf("Avec %ld * %ld \n", value_a, value_b);
	__m256i ia = _mm256_set1_epi64x(value_a);
	__m256i ib = _mm256_set1_epi64x(value_b);
	__m256i up_ab, lo_ab;

	avx_mul_aux(&up_ab, &lo_ab, ia, ib);
	print_alone_m256i(up_ab);
	print_alone_m256i(lo_ab);

	int64_t mask = ((int64_t)1 << 63) - 1;
	int128 prod = (int128)value_a * value_b;
	int64_t up_abt = (int64_t)((uint128)prod >> 63);
	int64_t lo_abt = (int64_t)prod & mask;
	//printf("we should get up --> %ld\n",up_abt);
	//printf("and lo --> %ld\n",lo_abt);

	printf("Test STANDART MUL : \n");

	int64_t rop5[NB_COEFF];
	int64_t rop6[NB_COEFF];

	mpz_t testint5, testint6, K2, K3;
	mpz_inits(testint5, testint6, K2, K3, NULL);

	//mpz_set_str(testint5,"15555555555555555555",10);
	//mpz_set_str(testint6,"1255555555555555555555555",10);

	mpz_set(testint5, A);
	mpz_set(testint6, B);

	gmp_printf("multiplication de deux entiers : %Zd * %Zd\n", testint5, testint6);

	//////////////////////////////////////////////////

	gmp_printf("RNS non vectorisée : \n");

	from_int_to_rns(rop5, &rns_a, testint5);
	from_int_to_rns(rop6, &rns_a, testint6);

	//print_RNS(&rns_a,rop5);

	int64_t result_mul[NB_COEFF];

	mul_rns_cr(result_mul, &rns_a, rop5, rop6);

	//print_RNS(&rns_a,result_mul);

	from_rns_to_int_crt(K2, &rns_a, result_mul);
	gmp_printf("résultat : %Zd\n", K2);

	//////////////////////////////////////////////////

	gmp_printf("RNS vectorisée : \n");

	from_int_to_rns(rop5, &rns_a, testint5);
	from_int_to_rns(rop6, &rns_a, testint6);

	__m256i rop5_aux[NB_COEFF / 4];
	__m256i rop6_aux[NB_COEFF / 4];

	from_rns_to_m256i(rop5_aux, &rns_a, rop5);
	from_rns_to_m256i(rop6_aux, &rns_a, rop6);

	//print_m256i(&rns_a,rop5_aux);
	//print_m256i(&rns_a,rop6_aux);

	int64_t result_mul_avx[NB_COEFF];
	__m256i result_mul_avx_aux[NB_COEFF / 4];

	avx_mul_rns_cr(result_mul_avx_aux, &rns_a, rop5_aux, rop6_aux);

	//print_m256i(&rns_a,result_mul_avx_aux);

	from_m256i_to_int64_t_rns(result_mul_avx, &rns_a, result_mul_avx_aux);

	from_rns_to_int_crt(K3, &rns_a, result_mul_avx);
	gmp_printf("résultat : %Zd\n", K3);

	///%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%///

	printf("\nTEST MOD MUL\n");

	printf("Test BASE CONVERSION \n");

	/*
  int64_t ttt[NB_COEFF];
  int64_t pa[NB_COEFF];
  int64_t pb[NB_COEFF];
  mpz_t Z1, Z2;
  mpz_inits(Z1,Z2,NULL);


  struct conv_base_t conv;
  conv.rns_a = &rns_a;
  conv.rns_b = &rns_b;
  initialize_inverses_base_conversion(&conv);		

  gmp_printf("B to be converted   : %Zd\n", B);

  from_int_to_rns(pa, &rns_a, B);

  base_conversion_cr(pb, &conv, pa, ttt);	
  //print_RNS(&rns_b,pb);	
  from_rns_to_int_crt(Z1, &rns_b, pb);
  gmp_printf("B ext (non vectorisé): %Zd\n", Z1);

  __m256i avx_pa[NB_COEFF/4];
  __m256i avx_pb[NB_COEFF/4];
  from_rns_to_m256i(avx_pa,&rns_a,pa);

  //avx_init_mrs(conv.avx_mrsa_to_b,&conv);
  avx_init_mrs(&conv);printf("dot\n");
  avx_base_conversion_cr(avx_pb, &conv, avx_pa, ttt);	printf("dot\n");
  from_m256i_to_int64_t_rns(pb,&rns_b,avx_pb);	printf("dot2\n");
  //print_RNS(&rns_b,pb);

  from_rns_to_int_crt(Z2, &rns_b, pb);printf("dot2\n");*/

	//TODO : test the modular multiplication algorithm

	/*
  	__m256i *tmp[3];  // RNS modular multiplication intermediate results
       				  // One more for the base convertion
      tmp[0]=(__m256i *)malloc(NB_COEFF*sizeof(__m256i)/4);
      tmp[1]=(__m256i *)malloc(NB_COEFF*sizeof(__m256i)/4);
      tmp[2]=(__m256i *)malloc(NB_COEFF*sizeof(__m256i)/4);

      int64_t a[NB_COEFF];

  	



  	free(tmp[0]);
  	free(tmp[1]);
  	free(tmp[2]);
  	//no need to free a*/

	printf("\nPERFORMANCE TESTS\n");

	int64_t pta[NB_COEFF];
	int64_t ptb[NB_COEFF];
	int64_t pc[NB_COEFF];

	__m256i avx_pta[NB_COEFF / 4];
	__m256i avx_ptb[NB_COEFF / 4];
	__m256i avx_pc[NB_COEFF / 4];

	mpz_t M;
	mpz_inits(M, NULL);
	mpz_set(M, rns_a.M);

	unsigned long long timer, meanTimer = 0, t1, t2;
	gmp_randstate_t state;
	gmp_randinit_default(state);

	/////////////////MUL/////////////////////////////////

	meanTimer = 0;

	// Heating caches
	for (int i = 0; i < NTEST; i++)
	{
		// appel de la fonction a mesurer a mettre ici
		// juste pour chauffer les caches
		mul_rns_cr(result_mul, &rns_a, rop5, rop6);
	}

	// timing
	for (int i = 0; i < NSAMPLES; i++)
	{

		// initialiser un nouveau jeu de donnees a tester
		mpz_urandomm(A, state, M); //Randomly generates A < M
		mpz_urandomm(B, state, M); //Randomly generates B < M
		from_int_to_rns(pta, &rns_a, A);
		from_int_to_rns(ptb, &rns_a, B);
		timer = (unsigned long long int)0x1 << 63;
		for (int j = 0; j < NTEST; j++)
		{
			t1 = cpucyclesStart();

			// appeler la fonction ici avec toujours le meme jeu de donnees
			mul_rns_cr(pc, &rns_a, pta, ptb);

			t2 = cpucyclesStop();

			if (timer > t2 - t1)
				timer = t2 - t1;
		}

		meanTimer += timer;
	}
	printf("\n RNS Standart multiplication : %lld CPU cycles  \n", meanTimer / NSAMPLES);

	meanTimer = 0;

	for (int i = 0; i < NTEST; i++)
	{
		// appel de la fonction a mesurer a mettre ici
		// juste pour chauffer les caches
		avx_mul_rns_cr(result_mul_avx_aux, &rns_a, rop5_aux, rop6_aux);
	}

	// timing
	for (int i = 0; i < NSAMPLES; i++)
	{

		// initialiser un nouveau jeu de donnees a tester
		mpz_urandomm(A, state, M); //Randomly generates A < M
		mpz_urandomm(B, state, M); //Randomly generates B < M
		from_int_to_rns(pta, &rns_a, A);
		from_int_to_rns(ptb, &rns_a, B);
		from_rns_to_m256i(avx_pta, &rns_a, pta);
		from_rns_to_m256i(avx_ptb, &rns_a, ptb);
		timer = (unsigned long long int)0x1 << 63;
		for (int j = 0; j < NTEST; j++)
		{
			t1 = cpucyclesStart();

			// appeler la fonction ici avec toujours le meme jeu de donnees
			avx_mul_rns_cr(avx_pc, &rns_a, avx_pta, avx_ptb);

			t2 = cpucyclesStop();

			if (timer > t2 - t1)
				timer = t2 - t1;
		}

		meanTimer += timer;
	}

	printf("\n RNS Standart multiplication vectored : %lld CPU cycles  \n", meanTimer / NSAMPLES);

	/////////////////ADD/////////////////////////////////

	meanTimer = 0;

	// Heating caches
	for (int i = 0; i < NTEST; i++)
	{
		// appel de la fonction a mesurer a mettre ici
		// juste pour chauffer les caches
		add_rns_cr(result, &rns_a, rop1, rop2);
	}

	// timing
	for (int i = 0; i < NSAMPLES; i++)
	{

		// initialiser un nouveau jeu de donnees a tester
		mpz_urandomm(A, state, M); //Randomly generates A < M
		mpz_urandomm(B, state, M); //Randomly generates B < M
		from_int_to_rns(pta, &rns_a, A);
		from_int_to_rns(ptb, &rns_a, B);
		timer = (unsigned long long int)0x1 << 63;
		for (int j = 0; j < NTEST; j++)
		{
			t1 = cpucyclesStart();

			// appeler la fonction ici avec toujours le meme jeu de donnees
			add_rns_cr(pc, &rns_a, pta, ptb);

			t2 = cpucyclesStop();

			if (timer > t2 - t1)
				timer = t2 - t1;
		}

		meanTimer += timer;
	}

	printf("\n RNS Addition: %lld CPU cycles  \n", meanTimer / NSAMPLES);

	meanTimer = 0;

	for (int i = 0; i < NTEST; i++)
	{
		// appel de la fonction a mesurer a mettre ici
		// juste pour chauffer les caches
		avx_add_rns_cr(result_avx_aux, &rns_a, rop1_aux, rop2_aux);
	}

	// timing
	for (int i = 0; i < NSAMPLES; i++)
	{

		// initialiser un nouveau jeu de donnees a tester
		mpz_urandomm(A, state, M); //Randomly generates A < M
		mpz_urandomm(B, state, M); //Randomly generates B < M
		from_int_to_rns(pta, &rns_a, A);
		from_int_to_rns(ptb, &rns_a, B);
		from_rns_to_m256i(avx_pta, &rns_a, pta);
		from_rns_to_m256i(avx_ptb, &rns_a, ptb);
		timer = (unsigned long long int)0x1 << 63;
		for (int j = 0; j < NTEST; j++)
		{
			t1 = cpucyclesStart();

			// appeler la fonction ici avec toujours le meme jeu de donnees
			avx_add_rns_cr(avx_pc, &rns_a, avx_pta, avx_ptb);

			t2 = cpucyclesStop();

			if (timer > t2 - t1)
				timer = t2 - t1;
		}

		meanTimer += timer;
	}

	printf("\n RNS Addition vectored: %lld CPU cycles  \n", meanTimer / NSAMPLES);

	return 0;
}