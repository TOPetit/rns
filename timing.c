#include <stdlib.h>
#include <stdio.h>
#include <gmp.h>

#include <math.h>

#include "rns.h"
#include "tests.c"

#include "rnsv.h"

int main(void){
	
	FILE* fpt;

	fpt = fopen("results/Results.json", "w+");

	fprintf(fpt, "{\n");

	printf("\nVectorized RNS timing :\n");

	printf("\n1. Generating data :\n");
	
	// Initializing random
	gmp_randstate_t state;
	gmp_randinit_default (state);
	// Timers
	unsigned long long timer, t1, t2;
	// Variables
	int64_t op1[NB_COEFF];
	int64_t op2[NB_COEFF];
	int64_t res[NB_COEFF];
	__m256i avx_op1[NB_COEFF/4];
	__m256i avx_op2[NB_COEFF/4];
	__m256i avx_res[NB_COEFF/4];

	mpz_t A, B, C;
	mpz_inits(A, B, C, NULL);

	// Base
	struct rns_base_t rns_a;
	rns_a.size = NB_COEFF;
	
	int64_t m_tmp[NB_COEFF] = {9223372036854775805,
9223372036854775801,
9223372036854775789,
9223372036854775783,
9223372036854775777,
9223372036854775769,
9223372036854775757,
9223372036854775747};
	rns_a.m = m_tmp;
	

	int k_tmp[NB_COEFF] = {3, 7, 19, 25, 31, 39, 51, 61};
	rns_a.k = k_tmp;

	init_rns(&rns_a);

	int64_t tmp_k[NB_COEFF];
    int j;
    for (j=0;j<NB_COEFF;j++){
    	tmp_k[j]=(int64_t)k_tmp[j];
    }

    __m256i avx_k[NB_COEFF/4];
    from_rns_to_m256i(avx_k,&rns_a,tmp_k);
    rns_a.avx_k = avx_k;

	
	mpz_t M;
	mpz_inits(M, NULL);
	mpz_set(M, rns_a.M); // Get M from the base
	unsigned long long int timing = ULLONG_MAX;
	unsigned long long int memory_cycles, memory_actual_cycles, memory_instructions, memory_ref;
	unsigned long before_cycles, after_cycles, cycles = ULONG_MAX;
	unsigned long before_instructions, after_instructions, instructions = ULONG_MAX;
	unsigned long before_ref, after_ref, ref = ULONG_MAX;

	
	printf("\n\tBase :\n");
	for (int i=0; i<NB_COEFF; i++) {printf("\t\t %ld\n", m_tmp[i]);}

	printf("\n2. Multiplication :\n");

	printf("\n\tHeating caches... ");
	mpz_urandomm(A, state, M); // Randomly generates A < M
	mpz_urandomm(B, state, M); // Randomly generated B < M
	from_int_to_rns(op1, &rns_a, A);
	from_int_to_rns(op2, &rns_a, B);

	for(int i=0; i<NTEST; i++){
		mul_rns_cr(res, &rns_a, op1, op2);
	}
	printf("Done.\n");
	
	printf("\tTesting... ");
	
	for(int i=0;i<NSAMPLES;i++) {
		mpz_urandomm (A, state, M);  //Randomly generates A < M
		mpz_urandomm (B, state, M);  //Randomly generates B < M
		from_int_to_rns(op1, &rns_a, A);
		from_int_to_rns(op2, &rns_a, B);
		unsigned long long int timer = 0;
		unsigned long int cycles_tmp = 0;
		unsigned long int instructions_tmp = 0;
		unsigned long int ref_tmp = 0;
		for(int j=0;j<NTEST;j++) {


			// RDTSC
			t1 = cpucyclesStart();

			mul_rns_cr(res, &rns_a, op1, op2);
			mul_rns_cr(res, &rns_a, op1, op2);
			mul_rns_cr(res, &rns_a, op1, op2);
			mul_rns_cr(res, &rns_a, op1, op2);
			mul_rns_cr(res, &rns_a, op1, op2);
			mul_rns_cr(res, &rns_a, op1, op2);
			mul_rns_cr(res, &rns_a, op1, op2);
			mul_rns_cr(res, &rns_a, op1, op2);
			mul_rns_cr(res, &rns_a, op1, op2);
			mul_rns_cr(res, &rns_a, op1, op2);

			t2 = cpucyclesStop();

			timer += (t2-t1)/10;

			// Instructions
			before_instructions = rdpmc_instructions();

			mul_rns_cr(res, &rns_a, op1, op2);
			mul_rns_cr(res, &rns_a, op1, op2);
			mul_rns_cr(res, &rns_a, op1, op2);
			mul_rns_cr(res, &rns_a, op1, op2);
			mul_rns_cr(res, &rns_a, op1, op2);
			mul_rns_cr(res, &rns_a, op1, op2);
			mul_rns_cr(res, &rns_a, op1, op2);
			mul_rns_cr(res, &rns_a, op1, op2);
			mul_rns_cr(res, &rns_a, op1, op2);
			mul_rns_cr(res, &rns_a, op1, op2);

			after_instructions = rdpmc_instructions();

			instructions_tmp += (after_instructions - before_instructions)/10;

			// actual cycles
			before_cycles = rdpmc_actual_cycles();

			mul_rns_cr(res, &rns_a, op1, op2);
			mul_rns_cr(res, &rns_a, op1, op2);
			mul_rns_cr(res, &rns_a, op1, op2);
			mul_rns_cr(res, &rns_a, op1, op2);
			mul_rns_cr(res, &rns_a, op1, op2);
			mul_rns_cr(res, &rns_a, op1, op2);
			mul_rns_cr(res, &rns_a, op1, op2);
			mul_rns_cr(res, &rns_a, op1, op2);
			mul_rns_cr(res, &rns_a, op1, op2);
			mul_rns_cr(res, &rns_a, op1, op2);

			after_cycles = rdpmc_actual_cycles();

			cycles_tmp += (after_cycles - before_cycles)/10;

			// reference cycles
			before_ref = rdpmc_reference_cycles();

			mul_rns_cr(res, &rns_a, op1, op2);
			mul_rns_cr(res, &rns_a, op1, op2);
			mul_rns_cr(res, &rns_a, op1, op2);
			mul_rns_cr(res, &rns_a, op1, op2);
			mul_rns_cr(res, &rns_a, op1, op2);
			mul_rns_cr(res, &rns_a, op1, op2);
			mul_rns_cr(res, &rns_a, op1, op2);
			mul_rns_cr(res, &rns_a, op1, op2);
			mul_rns_cr(res, &rns_a, op1, op2);
			mul_rns_cr(res, &rns_a, op1, op2);

			after_ref = rdpmc_reference_cycles();

			ref_tmp += (after_ref - before_ref)/10;
		}
		if (timing > timer/NTEST) timing = timer/NTEST;
		if (instructions > instructions_tmp/NTEST) instructions = instructions_tmp/NTEST;
		if (cycles > cycles_tmp/NTEST) cycles = cycles_tmp/NTEST;
		if (ref > ref_tmp/NTEST) ref = ref_tmp/NTEST;

	}

	fprintf(fpt, "\"multiplication\" :\n\t{\n");
	fprintf(fpt, "\t\"sequential\" :\n\t\t{\n");
	fprintf(fpt, "\t\t\"Cycles\" : %lld,\n\t\t\"Instructions\" : %ld,\n\t\t\"Actual cycles\" : %ld,\n\t\t\"Reference cycles\" : %ld\n\t\t},\n", timing, instructions, cycles, ref);

	printf("Done.\n");
	printf("\tRNS sequential multiplication : %lld CPU cycles.\n", timing);
	printf("\tRNS sequential multiplication : %ld instructions.\n", instructions);
	printf("\tRNS sequential multiplication : %ld actual CPU cycles.\n", cycles);
	printf("\tRNS sequential multiplication : %ld reference CPU cycles.\n", ref);
	
	memory_cycles = timing;
	memory_actual_cycles = cycles;
	memory_instructions = instructions;
	memory_ref = ref;

	timing = ULLONG_MAX;
	cycles = ULONG_MAX;
	instructions = ULONG_MAX;
	ref = ULONG_MAX;

	printf("\n\tHeating caches... ");
	mpz_urandomm (A, state, M);  //Randomly generates A < M
	mpz_urandomm (B, state, M);  //Randomly generates B < M
	
	from_int_to_rns(op1, &rns_a, A);
	from_int_to_rns(op2, &rns_a, B);
	from_rns_to_m256i(avx_op1,&rns_a,op1);
	from_rns_to_m256i(avx_op2,&rns_a,op2);
	
    for(int i=0;i<NTEST;i++)
	{
		
		avx_mul_rns_cr(avx_res,&rns_a,avx_op1, avx_op2);

	}
	
	printf("Done.\n");
	
	printf("\tTesting... ");
  
	// timing
	for(int i=0;i<NSAMPLES;i++)
	{
		mpz_urandomm (A, state, M);  //Randomly generates A < M
		mpz_urandomm (B, state, M);  //Randomly generates B < M
		from_int_to_rns(op1, &rns_a, A);
		from_int_to_rns(op2, &rns_a, B);
		from_rns_to_m256i(avx_op1,&rns_a,op1);
		from_rns_to_m256i(avx_op2,&rns_a,op2);
		unsigned long long int timer = 0;
		unsigned long int cycles_tmp = 0;
		unsigned long int instructions_tmp = 0;
		unsigned long int ref_tmp = 0;
		for(int j=0;j<NTEST;j++)
		{
			// RDTSC
			t1 = cpucyclesStart();

			avx_mul_rns_cr(avx_res, &rns_a, avx_op1, avx_op2);
			avx_mul_rns_cr(avx_res, &rns_a, avx_op1, avx_op2);
			avx_mul_rns_cr(avx_res, &rns_a, avx_op1, avx_op2);
			avx_mul_rns_cr(avx_res, &rns_a, avx_op1, avx_op2);
			avx_mul_rns_cr(avx_res, &rns_a, avx_op1, avx_op2);
			avx_mul_rns_cr(avx_res, &rns_a, avx_op1, avx_op2);
			avx_mul_rns_cr(avx_res, &rns_a, avx_op1, avx_op2);
			avx_mul_rns_cr(avx_res, &rns_a, avx_op1, avx_op2);
			avx_mul_rns_cr(avx_res, &rns_a, avx_op1, avx_op2);
			avx_mul_rns_cr(avx_res, &rns_a, avx_op1, avx_op2);
			
			t2 = cpucyclesStop();

			timer += (t2-t1)/10;

			// Instructions
			before_instructions = rdpmc_instructions();

			avx_mul_rns_cr(avx_res, &rns_a, avx_op1, avx_op2);
			avx_mul_rns_cr(avx_res, &rns_a, avx_op1, avx_op2);
			avx_mul_rns_cr(avx_res, &rns_a, avx_op1, avx_op2);
			avx_mul_rns_cr(avx_res, &rns_a, avx_op1, avx_op2);
			avx_mul_rns_cr(avx_res, &rns_a, avx_op1, avx_op2);
			avx_mul_rns_cr(avx_res, &rns_a, avx_op1, avx_op2);
			avx_mul_rns_cr(avx_res, &rns_a, avx_op1, avx_op2);
			avx_mul_rns_cr(avx_res, &rns_a, avx_op1, avx_op2);
			avx_mul_rns_cr(avx_res, &rns_a, avx_op1, avx_op2);
			avx_mul_rns_cr(avx_res, &rns_a, avx_op1, avx_op2);
			
			after_instructions = rdpmc_instructions();

			instructions_tmp += (after_instructions - before_instructions)/10;

			// actual cycles
			before_cycles = rdpmc_actual_cycles();

			avx_mul_rns_cr(avx_res, &rns_a, avx_op1, avx_op2);
			avx_mul_rns_cr(avx_res, &rns_a, avx_op1, avx_op2);
			avx_mul_rns_cr(avx_res, &rns_a, avx_op1, avx_op2);
			avx_mul_rns_cr(avx_res, &rns_a, avx_op1, avx_op2);
			avx_mul_rns_cr(avx_res, &rns_a, avx_op1, avx_op2);
			avx_mul_rns_cr(avx_res, &rns_a, avx_op1, avx_op2);
			avx_mul_rns_cr(avx_res, &rns_a, avx_op1, avx_op2);
			avx_mul_rns_cr(avx_res, &rns_a, avx_op1, avx_op2);
			avx_mul_rns_cr(avx_res, &rns_a, avx_op1, avx_op2);
			avx_mul_rns_cr(avx_res, &rns_a, avx_op1, avx_op2);
			
			after_cycles = rdpmc_actual_cycles();

			cycles_tmp += (after_cycles - before_cycles)/10;

			// reference cycles
			before_ref = rdpmc_reference_cycles();

			avx_mul_rns_cr(avx_res, &rns_a, avx_op1, avx_op2);
			avx_mul_rns_cr(avx_res, &rns_a, avx_op1, avx_op2);
			avx_mul_rns_cr(avx_res, &rns_a, avx_op1, avx_op2);
			avx_mul_rns_cr(avx_res, &rns_a, avx_op1, avx_op2);
			avx_mul_rns_cr(avx_res, &rns_a, avx_op1, avx_op2);
			avx_mul_rns_cr(avx_res, &rns_a, avx_op1, avx_op2);
			avx_mul_rns_cr(avx_res, &rns_a, avx_op1, avx_op2);
			avx_mul_rns_cr(avx_res, &rns_a, avx_op1, avx_op2);
			avx_mul_rns_cr(avx_res, &rns_a, avx_op1, avx_op2);
			avx_mul_rns_cr(avx_res, &rns_a, avx_op1, avx_op2);
			
			after_ref = rdpmc_reference_cycles();

			ref_tmp += (after_ref - before_ref)/10;
		}
		if (timing > timer/NTEST) timing = timer/NTEST;
		if (instructions > instructions_tmp/NTEST) instructions = instructions_tmp/NTEST;
		if (cycles > cycles_tmp/NTEST) cycles = cycles_tmp/NTEST;
		if (ref > ref_tmp/NTEST) ref = ref_tmp/NTEST;

	}

	fprintf(fpt, "\t\"vectorized\" :\n\t\t{\n");
	fprintf(fpt, "\t\t\"Cycles\" : %lld,\n\t\t\"Instructions\" : %ld,\n\t\t\"Actual cycles\" : %ld,\n\t\t\"Reference cycles\" : %ld\n\t\t}\n\t},\n", timing, instructions, cycles, ref);


	printf("Done.\n");
	printf("\tRNS vectorized multiplication : %lld CPU cycles.\n", timing);
	printf("\tRNS vectorized multiplication : %ld instructions.\n", instructions);
	printf("\tRNS vectorized multiplication : %ld actual CPU cycles.\n", cycles);
	printf("\tRNS vectorized multiplication : %ld reference CPU cycles.\n", ref);

	printf("\n\t -> %lld \%% cycles improvment.\n", 100 - (100 * timing) / memory_cycles);
	printf("\t -> %lld \%% instructions improvment.\n", 100 - (100 * instructions) / memory_instructions);
	printf("\t -> %lld \%% actual cycles improvment.\n", 100 - (100 * cycles) / memory_actual_cycles);
	printf("\t -> %lld \%% reference cycles improvment.\n", 100 - (100 * ref) / memory_ref);
	



	printf("\n\n2. Addition :\n");
	
	timing = ULLONG_MAX;
	cycles = ULONG_MAX;
	instructions = ULONG_MAX;
	ref = ULONG_MAX;

	printf("\n\tHeating caches... ");
	mpz_urandomm(A, state, M); // Randomly generates A < M
	mpz_urandomm(B, state, M); // Randomly generated B < M
	from_int_to_rns(op1, &rns_a, A);
	from_int_to_rns(op2, &rns_a, B);

	for(int i=0; i<NTEST; i++){
		add_rns_cr(res, &rns_a, op1, op2);
	}
	printf("Done.\n");
	
	printf("\tTesting... ");
	
	for(int i=0;i<NSAMPLES;i++) {
		mpz_urandomm (A, state, M);  //Randomly generates A < M
		mpz_urandomm (B, state, M);  //Randomly generates B < M
		from_int_to_rns(op1, &rns_a, A);
		from_int_to_rns(op2, &rns_a, B);
		unsigned long long int timer = 0;
		unsigned long int cycles_tmp = 0;
		unsigned long int instructions_tmp = 0;
		unsigned long int ref_tmp = 0;
		for(int j=0;j<NTEST;j++) {

			// RDTSC
			t1 = cpucyclesStart();

			add_rns_cr(res, &rns_a, op1, op2);
			add_rns_cr(res, &rns_a, op1, op2);
			add_rns_cr(res, &rns_a, op1, op2);
			add_rns_cr(res, &rns_a, op1, op2);
			add_rns_cr(res, &rns_a, op1, op2);
			add_rns_cr(res, &rns_a, op1, op2);
			add_rns_cr(res, &rns_a, op1, op2);
			add_rns_cr(res, &rns_a, op1, op2);
			add_rns_cr(res, &rns_a, op1, op2);
			add_rns_cr(res, &rns_a, op1, op2);
			
			t2 = cpucyclesStop();

			timer += (t2-t1)/10;

			// Instructions
			before_instructions = rdpmc_instructions();

			add_rns_cr(res, &rns_a, op1, op2);
			add_rns_cr(res, &rns_a, op1, op2);
			add_rns_cr(res, &rns_a, op1, op2);
			add_rns_cr(res, &rns_a, op1, op2);
			add_rns_cr(res, &rns_a, op1, op2);
			add_rns_cr(res, &rns_a, op1, op2);
			add_rns_cr(res, &rns_a, op1, op2);
			add_rns_cr(res, &rns_a, op1, op2);
			add_rns_cr(res, &rns_a, op1, op2);
			add_rns_cr(res, &rns_a, op1, op2);
			
			after_instructions = rdpmc_instructions();

			instructions_tmp += (after_instructions - before_instructions)/10;

			// actual cycles
			before_cycles = rdpmc_actual_cycles();

			add_rns_cr(res, &rns_a, op1, op2);
			add_rns_cr(res, &rns_a, op1, op2);
			add_rns_cr(res, &rns_a, op1, op2);
			add_rns_cr(res, &rns_a, op1, op2);
			add_rns_cr(res, &rns_a, op1, op2);
			add_rns_cr(res, &rns_a, op1, op2);
			add_rns_cr(res, &rns_a, op1, op2);
			add_rns_cr(res, &rns_a, op1, op2);
			add_rns_cr(res, &rns_a, op1, op2);
			add_rns_cr(res, &rns_a, op1, op2);
			
			after_cycles = rdpmc_actual_cycles();

			cycles_tmp += (after_cycles - before_cycles)/10;

			// reference cycles
			before_ref = rdpmc_reference_cycles();

			add_rns_cr(res, &rns_a, op1, op2);
			add_rns_cr(res, &rns_a, op1, op2);
			add_rns_cr(res, &rns_a, op1, op2);
			add_rns_cr(res, &rns_a, op1, op2);
			add_rns_cr(res, &rns_a, op1, op2);
			add_rns_cr(res, &rns_a, op1, op2);
			add_rns_cr(res, &rns_a, op1, op2);
			add_rns_cr(res, &rns_a, op1, op2);
			add_rns_cr(res, &rns_a, op1, op2);
			add_rns_cr(res, &rns_a, op1, op2);
			
			after_ref = rdpmc_reference_cycles();

			ref_tmp += (after_ref - before_ref)/10;
		}
		if (timing > timer/NTEST) timing = timer/NTEST;
		if (instructions > instructions_tmp/NTEST) instructions = instructions_tmp/NTEST;
		if (cycles > cycles_tmp/NTEST) cycles = cycles_tmp/NTEST;
		if (ref > ref_tmp/NTEST) ref = ref_tmp/NTEST;

	}

	fprintf(fpt, "\"addition\" :\n\t{\n");
	fprintf(fpt, "\t\"sequential\" :\n\t\t{\n");
	fprintf(fpt, "\t\t\"Cycles\" : %lld,\n\t\t\"Instructions\" : %ld,\n\t\t\"Actual cycles\" : %ld,\n\t\t\"Reference cycles\" : %ld\n\t\t},\n", timing, instructions, cycles, ref);


	printf("Done.\n");
	printf("\tRNS sequential addition : %lld CPU cycles.\n", timing);
	printf("\tRNS sequential addition : %ld instructions.\n", instructions);
	printf("\tRNS sequential addition : %ld actual CPU cycles.\n", cycles);
	printf("\tRNS sequential addition : %ld reference CPU cycles.\n", ref);

	memory_cycles = timing;
	memory_actual_cycles = cycles;
	memory_instructions = instructions;
	memory_ref = ref;

	timing = ULLONG_MAX;
	cycles = ULONG_MAX;
	instructions = ULONG_MAX;
	ref = ULONG_MAX;

	printf("\n\tHeating caches... ");
	mpz_urandomm (A, state, M);  //Randomly generates A < M
	mpz_urandomm (B, state, M);  //Randomly generates B < M
	
	from_int_to_rns(op1, &rns_a, A);
	from_int_to_rns(op2, &rns_a, B);
	from_rns_to_m256i(avx_op1,&rns_a,op1);
	from_rns_to_m256i(avx_op2,&rns_a,op2);
	
    for(int i=0;i<NTEST;i++)
	{
		
		avx_add_rns_cr(avx_res,&rns_a,avx_op1, avx_op2);

	}
	
	printf("Done.\n");
	
	printf("\tTesting... ");
  
	// timing
	for(int i=0;i<NSAMPLES;i++)
	{
		mpz_urandomm (A, state, M);  //Randomly generates A < M
		mpz_urandomm (B, state, M);  //Randomly generates B < M
		from_int_to_rns(op1, &rns_a, A);
		from_int_to_rns(op2, &rns_a, B);
		from_rns_to_m256i(avx_op1,&rns_a,op1);
		from_rns_to_m256i(avx_op2,&rns_a,op2);
		unsigned long long int timer = 0;
		unsigned long int cycles_tmp = 0;
		unsigned long int instructions_tmp = 0;
		unsigned long int ref_tmp = 0;
		for(int j=0;j<NTEST;j++)
		{
			// RDTSC
			t1 = cpucyclesStart();

			avx_add_rns_cr(avx_res, &rns_a, avx_op1, avx_op2);
			avx_add_rns_cr(avx_res, &rns_a, avx_op1, avx_op2);
			avx_add_rns_cr(avx_res, &rns_a, avx_op1, avx_op2);
			avx_add_rns_cr(avx_res, &rns_a, avx_op1, avx_op2);
			avx_add_rns_cr(avx_res, &rns_a, avx_op1, avx_op2);
			avx_add_rns_cr(avx_res, &rns_a, avx_op1, avx_op2);
			avx_add_rns_cr(avx_res, &rns_a, avx_op1, avx_op2);
			avx_add_rns_cr(avx_res, &rns_a, avx_op1, avx_op2);
			avx_add_rns_cr(avx_res, &rns_a, avx_op1, avx_op2);
			avx_add_rns_cr(avx_res, &rns_a, avx_op1, avx_op2);
			
			t2 = cpucyclesStop();

			timer += (t2-t1)/10;

			// Instructions
			before_instructions = rdpmc_instructions();

			avx_add_rns_cr(avx_res, &rns_a, avx_op1, avx_op2);
			avx_add_rns_cr(avx_res, &rns_a, avx_op1, avx_op2);
			avx_add_rns_cr(avx_res, &rns_a, avx_op1, avx_op2);
			avx_add_rns_cr(avx_res, &rns_a, avx_op1, avx_op2);
			avx_add_rns_cr(avx_res, &rns_a, avx_op1, avx_op2);
			avx_add_rns_cr(avx_res, &rns_a, avx_op1, avx_op2);
			avx_add_rns_cr(avx_res, &rns_a, avx_op1, avx_op2);
			avx_add_rns_cr(avx_res, &rns_a, avx_op1, avx_op2);
			avx_add_rns_cr(avx_res, &rns_a, avx_op1, avx_op2);
			avx_add_rns_cr(avx_res, &rns_a, avx_op1, avx_op2);
			
			after_instructions = rdpmc_instructions();

			instructions_tmp += (after_instructions - before_instructions)/10;

			// actual cycles
			before_cycles = rdpmc_actual_cycles();

			avx_add_rns_cr(avx_res, &rns_a, avx_op1, avx_op2);
			avx_add_rns_cr(avx_res, &rns_a, avx_op1, avx_op2);
			avx_add_rns_cr(avx_res, &rns_a, avx_op1, avx_op2);
			avx_add_rns_cr(avx_res, &rns_a, avx_op1, avx_op2);
			avx_add_rns_cr(avx_res, &rns_a, avx_op1, avx_op2);
			avx_add_rns_cr(avx_res, &rns_a, avx_op1, avx_op2);
			avx_add_rns_cr(avx_res, &rns_a, avx_op1, avx_op2);
			avx_add_rns_cr(avx_res, &rns_a, avx_op1, avx_op2);
			avx_add_rns_cr(avx_res, &rns_a, avx_op1, avx_op2);
			avx_add_rns_cr(avx_res, &rns_a, avx_op1, avx_op2);
			
			after_cycles = rdpmc_actual_cycles();

			cycles_tmp += (after_cycles - before_cycles)/10;

			// reference cycles
			before_ref = rdpmc_reference_cycles();

			avx_add_rns_cr(avx_res, &rns_a, avx_op1, avx_op2);
			avx_add_rns_cr(avx_res, &rns_a, avx_op1, avx_op2);
			avx_add_rns_cr(avx_res, &rns_a, avx_op1, avx_op2);
			avx_add_rns_cr(avx_res, &rns_a, avx_op1, avx_op2);
			avx_add_rns_cr(avx_res, &rns_a, avx_op1, avx_op2);
			avx_add_rns_cr(avx_res, &rns_a, avx_op1, avx_op2);
			avx_add_rns_cr(avx_res, &rns_a, avx_op1, avx_op2);
			avx_add_rns_cr(avx_res, &rns_a, avx_op1, avx_op2);
			avx_add_rns_cr(avx_res, &rns_a, avx_op1, avx_op2);
			avx_add_rns_cr(avx_res, &rns_a, avx_op1, avx_op2);
			
			after_ref = rdpmc_reference_cycles();

			ref_tmp += (after_ref - before_ref)/10;
		}
		if (timing > timer/NTEST) timing = timer/NTEST;
		if (instructions > instructions_tmp/NTEST) instructions = instructions_tmp/NTEST;
		if (cycles > cycles_tmp/NTEST) cycles = cycles_tmp/NTEST;
		if (ref > ref_tmp/NTEST) ref = ref_tmp/NTEST;

	}

	fprintf(fpt, "\t\"vectorized\" :\n\t\t{\n");
	fprintf(fpt, "\t\t\"Cycles\" : %lld,\n\t\t\"Instructions\" : %ld,\n\t\t\"Actual cycles\" : %ld,\n\t\t\"Reference cycles\" : %ld\n\t\t}\n\t},\n", timing, instructions, cycles, ref);


	printf("Done.\n");
	printf("\tRNS vectorized addition : %lld CPU cycles.\n", timing);
	printf("\tRNS vectorized addition : %ld instructions.\n", instructions);
	printf("\tRNS vectorized addition : %ld actual CPU cycles.\n", cycles);
	printf("\tRNS vectorized addition : %ld reference CPU cycles.\n", ref);

	printf("\n\t -> %lld \%% cycles improvment.\n", 100 - (100 * timing) / memory_cycles);
	printf("\t -> %lld \%% instructions improvment.\n", 100 - (100 * instructions) / memory_instructions);
	printf("\t -> %lld \%% actual cycles improvment.\n", 100 - (100 * cycles) / memory_actual_cycles);
	printf("\t -> %lld \%% reference cycles improvment.\n", 100 - (100 * ref) / memory_ref);


	printf("\n\n3. Substraction :\n");
	
	timing = ULLONG_MAX;
	cycles = ULONG_MAX;
	instructions = ULONG_MAX;
	ref = ULONG_MAX;

	printf("\n\tHeating caches... ");
	mpz_urandomm(A, state, M); // Randomly generates A < M
	mpz_urandomm(B, state, M); // Randomly generated B < M
	from_int_to_rns(op1, &rns_a, A);
	from_int_to_rns(op2, &rns_a, B);

	for(int i=0; i<NTEST; i++){
		sub_rns_cr(res, &rns_a, op1, op2);
	}
	printf("Done.\n");
	
	printf("\tTesting... ");
	
	for(int i=0;i<NSAMPLES;i++) {
		mpz_urandomm (A, state, M);  //Randomly generates A < M
		mpz_urandomm (B, state, M);  //Randomly generates B < M
		from_int_to_rns(op1, &rns_a, A);
		from_int_to_rns(op2, &rns_a, B);
		unsigned long long int timer = 0;
		unsigned long int cycles_tmp = 0;
		unsigned long int instructions_tmp = 0;
		unsigned long int ref_tmp = 0;
		for(int j=0;j<NTEST;j++) {

			// RDTSC
			t1 = cpucyclesStart();

			sub_rns_cr(res, &rns_a, op1, op2);
			sub_rns_cr(res, &rns_a, op1, op2);
			sub_rns_cr(res, &rns_a, op1, op2);
			sub_rns_cr(res, &rns_a, op1, op2);
			sub_rns_cr(res, &rns_a, op1, op2);
			sub_rns_cr(res, &rns_a, op1, op2);
			sub_rns_cr(res, &rns_a, op1, op2);
			sub_rns_cr(res, &rns_a, op1, op2);
			sub_rns_cr(res, &rns_a, op1, op2);
			sub_rns_cr(res, &rns_a, op1, op2);
			
			t2 = cpucyclesStop();

			timer += (t2-t1)/10;

			// Instructions
			before_instructions = rdpmc_instructions();

			sub_rns_cr(res, &rns_a, op1, op2);
			sub_rns_cr(res, &rns_a, op1, op2);
			sub_rns_cr(res, &rns_a, op1, op2);
			sub_rns_cr(res, &rns_a, op1, op2);
			sub_rns_cr(res, &rns_a, op1, op2);
			sub_rns_cr(res, &rns_a, op1, op2);
			sub_rns_cr(res, &rns_a, op1, op2);
			sub_rns_cr(res, &rns_a, op1, op2);
			sub_rns_cr(res, &rns_a, op1, op2);
			sub_rns_cr(res, &rns_a, op1, op2);
			
			after_instructions = rdpmc_instructions();

			instructions_tmp += (after_instructions - before_instructions)/10;

			// actual cycles
			before_cycles = rdpmc_actual_cycles();

			sub_rns_cr(res, &rns_a, op1, op2);
			sub_rns_cr(res, &rns_a, op1, op2);
			sub_rns_cr(res, &rns_a, op1, op2);
			sub_rns_cr(res, &rns_a, op1, op2);
			sub_rns_cr(res, &rns_a, op1, op2);
			sub_rns_cr(res, &rns_a, op1, op2);
			sub_rns_cr(res, &rns_a, op1, op2);
			sub_rns_cr(res, &rns_a, op1, op2);
			sub_rns_cr(res, &rns_a, op1, op2);
			sub_rns_cr(res, &rns_a, op1, op2);
			
			after_cycles = rdpmc_actual_cycles();

			cycles_tmp += (after_cycles - before_cycles)/10;

			// reference cycles
			before_ref = rdpmc_reference_cycles();

			sub_rns_cr(res, &rns_a, op1, op2);
			sub_rns_cr(res, &rns_a, op1, op2);
			sub_rns_cr(res, &rns_a, op1, op2);
			sub_rns_cr(res, &rns_a, op1, op2);
			sub_rns_cr(res, &rns_a, op1, op2);
			sub_rns_cr(res, &rns_a, op1, op2);
			sub_rns_cr(res, &rns_a, op1, op2);
			sub_rns_cr(res, &rns_a, op1, op2);
			sub_rns_cr(res, &rns_a, op1, op2);
			sub_rns_cr(res, &rns_a, op1, op2);
			
			after_ref = rdpmc_reference_cycles();

			ref_tmp += (after_ref - before_ref)/10;
		}
		if (timing > timer/NTEST) timing = timer/NTEST;
		if (instructions > instructions_tmp/NTEST) instructions = instructions_tmp/NTEST;
		if (cycles > cycles_tmp/NTEST) cycles = cycles_tmp/NTEST;
		if (ref > ref_tmp/NTEST) ref = ref_tmp/NTEST;

	}

	fprintf(fpt, "\"substraction\" :\n\t{\n");
	fprintf(fpt, "\t\"sequential\" :\n\t\t{\n");
	fprintf(fpt, "\t\t\"Cycles\" : %lld,\n\t\t\"Instructions\" : %ld,\n\t\t\"Actual cycles\" : %ld,\n\t\t\"Reference cycles\" : %ld\n\t\t},\n", timing, instructions, cycles, ref);

	printf("Done.\n");
	printf("\tRNS sequential substraction : %lld CPU cycles.\n", timing);
	printf("\tRNS sequential substraction : %ld instructions.\n", instructions);
	printf("\tRNS sequential substraction : %ld actual CPU cycles.\n", cycles);
	printf("\tRNS sequential substraction : %ld reference CPU cycles.\n", ref);

	memory_cycles = timing;
	memory_actual_cycles = cycles;
	memory_instructions = instructions;
	memory_ref = ref;

	timing = ULLONG_MAX;
	cycles = ULONG_MAX;
	instructions = ULONG_MAX;
	ref = ULONG_MAX;

	printf("\n\tHeating caches... ");
	mpz_urandomm (A, state, M);  //Randomly generates A < M
	mpz_urandomm (B, state, M);  //Randomly generates B < M
	
	from_int_to_rns(op1, &rns_a, A);
	from_int_to_rns(op2, &rns_a, B);
	from_rns_to_m256i(avx_op1,&rns_a,op1);
	from_rns_to_m256i(avx_op2,&rns_a,op2);
	
    for(int i=0;i<NTEST;i++)
	{
		
		avx_sub_rns_cr(avx_res,&rns_a,avx_op1, avx_op2);

	}
	
	printf("Done.\n");
	
	printf("\tTesting... ");
  
	// timing
	for(int i=0;i<NSAMPLES;i++)
	{
		mpz_urandomm (A, state, M);  //Randomly generates A < M
		mpz_urandomm (B, state, M);  //Randomly generates B < M
		from_int_to_rns(op1, &rns_a, A);
		from_int_to_rns(op2, &rns_a, B);
		from_rns_to_m256i(avx_op1,&rns_a,op1);
		from_rns_to_m256i(avx_op2,&rns_a,op2);
		unsigned long long int timer = 0;
		unsigned long int cycles_tmp = 0;
		unsigned long int instructions_tmp = 0;
		unsigned long int ref_tmp = 0;
		for(int j=0;j<NTEST;j++)
		{
			// RDTSC
			t1 = cpucyclesStart();

			avx_sub_rns_cr(avx_res, &rns_a, avx_op1, avx_op2);
			avx_sub_rns_cr(avx_res, &rns_a, avx_op1, avx_op2);
			avx_sub_rns_cr(avx_res, &rns_a, avx_op1, avx_op2);
			avx_sub_rns_cr(avx_res, &rns_a, avx_op1, avx_op2);
			avx_sub_rns_cr(avx_res, &rns_a, avx_op1, avx_op2);
			avx_sub_rns_cr(avx_res, &rns_a, avx_op1, avx_op2);
			avx_sub_rns_cr(avx_res, &rns_a, avx_op1, avx_op2);
			avx_sub_rns_cr(avx_res, &rns_a, avx_op1, avx_op2);
			avx_sub_rns_cr(avx_res, &rns_a, avx_op1, avx_op2);
			avx_sub_rns_cr(avx_res, &rns_a, avx_op1, avx_op2);
			
			t2 = cpucyclesStop();

			timer += (t2-t1)/10;

			// Instructions
			before_instructions = rdpmc_instructions();

			avx_sub_rns_cr(avx_res, &rns_a, avx_op1, avx_op2);
			avx_sub_rns_cr(avx_res, &rns_a, avx_op1, avx_op2);
			avx_sub_rns_cr(avx_res, &rns_a, avx_op1, avx_op2);
			avx_sub_rns_cr(avx_res, &rns_a, avx_op1, avx_op2);
			avx_sub_rns_cr(avx_res, &rns_a, avx_op1, avx_op2);
			avx_sub_rns_cr(avx_res, &rns_a, avx_op1, avx_op2);
			avx_sub_rns_cr(avx_res, &rns_a, avx_op1, avx_op2);
			avx_sub_rns_cr(avx_res, &rns_a, avx_op1, avx_op2);
			avx_sub_rns_cr(avx_res, &rns_a, avx_op1, avx_op2);
			avx_sub_rns_cr(avx_res, &rns_a, avx_op1, avx_op2);
			
			after_instructions = rdpmc_instructions();

			instructions_tmp += (after_instructions - before_instructions)/10;

			// actual cycles
			before_cycles = rdpmc_actual_cycles();

			avx_sub_rns_cr(avx_res, &rns_a, avx_op1, avx_op2);
			avx_sub_rns_cr(avx_res, &rns_a, avx_op1, avx_op2);
			avx_sub_rns_cr(avx_res, &rns_a, avx_op1, avx_op2);
			avx_sub_rns_cr(avx_res, &rns_a, avx_op1, avx_op2);
			avx_sub_rns_cr(avx_res, &rns_a, avx_op1, avx_op2);
			avx_sub_rns_cr(avx_res, &rns_a, avx_op1, avx_op2);
			avx_sub_rns_cr(avx_res, &rns_a, avx_op1, avx_op2);
			avx_sub_rns_cr(avx_res, &rns_a, avx_op1, avx_op2);
			avx_sub_rns_cr(avx_res, &rns_a, avx_op1, avx_op2);
			avx_sub_rns_cr(avx_res, &rns_a, avx_op1, avx_op2);
			
			after_cycles = rdpmc_actual_cycles();

			cycles_tmp += (after_cycles - before_cycles)/10;

			// reference cycles
			before_ref = rdpmc_reference_cycles();

			avx_sub_rns_cr(avx_res, &rns_a, avx_op1, avx_op2);
			avx_sub_rns_cr(avx_res, &rns_a, avx_op1, avx_op2);
			avx_sub_rns_cr(avx_res, &rns_a, avx_op1, avx_op2);
			avx_sub_rns_cr(avx_res, &rns_a, avx_op1, avx_op2);
			avx_sub_rns_cr(avx_res, &rns_a, avx_op1, avx_op2);
			avx_sub_rns_cr(avx_res, &rns_a, avx_op1, avx_op2);
			avx_sub_rns_cr(avx_res, &rns_a, avx_op1, avx_op2);
			avx_sub_rns_cr(avx_res, &rns_a, avx_op1, avx_op2);
			avx_sub_rns_cr(avx_res, &rns_a, avx_op1, avx_op2);
			avx_sub_rns_cr(avx_res, &rns_a, avx_op1, avx_op2);
			
			after_ref = rdpmc_reference_cycles();

			ref_tmp += (after_ref - before_ref)/10;
		}
		if (timing > timer/NTEST) timing = timer/NTEST;
		if (instructions > instructions_tmp/NTEST) instructions = instructions_tmp/NTEST;
		if (cycles > cycles_tmp/NTEST) cycles = cycles_tmp/NTEST;
		if (ref > ref_tmp/NTEST) ref = ref_tmp/NTEST;

	}

	fprintf(fpt, "\t\"vectorized\" :\n\t\t{\n");
	fprintf(fpt, "\t\t\"Cycles\" : %lld,\n\t\t\"Instructions\" : %ld,\n\t\t\"Actual cycles\" : %ld,\n\t\t\"Reference cycles\" : %ld\n\t\t}\n\t},\n", timing, instructions, cycles, ref);


	printf("Done.\n");
	printf("\tRNS vectorized substraction : %lld CPU cycles.\n", timing);
	printf("\tRNS vectorized substraction : %ld instructions.\n", instructions);
	printf("\tRNS vectorized substraction : %ld actual CPU cycles.\n", cycles);
	printf("\tRNS vectorized substraction : %ld reference CPU cycles.\n", ref);

	printf("\n\t -> %lld \%% cycles improvment.\n", 100 - (100 * timing) / memory_cycles);
	printf("\t -> %lld \%% instructions improvment.\n", 100 - (100 * instructions) / memory_instructions);
	printf("\t -> %lld \%% actual cycles improvment.\n", 100 - (100 * cycles) / memory_actual_cycles);
	printf("\t -> %lld \%% reference cycles improvment.\n", 100 - (100 * ref) / memory_ref);

	


	// Second Base
	struct rns_base_t rns_b;

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

    int64_t tmp_k2[NB_COEFF];
    for (j=0;j<NB_COEFF;j++){
    	tmp_k2[j]=(int64_t)k2_bis[j];
    }
    __m256i avx_k2[NB_COEFF/4];
    from_rns_to_m256i(avx_k2,&rns_b,tmp_k2); 
    rns_b.avx_k = avx_k2;


	// Base conversion
	struct conv_base_t conv;
	conv.rns_a = &rns_a;
	conv.rns_b = &rns_b;
	initialize_inverses_base_conversion(&conv);


	int64_t ttt[NB_COEFF];

	timing = ULLONG_MAX;
	cycles = ULONG_MAX;
	instructions = ULONG_MAX;
	ref = ULONG_MAX;

	printf("\n\n4. Base conversion :\n");

	printf("\n\tHeating caches... ");
	mpz_urandomm(A, state, M); // Randomly generates A < M
	mpz_urandomm(B, state, M); // Randomly generated B < M
	from_int_to_rns(op1, &rns_a, A);
	from_int_to_rns(op2, &rns_a, B);

	for(int i=0; i<NTEST; i++){
		sub_rns_cr(res, &rns_a, op1, op2);
	}
	printf("Done.\n");
	
	printf("\tTesting... ");
	
	for(int i=0;i<NSAMPLES;i++) {
		mpz_urandomm (A, state, M);  //Randomly generates A < M
		mpz_urandomm (B, state, M);  //Randomly generates B < M
		from_int_to_rns(op1, &rns_a, A);
		unsigned long long int timer = 0;
		unsigned long int cycles_tmp = 0;
		unsigned long int instructions_tmp = 0;
		unsigned long int ref_tmp = 0;
		for(int j=0;j<NTEST;j++) {

			// RDTSC
			t1 = cpucyclesStart();

			base_conversion_cr(op2, &conv, op1, ttt);
			base_conversion_cr(op2, &conv, op1, ttt);
			base_conversion_cr(op2, &conv, op1, ttt);
			base_conversion_cr(op2, &conv, op1, ttt);
			base_conversion_cr(op2, &conv, op1, ttt);
			base_conversion_cr(op2, &conv, op1, ttt);
			base_conversion_cr(op2, &conv, op1, ttt);
			base_conversion_cr(op2, &conv, op1, ttt);
			base_conversion_cr(op2, &conv, op1, ttt);
			base_conversion_cr(op2, &conv, op1, ttt);
			
			t2 = cpucyclesStop();

			timer += (t2-t1)/10;

			// Instructions
			before_instructions = rdpmc_instructions();

			base_conversion_cr(op2, &conv, op1, ttt);
			base_conversion_cr(op2, &conv, op1, ttt);
			base_conversion_cr(op2, &conv, op1, ttt);
			base_conversion_cr(op2, &conv, op1, ttt);
			base_conversion_cr(op2, &conv, op1, ttt);
			base_conversion_cr(op2, &conv, op1, ttt);
			base_conversion_cr(op2, &conv, op1, ttt);
			base_conversion_cr(op2, &conv, op1, ttt);
			base_conversion_cr(op2, &conv, op1, ttt);
			base_conversion_cr(op2, &conv, op1, ttt);
			
			after_instructions = rdpmc_instructions();

			instructions_tmp += (after_instructions - before_instructions)/10;
			
			// actual cycles
			before_cycles = rdpmc_actual_cycles();

			base_conversion_cr(op2, &conv, op1, ttt);
			base_conversion_cr(op2, &conv, op1, ttt);
			base_conversion_cr(op2, &conv, op1, ttt);
			base_conversion_cr(op2, &conv, op1, ttt);
			base_conversion_cr(op2, &conv, op1, ttt);
			base_conversion_cr(op2, &conv, op1, ttt);
			base_conversion_cr(op2, &conv, op1, ttt);
			base_conversion_cr(op2, &conv, op1, ttt);
			base_conversion_cr(op2, &conv, op1, ttt);
			base_conversion_cr(op2, &conv, op1, ttt);
			
			after_cycles = rdpmc_actual_cycles();

			cycles_tmp += (after_cycles - before_cycles)/10;

			// reference cycles
			before_ref = rdpmc_reference_cycles();

			base_conversion_cr(op2, &conv, op1, ttt);
			base_conversion_cr(op2, &conv, op1, ttt);
			base_conversion_cr(op2, &conv, op1, ttt);
			base_conversion_cr(op2, &conv, op1, ttt);
			base_conversion_cr(op2, &conv, op1, ttt);
			base_conversion_cr(op2, &conv, op1, ttt);
			base_conversion_cr(op2, &conv, op1, ttt);
			base_conversion_cr(op2, &conv, op1, ttt);
			base_conversion_cr(op2, &conv, op1, ttt);
			base_conversion_cr(op2, &conv, op1, ttt);
			
			after_ref = rdpmc_reference_cycles();

			ref_tmp += (after_ref - before_ref)/10;
			;
		}
		if (timing > timer/NTEST) timing = timer/NTEST;
		if (instructions > instructions_tmp/NTEST) instructions = instructions_tmp/NTEST;
		if (cycles > cycles_tmp/NTEST) cycles = cycles_tmp/NTEST;
		if (ref > ref_tmp/NTEST) ref = ref_tmp/NTEST;

	}

	fprintf(fpt, "\"base_conversion\" :\n\t{\n");
	fprintf(fpt, "\t\"sequential\" :\n\t\t{\n");
	fprintf(fpt, "\t\t\"Cycles\" : %lld,\n\t\t\"Instructions\" : %ld,\n\t\t\"Actual cycles\" : %ld,\n\t\t\"Reference cycles\" : %ld\n\t\t},\n", timing, instructions, cycles, ref);


	printf("Done.\n");
	printf("\tRNS sequential base conversion : %lld CPU cycles.\n", timing);
	printf("\tRNS sequential base conversion : %ld instructions.\n", instructions);
	printf("\tRNS sequential base conversion : %ld actual CPU cycles.\n", cycles);
	printf("\tRNS sequential base conversion : %ld reference CPU cycles.\n", ref);


	memory_cycles = timing;
	memory_actual_cycles = cycles;
	memory_instructions = instructions;
	memory_ref = ref;

	timing = ULLONG_MAX;
	cycles = ULONG_MAX;
	instructions = ULONG_MAX;
	ref = ULONG_MAX;

	printf("\n\tHeating caches... ");
	mpz_urandomm (A, state, M);  //Randomly generates A < M
	
	from_int_to_rns(op1, &rns_a, A);
	from_rns_to_m256i(avx_op1,&rns_a,op1);
	
    for(int i=0;i<NTEST;i++)
	{
		
		avx_base_conversion_cr(avx_op2, &conv, avx_op1, ttt);

	}
	
	printf("Done.\n");
	
	printf("\tTesting... ");
  
	// timing
	for(int i=0;i<NSAMPLES;i++)
	{
		mpz_urandomm (A, state, M);  //Randomly generates A < M
		from_int_to_rns(op1, &rns_a, A);
		from_rns_to_m256i(avx_op1,&rns_a,op1);
		unsigned long long int timer = 0;
		unsigned long int cycles_tmp = 0;
		unsigned long int instructions_tmp = 0;
		unsigned long int ref_tmp = 0;
		for(int j=0;j<NTEST;j++)
		{
			// RDTSC
			t1 = cpucyclesStart();

			avx_base_conversion_cr(avx_op2, &conv, avx_op1, ttt);
			avx_base_conversion_cr(avx_op2, &conv, avx_op1, ttt);
			avx_base_conversion_cr(avx_op2, &conv, avx_op1, ttt);
			avx_base_conversion_cr(avx_op2, &conv, avx_op1, ttt);
			avx_base_conversion_cr(avx_op2, &conv, avx_op1, ttt);
			avx_base_conversion_cr(avx_op2, &conv, avx_op1, ttt);
			avx_base_conversion_cr(avx_op2, &conv, avx_op1, ttt);
			avx_base_conversion_cr(avx_op2, &conv, avx_op1, ttt);
			avx_base_conversion_cr(avx_op2, &conv, avx_op1, ttt);
			avx_base_conversion_cr(avx_op2, &conv, avx_op1, ttt);
			
			t2 = cpucyclesStop();

			timer += (t2-t1)/10;

			// Instructions
			before_instructions = rdpmc_instructions();

			avx_base_conversion_cr(avx_op2, &conv, avx_op1, ttt);
			avx_base_conversion_cr(avx_op2, &conv, avx_op1, ttt);
			avx_base_conversion_cr(avx_op2, &conv, avx_op1, ttt);
			avx_base_conversion_cr(avx_op2, &conv, avx_op1, ttt);
			avx_base_conversion_cr(avx_op2, &conv, avx_op1, ttt);
			avx_base_conversion_cr(avx_op2, &conv, avx_op1, ttt);
			avx_base_conversion_cr(avx_op2, &conv, avx_op1, ttt);
			avx_base_conversion_cr(avx_op2, &conv, avx_op1, ttt);
			avx_base_conversion_cr(avx_op2, &conv, avx_op1, ttt);
			avx_base_conversion_cr(avx_op2, &conv, avx_op1, ttt);
			
			after_instructions = rdpmc_instructions();

			instructions_tmp += (after_instructions - before_instructions)/10;

			// actual cycles
			before_cycles = rdpmc_actual_cycles();

			avx_base_conversion_cr(avx_op2, &conv, avx_op1, ttt);
			avx_base_conversion_cr(avx_op2, &conv, avx_op1, ttt);
			avx_base_conversion_cr(avx_op2, &conv, avx_op1, ttt);
			avx_base_conversion_cr(avx_op2, &conv, avx_op1, ttt);
			avx_base_conversion_cr(avx_op2, &conv, avx_op1, ttt);
			avx_base_conversion_cr(avx_op2, &conv, avx_op1, ttt);
			avx_base_conversion_cr(avx_op2, &conv, avx_op1, ttt);
			avx_base_conversion_cr(avx_op2, &conv, avx_op1, ttt);
			avx_base_conversion_cr(avx_op2, &conv, avx_op1, ttt);
			avx_base_conversion_cr(avx_op2, &conv, avx_op1, ttt);
			
			after_cycles = rdpmc_actual_cycles();

			cycles_tmp += (after_cycles - before_cycles)/10;

			// reference cycles
			before_ref = rdpmc_reference_cycles();

			avx_base_conversion_cr(avx_op2, &conv, avx_op1, ttt);
			avx_base_conversion_cr(avx_op2, &conv, avx_op1, ttt);
			avx_base_conversion_cr(avx_op2, &conv, avx_op1, ttt);
			avx_base_conversion_cr(avx_op2, &conv, avx_op1, ttt);
			avx_base_conversion_cr(avx_op2, &conv, avx_op1, ttt);
			avx_base_conversion_cr(avx_op2, &conv, avx_op1, ttt);
			avx_base_conversion_cr(avx_op2, &conv, avx_op1, ttt);
			avx_base_conversion_cr(avx_op2, &conv, avx_op1, ttt);
			avx_base_conversion_cr(avx_op2, &conv, avx_op1, ttt);
			avx_base_conversion_cr(avx_op2, &conv, avx_op1, ttt);
			
			after_ref = rdpmc_reference_cycles();

			ref_tmp += (after_ref - before_ref)/10;
		}
		if (timing > timer/NTEST) timing = timer/NTEST;
		if (instructions > instructions_tmp/NTEST) instructions = instructions_tmp/NTEST;
		if (cycles > cycles_tmp/NTEST) cycles = cycles_tmp/NTEST;
		if (ref > ref_tmp/NTEST) ref = ref_tmp/NTEST;

	}

	fprintf(fpt, "\t\"vectorized\" :\n\t\t{\n");
	fprintf(fpt, "\t\t\"Cycles\" : %lld,\n\t\t\"Instructions\" : %ld,\n\t\t\"Actual cycles\" : %ld,\n\t\t\"Reference cycles\" : %ld\n\t\t}\n\t}\n}", timing, instructions, cycles, ref);


	printf("Done.\n");
	printf("\tRNS vectorized base conversion : %lld CPU cycles.\n", timing);
	printf("\tRNS vectorized base conversion : %ld instructions.\n", instructions);
	printf("\tRNS vectorized base conversion : %ld actual CPU cycles.\n", cycles);
	printf("\tRNS vectorized base conversion : %ld reference CPU cycles.\n", ref);

	printf("\n\t -> %lld \%% cycles improvment.\n", 100 - (100 * timing) / memory_cycles);
	printf("\t -> %lld \%% instructions improvment.\n", 100 - (100 * instructions) / memory_instructions);
	printf("\t -> %lld \%% actual cycles improvment.\n", 100 - (100 * cycles) / memory_actual_cycles);
	printf("\t -> %lld \%% reference cycles improvment.\n", 100 - (100 * ref) / memory_ref);

	fclose(fpt);

	return 0;
}

