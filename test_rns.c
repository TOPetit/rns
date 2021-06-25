/////////////////////////////////////////////////////////////////////
///
/// ./rns P nb_base1 m1 m2 .... nb_base2 m'A m'2 ....
///
/////////////////////////////////////////////////////////////////////
#include <stdlib.h>

#include <stdio.h>

#include <gmp.h>

#include <stdint.h>

#include "rns.h"

#include "tests.c"
 //~ The module p
mpz_t modul_p;

#define N 63

typedef __int128 int128;

int main(int argc, char ** argv) {
  int i;
  int64_t * pa;
  int64_t * pb;
  int64_t * pab;
  int64_t * pbb;
  int64_t * pc;
  int64_t * pp1;
  int64_t * pp2;
  int64_t * pp3;
  int128 temp;

  mpz_t modul_p;

  mpz_t A, B, C, E, F, inv_p_modM, inv_M_modMp;
  mpz_t Arns, Brns;
  mpz_t g, s, t;
  mpz_init(modul_p);
  mpz_inits(A, B, C, E, F, inv_p_modM, inv_M_modMp, NULL);
  mpz_inits(Arns, Brns, NULL);
  mpz_inits(g, s, t, NULL);

  //printf("Initializations done \n");
  //printf("first parameter %s \n", argv[1]);

  mpz_set_str(modul_p, argv[1], 10); // Load P

  gmp_printf("%Zd ", modul_p);

  //printf("First load done\n");

  // Load RNS base
  struct rns_base_t rns_a;
  int64_t * base;
  int * k;
  int size_base = atoi(argv[2]);
  base = (int64_t * ) malloc(size_base * sizeof(int64_t));
  k = (int * ) malloc(size_base * sizeof(int));

  printf("%d ", size_base);

  for (i = 3; i < size_base + 3; i++) {

    // printf("%s \n",argv[i]);

    k[i - 3] = atoi(argv[i]);
    // printf("%i \n", k[i-3]);

    temp = ((int128) 1 << N) - k[i - 3];
    base[i - 3] = (int64_t) temp;
  }
  rns_a.size = size_base;
  rns_a.m = base;
  rns_a.k = k;

  //printf("before RNS \n");

  init_rns( & rns_a);

  //printf("first base done\n");

  //Compute P^{-1} mod M
  // mpz_set_str (inv_p_modM, argv[size_base+3], 10); 

  mpz_gcdext(g, inv_p_modM, t, modul_p, rns_a.M);

  // gmp_printf("M   : %Zd\n", rns_a.M);
  // for(i=0; i<rns_a.size;i++){
  // 	printf("%lu \n", rns_a.m[i]);
  // 	printf("%i \n", rns_a.k[i]);
  // }
  // gmp_printf("inv_p_modM   : %Zd\n", inv_p_modM);
  // gmp_printf("inv_p_modM   : %Zd\n", s);
  // gmp_printf("inv_p_modM   : %Zd\n", t);

  // Load RNS base2
  struct rns_base_t rns_b;
  int64_t * base2;
  int * k2;
  int size_base2 = atoi(argv[size_base + 3]);
  base2 = (int64_t * ) malloc(size_base2 * sizeof(int64_t));
  k2 = (int * ) malloc(size_base2 * sizeof(int));

  printf("%d ", size_base2);

  for (i = size_base + 4; i < size_base2 + size_base + 4; i++) {
    k2[i - size_base - 4] = atoi(argv[i]);
    temp = ((int128) 1 << N) - k2[i - size_base - 4];
    base2[i - size_base - 4] = (int64_t) temp;
  }
  rns_b.size = size_base2;
  rns_b.m = base2;
  rns_b.k = k2;

  init_rns( & rns_b);

  //Compute M^{-1} mod M2
  //	mpz_set_str (inv_M_modMp, argv[size_base2+size_base+4], 10); 

  mpz_gcdext(g, s, t, rns_a.M, rns_b.M);

  // Initialize conversion constants
  int max_size = size_base > size_base2 ? size_base : size_base2;
  pa = (int64_t * ) malloc(max_size * sizeof(int64_t));
  pb = (int64_t * ) malloc(max_size * sizeof(int64_t));
  pab = (int64_t * ) malloc(max_size * sizeof(int64_t));
  pbb = (int64_t * ) malloc(max_size * sizeof(int64_t));
  pc = (int64_t * ) malloc(max_size * sizeof(int64_t));
  pp1 = (int64_t * ) malloc(max_size * sizeof(int64_t));
  pp2 = (int64_t * ) malloc(max_size * sizeof(int64_t));
  pp3 = (int64_t * ) malloc(max_size * sizeof(int64_t));

  struct conv_base_t conv;

  conv.rns_a = & rns_a;
  conv.rns_b = & rns_b;

  initialize_inverses_base_conversion( & conv);

  // Initialize multiplication constants
  struct mod_mul_t mult;
  from_int_to_rns(pa, & rns_a, A);
  from_int_to_rns(pb, & rns_a, B);
  from_int_to_rns(pab, & rns_b, A);
  from_int_to_rns(pb, & rns_b, B);
  from_int_to_rns(pp1, & rns_a, inv_p_modM); //(-P)^-1 mod Ma
  from_int_to_rns(pp2, & rns_b, modul_p); // P mod Mb
  from_int_to_rns(pp3, & rns_b, inv_M_modMp); // Ma^{-1} mod Mb

  mult.inv_p_modMa = pp1;
  mult.p_modMb = pp2;
  mult.inv_Ma_modMb = pp3;
  mult.conv = & conv;

  // Initialize temporary sapce for RNs computations
  int64_t * tmp[4]; // RNS modular multiplication intermediate results
  // One more for the base convertion
  tmp[0] = (int64_t * ) malloc(max_size * sizeof(int64_t)); //////////////////////////////
  tmp[1] = (int64_t * ) malloc(max_size * sizeof(int64_t)); //////////////////////////////
  tmp[2] = (int64_t * ) malloc(max_size * sizeof(int64_t)); //////////////////////////////
  tmp[3] = (int64_t * ) malloc(max_size * sizeof(int64_t)); //////////////////////////////

  // Timing tests
  gmp_randstate_t state;
  gmp_randinit_default(state);
  unsigned long long timer, meanTimer = 0, t1, t2;

  meanTimer = 0;
  // Heating caches
  for (int i = 0; i < NTEST; i++) {
    // appel de la fonction a mesurer a mettre ici
    // juste pour chauffer les caches
    mult_mod_rns(pc, pa, pab, pb, pbb, & mult, tmp);

  }

  // timing
  for (int i = 0; i < NSAMPLES; i++) {

    // initialiser un nouveau jeu de donnees a tester
    mpz_urandomm(A, state, modul_p); //Randomly generates A < P
    mpz_urandomm(B, state, modul_p); //Randomly generates B < P
    from_int_to_rns(pa, & rns_a, A);
    from_int_to_rns(pb, & rns_a, B);
    from_int_to_rns(pab, & rns_b, A);
    from_int_to_rns(pbb, & rns_b, B);
    timer = (unsigned long long int) 0x1 << 63;
    for (int j = 0; j < NTEST; j++) {
      t1 = cpucyclesStart();
      // appeler la fonction ici avec toujours le meme jeu de donnees
      mult_mod_rns(pc, pa, pab, pb, pbb, & mult, tmp);

      t2 = cpucyclesStop();
      if (timer > t2 - t1) timer = t2 - t1;
    }

    meanTimer += timer;
  }

  printf("%lld ", meanTimer / NSAMPLES);

  meanTimer = 0;
  // Heating caches
  for (int i = 0; i < NTEST; i++) {
    // appel de la fonction a mesurer a mettre ici
    // juste pour chauffer les caches
    mult_mod_rns_cr(pc, pa, pab, pb, pbb, & mult, tmp);

  }

  // timing
  for (int i = 0; i < NSAMPLES; i++) {

    // initialiser un nouveau jeu de donnees a tester
    mpz_urandomm(A, state, modul_p); //Randomly generates A < P
    mpz_urandomm(B, state, modul_p); //Randomly generates B < P
    from_int_to_rns(pa, & rns_a, A);
    from_int_to_rns(pb, & rns_a, B);
    timer = (unsigned long long int) 0x1 << 63;
    for (int j = 0; j < NTEST; j++) {
      t1 = cpucyclesStart();
      // appeler la fonction ici avec toujours le meme jeu de donnees
      mult_mod_rns_cr(pc, pa, pab, pb, pbb, & mult, tmp);

      t2 = cpucyclesStop();
      if (timer > t2 - t1) timer = t2 - t1;
    }

    meanTimer += timer;
  }

  printf("%lld \n", meanTimer / NSAMPLES);

}