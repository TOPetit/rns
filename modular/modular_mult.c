#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include "rns.h"
#include "tests.c"

#define W 32
#define WORD_MAX UINT32_MAX
typedef uint32_t word;
typedef uint64_t double_word;

word gcdExtended(word a, word b, word *x, word *y);
static inline word montgommery(word a, word b, word p, word r);
static inline word barrett(word a, word b, word p, word r);
static inline word mersenne_32(word a, word b);
static inline word pseudo_mersenne(word a, word b, word k, word p);
static inline word generalized_mersenne(word a, word b, word k, word p);
static inline word montgommery_friendly(word a,
                                        word b,
                                        word k,
                                        word e,
                                        word p,
                                        double_word mask);
static inline word plantard(word a, word b, word p, double_word r);

word gcdExtended(word a, word b, word *x, word *y)
{
  // Base Case
  if (a == 0)
  {
    *x = 0;
    *y = 1;
    return b;
  }

  word x1, y1; // To store results of recursive call
  word gcd = gcdExtended(b % a, a, &x1, &y1);

  // Update x and y using results of recursive
  // call
  *x = y1 - (b / a) * x1;
  *y = x1;

  return gcd;
}

//////////////////////////////////////
// Montgommery word multiplication
//////////////////////////////////////
// p : the modulus
// r : (-p^{-1} mod 2^n)
// output : a*b*2^{-n} mod p
//////////////////////////////////////
static inline word montgommery(word a, word b, word p, word r)
{
  double_word c;
  //	double_word tmp1_dw, tmp2_dw;
  word tmp1_w, tmp2_w, tmp3_w;

  c = (double_word)a * b;
  tmp1_w = (word)(c * r);
  tmp2_w = (c + tmp1_w * p) >> W;
  if (tmp2_w > p)
    return tmp2_w - p;
  return tmp2_w;
}

//////////////////////////////////////
// Barrett word multiplication
//////////////////////////////////////
// p : the modulus
// r : ((double_word)1<<W)/p  ie :int(2^n/p)
// output : a*b mod p
//////////////////////////////////////
static inline word barrett(word a, word b, word p, word r)
{
  double_word c;
  word tmp1_w, tmp2_w, tmp3_w;
  double_word tmp1_dw, tmp2_dw;

  c = (double_word)a * b;
  tmp1_dw = ((c >> (W)) * r);
  tmp2_w = c - tmp1_dw * p;

  if (tmp2_w >= 2 * WORD_MAX)
    return tmp2_w - 2 * p;
  if (tmp2_w >= WORD_MAX)
    return tmp2_w - p;
  return tmp2_w;
}

//////////////////////////////////////
// Mersenne word multiplication
//////////////////////////////////////
// p is implicitly 2^n-1
// output : a*b mod p
//////////////////////////////////////
static inline word mersenne_32(word a, word b)
{
  word tmp1_w;
  double_word c;

  c = (double_word)a * b;
  tmp1_w = (word)c + (c >> W);
  if (tmp1_w >= WORD_MAX) // Inultile pour WORD_MAX ???
    return (word)(c - WORD_MAX);
  else
    return (word)c;
}

//////////////////////////////////////
// Pseudo Mersenne word multiplication
//////////////////////////////////////
// p : the modulus p = 2^n-k
// k small enough
// output : a*b mod p
//////////////////////////////////////
static inline word pseudo_mersenne(word a, word b, word k, word p)
{
  double_word c;
  double_word tmp1_dw;
  word tmp2_w;

  c = (double_word)a * b;
  tmp1_dw = (c >> W) * k + (word)c;
  //	printf(" pseudo %lu ", tmp1_dw);
  tmp2_w = (tmp1_dw >> W) * k + (word)tmp1_dw;
  if (tmp2_w >= p)
    return tmp2_w - p;
  return tmp2_w;
}

//////////////////////////////////////
// Generalized Mersenne word multiplication
//////////////////////////////////////
// p : the modulus p = 2^n-2^k-1
// k small enough
// output : a*b mod p
//////////////////////////////////////
static inline word generalized_mersenne(word a, word b, word k, word p)
{
  double_word c;
  double_word tmp1_dw;
  word tmp2_w;

  c = (double_word)a * b;
  tmp1_dw = (c >> W) + ((c >> W) << k) + (word)c;
  //	printf(" generalized %lu ", tmp1_dw);
  tmp2_w = (tmp1_dw >> W) + ((tmp1_dw >> W) << k) + (word)tmp1_dw;
  if (tmp2_w >= p)
    return tmp2_w - p;
  return tmp2_w;
}

// result = a.b.2^{-2e} mod p
static inline word montgommery_friendly(word a,
                                        word b,
                                        word k,
                                        word e,
                                        word p,
                                        double_word mask)
{
  double_word c;
  double_word tmp1_dw;
  word tmp2_w;

  c = (double_word)a * b;
  tmp1_dw = (c & mask) * k + (c >> e);
  tmp2_w = (tmp1_dw & mask) * k + (tmp1_dw >> e);
  if (tmp2_w >= p)
  {
    return (tmp2_w - p);
  }
  return tmp2_w;
}

static inline word plantard(word a, word b, word p, double_word r)
{
  double_word c;
  word tmp1_w;
  word tmp2_w;

  c = (double_word)a * b;
  tmp1_w = ((c * r) >> W) + 1;
  tmp2_w = ((double_word)tmp1_w * p) >> W;
  if (tmp2_w == p)
    return 0;
  return tmp2_w;
}

int main()
{
  word a, b, u, v, p, mp, r;
  word g;
  word result;
  double_word r_dw;

  a = 2147483648;
  b = 214748364;
  // p = 18446744073709551599; // 2^64-17 thus -p mod 2^64 = 17
  p = 4294967279; // 2^32-17 = 2^32-2^4-1
  r = 4042322161; // (-p^{-1} mod 2^n)
  mp = WORD_MAX - p + 1;
  // r = 17361641481138401521; // -p^{-1} mod 2^64

  g = gcdExtended(a, b, &u, &v);

  printf("%u %u %u %u \n", g, u, v, (a * u + b * v)); // Test gcdext OK

  result = montgommery(a, b, p, r);
  printf("montgommery %u \n", result);
  r = ((double_word)1 << W) / p;
  // printf("R Barrett %u \n", r);
  result = barrett(a, b, p, r);
  printf("Barrett %u \n", result);
  result = mersenne_32(a, b);
  printf("Mersenne %u \n", result);
  result = pseudo_mersenne(a, b, 17, p);
  printf("Pseudo-Mersenne %u \n", result);
  result = generalized_mersenne(a, b, 4, p);
  printf("Generalized-Mersenne %u \n", result);

  //////// A vérifier ce qui suit  !!!!!
  p = (double_word)17 * (1 << 27) - 1;
  printf("p = %u \n ", p);
  result = montgommery_friendly(a, b, 17, 27, p, ((double_word)1 << 27) - 1);
  printf("montgommery-Friendly %u \n", result);

  p = 4294967279; // 2^32-17 = 2^32-2^4-1

  ///// Plantard OK
  p = 2027808485;             // un peu moins de 2^32/phi
  r_dw = 8520319342442342125; // (p^{-1} mod 2^(2n)
  result = plantard(a, b, p, r_dw);
  printf("plantard %u \n", result);

  //////////////////////////////////////////////////
  // Timings
  //////////////////////////////////////////////////
  unsigned long long meanTimer;
  unsigned long long t1, t2, timer, timer2;
  unsigned long long nb1, nb2;

  //////////////////////////////////////////////////
  // Montgommery
  //////////////////////////////////////////////////
  meanTimer = 0;
  // Heating caches
  for (int i = 0; i < NTEST; i++)
  {
    // appel de la fonction a mesurer a mettre ici
    // juste pour chauffer les caches
    result = montgommery(a, b, p, r);
  }
  r = 4042322161;
  p = 4294967279;
  // timing
  for (int i = 0; i < NSAMPLES; i++)
  {
    // initialiser un nouveau jeu de donnees a tester
    a = rand();
    b = rand();
    timer = (unsigned long long int)0x1 << 63;
    timer2 = (unsigned long long int)0x1 << 63;
    for (int j = 0; j < NTEST; j++)
    {
      t1 = cpucyclesStart();

      // t1 = rdpmc_actual_cycles();
      // nb1 = rdpmc_instructions();
      // appeler la fonction ici avec toujours le meme jeu de donnees
      for (int k = 0; k < 100; k++)
      {
        result = montgommery(a, b, p, r);
      }

      // nb2= rdpmc_instructions();
      // t2 = rdpmc_actual_cycles();
      t2 = cpucyclesStop();
      // if(timer2>nb2-nb1) timer2 = nb2-nb1;
      if (timer > t2 - t1)
        timer = t2 - t1;
    }
    meanTimer += timer;
    // meanTimer2 += timer2;
  }
  printf("\n Montgommery : %lld CPU cycles  \n", meanTimer / (NSAMPLES * 100));

  meanTimer = 0;
  // Heating caches
  for (int i = 0; i < NTEST; i++)
  {
    // appel de la fonction a mesurer a mettre ici
    // juste pour chauffer les caches
    result = montgommery(a, b, p, r);
  }
  r = 4042322161;
  p = 4294967279;
  // timing
  for (int i = 0; i < NSAMPLES; i++)
  {
    // initialiser un nouveau jeu de donnees a tester
    a = rand();
    b = rand();
    timer = (unsigned long long int)0x1 << 63;
    timer2 = (unsigned long long int)0x1 << 63;
    for (int j = 0; j < NTEST; j++)
    {
      // t1 = cpucyclesStart();

      // t1 = rdpmc_actual_cycles();
      nb1 = rdpmc_instructions();
      // appeler la fonction ici avec toujours le meme jeu de donnees
      for (int k = 0; k < 100; k++)
      {
        result = montgommery(a, b, p, r);
      }

      nb2 = rdpmc_instructions();
      // t2 = rdpmc_actual_cycles();
      t2 = cpucyclesStop();
      if (timer2 > nb2 - nb1)
        timer2 = nb2 - nb1;
      // if(timer>t2-t1) timer = t2-t1;
    }
    meanTimer += timer2;
    // meanTimer2 += timer2;
  }
  printf(" Montgommery : %lld instructions  \n", meanTimer / (NSAMPLES * 100));

  //////////////////////////////////////////////////
  // Barrett
  //////////////////////////////////////////////////
  meanTimer = 0;
  // Heating caches
  for (int i = 0; i < NTEST; i++)
  {
    // appel de la fonction a mesurer a mettre ici
    // juste pour chauffer les caches
    result = barrett(a, b, p, r);
  }
  r = ((double_word)1 << 32) / p;
  p = 4294967279;
  // timing
  for (int i = 0; i < NSAMPLES; i++)
  {
    // initialiser un nouveau jeu de donnees a tester
    a = rand();
    b = rand();
    timer = (unsigned long long int)0x1 << 63;
    // timer2 = (unsigned long long int)0x1<<63;
    for (int j = 0; j < NTEST; j++)
    {
      t1 = cpucyclesStart();

      // t1 = rdpmc_actual_cycles();
      // nb1 = rdpmc_instructions();
      // appeler la fonction ici avec toujours le meme jeu de donnees
      for (int k = 0; k < 100; k++)
      {
        result = barrett(a, b, p, r);
      }

      // nb2= rdpmc_instructions();
      // t2 = rdpmc_actual_cycles();
      t2 = cpucyclesStop();
      // if(timer2>nb2-nb1) timer2 = nb2-nb1;
      if (timer > t2 - t1)
        timer = t2 - t1;
    }
    meanTimer += timer;
    // meanTimer2 += timer2;
  }
  printf("\n Barrett : %lld CPU cycles  \n", meanTimer / (NSAMPLES * 100));

  meanTimer = 0;
  // Heating caches
  for (int i = 0; i < NTEST; i++)
  {
    // appel de la fonction a mesurer a mettre ici
    // juste pour chauffer les caches
    result = barrett(a, b, p, r);
  }
  r = ((double_word)1 << 32) / p;
  p = 4294967279;
  // timing
  for (int i = 0; i < NSAMPLES; i++)
  {
    // initialiser un nouveau jeu de donnees a tester
    a = rand();
    b = rand();
    // timer = (unsigned long long int)0x1<<63;
    timer2 = (unsigned long long int)0x1 << 63;
    for (int j = 0; j < NTEST; j++)
    {
      // t1 = cpucyclesStart();

      // t1 = rdpmc_actual_cycles();
      nb1 = rdpmc_instructions();
      // appeler la fonction ici avec toujours le meme jeu de donnees
      for (int k = 0; k < 100; k++)
      {
        result = barrett(a, b, p, r);
      }

      nb2 = rdpmc_instructions();
      // t2 = rdpmc_actual_cycles();
      // t2 = cpucyclesStop();
      if (timer2 > nb2 - nb1)
        timer2 = nb2 - nb1;
      // if(timer>t2-t1) timer = t2-t1;
    }
    meanTimer += timer2;
    // meanTimer2 += timer2;
  }
  printf(" Barrett : %lld instruction  \n", meanTimer / (NSAMPLES * 100));

  //////////////////////////////////////////////////
  // Mersennes
  //////////////////////////////////////////////////
  meanTimer = 0;
  // Heating caches
  for (int i = 0; i < NTEST; i++)
  {
    // appel de la fonction a mesurer a mettre ici
    // juste pour chauffer les caches
    result = mersenne_32(a, b);
  }
  // timing
  for (int i = 0; i < NSAMPLES; i++)
  {
    // initialiser un nouveau jeu de donnees a tester
    a = rand();
    b = rand();
    timer = (unsigned long long int)0x1 << 63;
    // timer2 = (unsigned long long int)0x1<<63;
    for (int j = 0; j < NTEST; j++)
    {
      t1 = cpucyclesStart();

      // t1 = rdpmc_actual_cycles();
      // nb1 = rdpmc_instructions();
      // appeler la fonction ici avec toujours le meme jeu de donnees
      for (int k = 0; k < 100; k++)
      {
        result = mersenne_32(a, b);
      }

      // nb2= rdpmc_instructions();
      // t2 = rdpmc_actual_cycles();
      t2 = cpucyclesStop();
      // if(timer2>nb2-nb1) timer2 = nb2-nb1;
      if (timer > t2 - t1)
        timer = t2 - t1;
    }
    meanTimer += timer;
    // meanTimer2 += timer2;
  }
  printf("\n Mersenne : %lld CPU cycles  \n", meanTimer / (NSAMPLES * 100));

  meanTimer = 0;
  // Heating caches
  for (int i = 0; i < NTEST; i++)
  {
    // appel de la fonction a mesurer a mettre ici
    // juste pour chauffer les caches
    result = mersenne_32(a, b);
  }
  // timing
  for (int i = 0; i < NSAMPLES; i++)
  {
    // initialiser un nouveau jeu de donnees a tester
    a = rand();
    b = rand();
    // timer = (unsigned long long int)0x1<<63;
    timer2 = (unsigned long long int)0x1 << 63;
    for (int j = 0; j < NTEST; j++)
    {
      // t1 = cpucyclesStart();

      // t1 = rdpmc_actual_cycles();
      nb1 = rdpmc_instructions();
      // appeler la fonction ici avec toujours le meme jeu de donnees
      for (int k = 0; k < 100; k++)
      {
        result = mersenne_32(a, b);
      }

      nb2 = rdpmc_instructions();
      // t2 = rdpmc_actual_cycles();
      // t2 = cpucyclesStop();
      if (timer2 > nb2 - nb1)
        timer2 = nb2 - nb1;
      // if(timer>t2-t1) timer = t2-t1;
    }
    meanTimer += timer2;
    // meanTimer2 += timer2;
  }
  printf(" Mersenne : %lld instructions  \n", meanTimer / (NSAMPLES * 100));

  //////////////////////////////////////////////////
  // Pseudo-Mersennes
  //////////////////////////////////////////////////
  meanTimer = 0;
  // Heating caches
  for (int i = 0; i < NTEST; i++)
  {
    // appel de la fonction a mesurer a mettre ici
    // juste pour chauffer les caches
    result = pseudo_mersenne(a, b, 17, p);
  }
  p = 4294967279;
  // timing
  for (int i = 0; i < NSAMPLES; i++)
  {
    // initialiser un nouveau jeu de donnees a tester
    a = rand();
    b = rand();
    timer = (unsigned long long int)0x1 << 63;
    // timer2 = (unsigned long long int)0x1<<63;
    for (int j = 0; j < NTEST; j++)
    {
      t1 = cpucyclesStart();

      // t1 = rdpmc_actual_cycles();
      // nb1 = rdpmc_instructions();
      // appeler la fonction ici avec toujours le meme jeu de donnees
      for (int k = 0; k < 100; k++)
      {
        result = pseudo_mersenne(a, b, 17, p);
      }

      // nb2= rdpmc_instructions();
      // t2 = rdpmc_actual_cycles();
      t2 = cpucyclesStop();
      // if(timer2>nb2-nb1) timer2 = nb2-nb1;
      if (timer > t2 - t1)
        timer = t2 - t1;
    }
    meanTimer += timer;
    // meanTimer2 += timer2;
  }
  printf("\n Pseudo-Mersennes : %lld CPU cycles  \n",
         meanTimer / (NSAMPLES * 100));

  meanTimer = 0;
  // Heating caches
  for (int i = 0; i < NTEST; i++)
  {
    // appel de la fonction a mesurer a mettre ici
    // juste pour chauffer les caches
    result = pseudo_mersenne(a, b, 17, p);
  }
  p = 4294967279;
  // timing
  for (int i = 0; i < NSAMPLES; i++)
  {
    // initialiser un nouveau jeu de donnees a tester
    a = rand();
    b = rand();
    // timer = (unsigned long long int)0x1<<63;
    timer2 = (unsigned long long int)0x1 << 63;
    for (int j = 0; j < NTEST; j++)
    {
      // t1 = cpucyclesStart();

      // t1 = rdpmc_actual_cycles();
      nb1 = rdpmc_instructions();
      // appeler la fonction ici avec toujours le meme jeu de donnees
      for (int k = 0; k < 100; k++)
      {
        result = pseudo_mersenne(a, b, 17, p);
      }

      nb2 = rdpmc_instructions();
      // t2 = rdpmc_actual_cycles();
      // t2 = cpucyclesStop();
      if (timer2 > nb2 - nb1)
        timer2 = nb2 - nb1;
      // if(timer>t2-t1) timer = t2-t1;
    }
    meanTimer += timer2;
    // meanTimer2 += timer2;
  }
  printf(" Pseudo-Mersennes : %lld instructions  \n",
         meanTimer / (NSAMPLES * 100));

  //////////////////////////////////////////////////
  // Generalized-Mersennes
  //////////////////////////////////////////////////
  meanTimer = 0;
  // Heating caches
  for (int i = 0; i < NTEST; i++)
  {
    // appel de la fonction a mesurer a mettre ici
    // juste pour chauffer les caches
    result = generalized_mersenne(a, b, 4, p);
  }
  r = ((double_word)1 << 32) / p;
  p = 4294967279;
  // timing
  for (int i = 0; i < NSAMPLES; i++)
  {
    // initialiser un nouveau jeu de donnees a tester
    a = rand();
    b = rand();
    timer = (unsigned long long int)0x1 << 63;
    // timer2 = (unsigned long long int)0x1<<63;
    for (int j = 0; j < NTEST; j++)
    {
      t1 = cpucyclesStart();

      // t1 = rdpmc_actual_cycles();
      // nb1 = rdpmc_instructions();
      // appeler la fonction ici avec toujours le meme jeu de donnees
      for (int k = 0; k < 100; k++)
      {
        result = generalized_mersenne(a, b, 4, p);
      }

      // nb2= rdpmc_instructions();
      // t2 = rdpmc_actual_cycles();
      t2 = cpucyclesStop();
      // if(timer2>nb2-nb1) timer2 = nb2-nb1;
      if (timer > t2 - t1)
        timer = t2 - t1;
    }
    meanTimer += timer;
    // meanTimer2 += timer2;
  }
  printf("\n Generalized-Mersennes : %lld CPU cycles  \n",
         meanTimer / (NSAMPLES * 100));

  meanTimer = 0;
  // Heating caches
  for (int i = 0; i < NTEST; i++)
  {
    // appel de la fonction a mesurer a mettre ici
    // juste pour chauffer les caches
    result = generalized_mersenne(a, b, 4, p);
  }
  r = ((double_word)1 << 32) / p;
  p = 4294967279;
  // timing
  for (int i = 0; i < NSAMPLES; i++)
  {
    // initialiser un nouveau jeu de donnees a tester
    a = rand();
    b = rand();
    // timer = (unsigned long long int)0x1<<63;
    timer2 = (unsigned long long int)0x1 << 63;
    for (int j = 0; j < NTEST; j++)
    {
      // t1 = cpucyclesStart();

      // t1 = rdpmc_actual_cycles();
      nb1 = rdpmc_instructions();
      // appeler la fonction ici avec toujours le meme jeu de donnees
      for (int k = 0; k < 100; k++)
      {
        result = generalized_mersenne(a, b, 4, p);
      }

      nb2 = rdpmc_instructions();
      // t2 = rdpmc_actual_cycles();
      // t2 = cpucyclesStop();
      if (timer2 > nb2 - nb1)
        timer2 = nb2 - nb1;
      // if(timer>t2-t1) timer = t2-t1;
    }
    meanTimer += timer2;
    // meanTimer2 += timer2;
  }
  printf(" Generalized-Mersennes : %lld instructions  \n",
         meanTimer / (NSAMPLES * 100));

  //////////////////////////////////////////////////
  // Montgommery friendly
  //////////////////////////////////////////////////
  meanTimer = 0;
  // Heating caches
  for (int i = 0; i < NTEST; i++)
  {
    // appel de la fonction a mesurer a mettre ici
    // juste pour chauffer les caches
    result = montgommery_friendly(a, b, 17, 27, p, ((double_word)1 << 27) - 1);
  }
  p = (double_word)17 * (1 << 27) - 1;
  result = montgommery_friendly(a, b, 17, 27, p, ((double_word)1 << 27) - 1);

  // timing
  for (int i = 0; i < NSAMPLES; i++)
  {
    // initialiser un nouveau jeu de donnees a tester
    a = rand();
    b = rand();
    timer = (unsigned long long int)0x1 << 63;
    // timer2 = (unsigned long long int)0x1<<63;
    for (int j = 0; j < NTEST; j++)
    {
      t1 = cpucyclesStart();

      // t1 = rdpmc_actual_cycles();
      // nb1 = rdpmc_instructions();
      // appeler la fonction ici avec toujours le meme jeu de donnees
      for (int k = 0; k < 100; k++)
      {
        result =
            montgommery_friendly(a, b, 17, 27, p, ((double_word)1 << 27) - 1);
      }

      // nb2= rdpmc_instructions();
      // t2 = rdpmc_actual_cycles();
      t2 = cpucyclesStop();
      // if(timer2>nb2-nb1) timer2 = nb2-nb1;
      if (timer > t2 - t1)
        timer = t2 - t1;
    }
    meanTimer += timer;
    // meanTimer2 += timer2;
  }
  printf("\n Montgommery friendly : %lld CPU cycles  \n",
         meanTimer / (NSAMPLES * 100));

  meanTimer = 0;
  // Heating caches
  for (int i = 0; i < NTEST; i++)
  {
    // appel de la fonction a mesurer a mettre ici
    // juste pour chauffer les caches
    result = montgommery_friendly(a, b, 17, 27, p, ((double_word)1 << 27) - 1);
  }
  p = (double_word)17 * (1 << 27) - 1;
  result = montgommery_friendly(a, b, 17, 27, p, ((double_word)1 << 27) - 1);

  // timing
  for (int i = 0; i < NSAMPLES; i++)
  {
    // initialiser un nouveau jeu de donnees a tester
    a = rand();
    b = rand();
    // timer = (unsigned long long int)0x1<<63;
    timer2 = (unsigned long long int)0x1 << 63;
    for (int j = 0; j < NTEST; j++)
    {
      // t1 = cpucyclesStart();

      // t1 = rdpmc_actual_cycles();
      nb1 = rdpmc_instructions();
      // appeler la fonction ici avec toujours le meme jeu de donnees
      for (int k = 0; k < 100; k++)
      {
        result =
            montgommery_friendly(a, b, 17, 27, p, ((double_word)1 << 27) - 1);
      }

      nb2 = rdpmc_instructions();
      // t2 = rdpmc_actual_cycles();
      // t2 = cpucyclesStop();
      if (timer2 > nb2 - nb1)
        timer2 = nb2 - nb1;
      // if(timer>t2-t1) timer = t2-t1;
    }
    meanTimer += timer2;
    // meanTimer2 += timer2;
  }
  printf(" Montgommery friendly : %lld instructions  \n",
         meanTimer / (NSAMPLES * 100));

  //////////////////////////////////////////////////
  // Plantard
  //////////////////////////////////////////////////
  meanTimer = 0;
  // Heating caches
  for (int i = 0; i < NTEST; i++)
  {
    // appel de la fonction a mesurer a mettre ici
    // juste pour chauffer les caches
    result = plantard(a, b, p, r);
  }
  p = 2027808485;             // un peu moins de 2^32/phi
  r_dw = 8520319342442342125; // (p^{-1} mod 2^(2n)

  // timing
  for (int i = 0; i < NSAMPLES; i++)
  {
    // initialiser un nouveau jeu de donnees a tester
    a = rand();
    b = rand();
    timer = (unsigned long long int)0x1 << 63;
    // timer2 = (unsigned long long int)0x1<<63;
    for (int j = 0; j < NTEST; j++)
    {
      t1 = cpucyclesStart();

      // t1 = rdpmc_actual_cycles();
      // nb1 = rdpmc_instructions();
      // appeler la fonction ici avec toujours le meme jeu de donnees
      for (int k = 0; k < 100; k++)
      {
        result = plantard(a, b, p, r);
      }

      // nb2= rdpmc_instructions();
      // t2 = rdpmc_actual_cycles();
      t2 = cpucyclesStop();
      // if(timer2>nb2-nb1) timer2 = nb2-nb1;
      if (timer > t2 - t1)
        timer = t2 - t1;
    }
    meanTimer += timer;
    // meanTimer2 += timer2;
  }
  printf("\n Plantard : %lld CPU cycles  \n", meanTimer / (NSAMPLES * 100));

  meanTimer = 0;
  // Heating caches
  for (int i = 0; i < NTEST; i++)
  {
    // appel de la fonction a mesurer a mettre ici
    // juste pour chauffer les caches
    result = plantard(a, b, p, r_dw);
  }
  // timing
  for (int i = 0; i < NSAMPLES; i++)
  {
    // initialiser un nouveau jeu de donnees a tester
    a = rand();
    b = rand();
    // timer = (unsigned long long int)0x1<<63;
    timer2 = (unsigned long long int)0x1 << 63;
    for (int j = 0; j < NTEST; j++)
    {
      // t1 = cpucyclesStart();

      // t1 = rdpmc_actual_cycles();
      nb1 = rdpmc_instructions();
      // appeler la fonction ici avec toujours le meme jeu de donnees
      for (int k = 0; k < 100; k++)
      {
        result = plantard(a, b, p, r_dw);
      }

      nb2 = rdpmc_instructions();
      // t2 = rdpmc_actual_cycles();
      // t2 = cpucyclesStop();
      if (timer2 > nb2 - nb1)
        timer2 = nb2 - nb1;
      // if(timer>t2-t1) timer = t2-t1;
    }
    meanTimer += timer2;
    // meanTimer2 += timer2;
  }
  printf(" Plantard : %lld instructions  \n", meanTimer / (NSAMPLES * 100));
}
