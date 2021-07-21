# Precise timing of vectorized RNS operations

## Base conversion with Crandall's moduli

### Int64_t version

##### Full function

```C
void base_conversion_cr(int64_t *rop, struct conv_base_t *conv_base, int64_t *op, int64_t *a)
{
	int i, j;
	int64_t tmp;
	int128 tmp2;
	int128 tmp3;
	int64_t up, up2, lo, lo2;
	int64_t mask = ((int64_t)1 << 63) - 1;
	int size = conv_base->rns_a->size;

	for (j = 0; j < size; j++)
		a[j] = op[j];
	for (i = 0; i < size - 1; i++)
	{
		for (j = i + 1; j < size; j++)
		{
			tmp = a[j] - a[i];
			a[j] = mul_mod_cr_t(tmp, conv_base->inva_to_b[i][j], conv_base->rns_a->k[j]);
		}
	}

	for (j = 0; j < size; j++)
		rop[j] = (a[0] > conv_base->rns_b->m[j]) ? a[0] - conv_base->rns_b->m[j] : a[0];
	for (j = 0; j < size; j++)
	{
		for (i = 1; i < size; i++)
		{
			tmp = mul_mod_cr(a[i], conv_base->mrsa_to_b[i - 1][j], conv_base->rns_b->k[j]);
			rop[j] = add_mod_cr(rop[j], tmp, conv_base->rns_b->k[j]);
		}
	}
}
```

##### Initialisation
***172 instructions***
```C
int i, j;
int64_t tmp;
int128 tmp2;
int128 tmp3;
int64_t up, up2, lo, lo2;
int64_t mask = ((int64_t)1 << 63) - 1;
int size = conv_base->rns_a->size;
for (j = 0; j < size; j++)
	a[j] = op[j];
```

##### Mixed Radix
***1352 instructions***
```C
for (i = 0; i < size - 1; i++)
{
	for (j = i + 1; j < size; j++)
	{
		tmp = a[j] - a[i];
		a[j] = mul_mod_cr_t(tmp, conv_base->inva_to_b[i][j], conv_base->rns_a->k[j]);
	}
}
```

##### Operations
***2755 instructions***
```C
for (j = 0; j < size; j++)
{
	for (i = 1; i < size; i++)
	{
		tmp = mul_mod_cr(a[i], conv_base->mrsa_to_b[i - 1][j], conv_base->rns_b->k[j]);
		rop[j] = add_mod_cr(rop[j], tmp, conv_base->rns_b->k[j]);
	}
}
```

### AVX-2 version

```C
inline void avx_base_conversion_cr(__m256i *rop, struct conv_base_t *conv_base, __m256i *op, int64_t *a)
{
	int i, j;
	__m256i avx_tmp;
	int64_t tmp;
	int size = conv_base->rns_a->size;

	for (i = 0; i < size - 1; i++)
	{
		for (j = i + 1; j < size; j++)
		{
			tmp = a[j] - a[i];
			a[j] = mul_mod_cr(tmp, conv_base->inva_to_b[i][j], conv_base->rns_a->k[j]);
		}
	}

	__m256i a0_256 = _mm256_set1_epi64x(a[0]);

	for (j = 0; j < size / 4; j++)
	{
		__m256i b256 = _mm256_lddqu_si256((__m256i *)&conv_base->rns_b->m[j >> 2]);
		__m256i cmp256 = _mm256_cmpgt_epi64(a0_256, b256);
		rop[j] = _mm256_sub_epi64(a0_256, _mm256_and_si256(cmp256, b256));
	}

	for (j = 0; j < size / 4; j++)
	{
		for (i = 1; i < size; i++)
		{
			avx_tmp = avx_mul_mod_cr(_mm256_set1_epi64x(a[i]), conv_base->avx_mrsa_to_b[i - 1][j], conv_base->rns_b->avx_k[j]);
			rop[j] = avx_add_mod_cr(rop[j], avx_tmp, conv_base->rns_b->avx_k[j]);
		}
	}
}
```
##### Initialisation
***36 instructions***
```C
int i, j;
__m256i avx_tmp;
int64_t tmp;
int size = conv_base->rns_a->size;
```

##### Mixed Radix
***1386 instructions***
```C
for (i = 0; i < size - 1; i++)
{
	for (j = i + 1; j < size; j++)
	{
		tmp = a[j] - a[i];
		a[j] = mul_mod_cr_t(tmp, conv_base->inva_to_b[i][j], conv_base->rns_a->k[j]);
	}
}
```

##### Operations
***1443 instructions***
```C
for (j = 0; j < size; j++)
{
	for (i = 1; i < size; i++)
	{
		tmp = mul_mod_cr(a[i], conv_base->mrsa_to_b[i - 1][j], conv_base->rns_b->k[j]);
		rop[j] = add_mod_cr(rop[j], tmp, conv_base->rns_b->k[j]);
	}
}
```

## Base conversion with Cox-Rower method

### Int64_t version

##### Full function

```C
void base_conversion_cox(int64_t *rop, struct conv_base_t *conv_base, int64_t *op, int r, int q, int alpha)
{
	int i, j;
	int n, k_i;
	int sigma;
	int64_t mask, mask2, xhi, trunk;
	int size = conv_base->rns_a->size;
	int64_t tmp, tmp2, tmp3;

	r = 63;							 
	q = 7;							 
	alpha = ((int64_t)1 << (q - 1));

	mask = ((int64_t)1 << r) - ((int64_t)1 << (r - q));
	mask2 = ((int64_t)1 << q);
	n = conv_base->rns_a->size;

	sigma = alpha;
	
	for (i = 0; i < n; i++)
		rop[i] = 0;

	for (i = 0; i < n; i++)
	{

		xhi = mul_mod_cr(op[i], conv_base->rns_a->int_inv_Mi[i], conv_base->rns_a->k[i]);
		trunk = xhi & mask;

		sigma += trunk >> (r - q);
		k_i = sigma & mask2;

		sigma[%d]=%d   k[%d]=%d \n ", i, op[i], i, conv_base->rns_a->int_inv_Mi[i], i, xhi, i, trunk, i, sigma, i, k_i);

		sigma -= k_i;
		k_i = k_i >> q;

		for (j = 0; j < size; j++)
		{
			tmp = mul_mod_cr(xhi, conv_base->Mi_modPi[i][j], conv_base->rns_b->k[j]);
			tmp2 = conv_base->invM_modPi[j] * k_i;
			tmp3 = add_mod_cr(tmp, tmp2, conv_base->rns_b->k[j]);
			rop[j] = add_mod_cr(rop[j], tmp3, conv_base->rns_b->k[j]);
		}
	}
}
```

##### Initialisation
***160 instructions***
```C
int i, j;
int n, k_i;
int sigma;
int64_t mask, mask2, xhi, trunk;
int size = conv_base->rns_a->size;
int64_t tmp, tmp2, tmp3;

r = 63;							 
q = 7;							 
alpha = ((int64_t)1 << (q - 1));

mask = ((int64_t)1 << r) - ((int64_t)1 << (r - q));
mask2 = ((int64_t)1 << q);
n = conv_base->rns_a->size;

sigma = alpha;

for (i = 0; i < n; i++)
	rop[i] = 0;
```

##### Operations
***3780 instructions***
```C
for (i = 0; i < n; i++)
	{

	xhi = mul_mod_cr(op[i], conv_base->rns_a->int_inv_Mi[i], conv_base->rns_a->k[i]);
	trunk = xhi & mask;

	sigma += trunk >> (r - q);
	k_i = sigma & mask2;

	sigma[%d]=%d   k[%d]=%d \n ", i, op[i], i, conv_base->rns_a->int_inv_Mi[i], i, xhi, i, trunk, i, sigma, i, k_i);

	sigma -= k_i;
	k_i = k_i >> q;

	for (j = 0; j < size; j++)
	{
		tmp = mul_mod_cr(xhi, conv_base->Mi_modPi[i][j], conv_base->rns_b->k[j]);
		tmp2 = conv_base->invM_modPi[j] * k_i;
		tmp3 = add_mod_cr(tmp, tmp2, conv_base->rns_b->k[j]);
		rop[j] = add_mod_cr(rop[j], tmp3, conv_base->rns_b->k[j]);
	}
}
```

### AVX-2 version

```C
static __m256i cox_mask = (__m256i){0x7F00000000000000UL, 0x7F00000000000000UL, 0x7F00000000000000UL, 0x7F00000000000000UL};
static __m256i cox_mask2 = (__m256i){0x80UL, 0x80UL, 0x80UL, 0x80UL};
static __m256i cox_sigma = (__m256i){0x40UL, 0x40UL, 0x40UL, 0x40UL};

inline void avx_base_conversion_cox(__m256i *rop, struct conv_base_t *conv_base, __m256i *op, int64_t *a)
{
	int i, j;
	int size = conv_base->rns_a->size;

	int r = 63;
	int q = 7;

	__m256i sigma = cox_sigma;

	for (i = 0; i < size / 4; i++)
	{
		rop[i] = _mm256_set1_epi64x(0);
	}

	__m256i xhi, trunk, k_i;

	__m256i tmp0, tmp1, tmp2;

	for (i = 0; i < size; i++)
	{
		xhi = _mm256_set1_epi64x(mul_mod_cr(a[i], conv_base->rns_a->int_inv_Mi[i], conv_base->rns_a->k[i]));

		trunk = xhi & cox_mask;
		sigma = _mm256_add_epi64(sigma, _mm256_srli_epi64(trunk, r - q));
		k_i = sigma & cox_mask2;

		sigma = _mm256_sub_epi64(sigma, k_i);
		k_i = _mm256_srli_epi64(k_i, q);

		for (j = 0; j < size / 4; j++)
		{
			tmp0 = avx_mul_mod_cr(xhi, conv_base->avx_Mi_modPi[i][j], conv_base->rns_b->avx_k[j]);
			tmp1 = conv_base->avx_invM_modPi[j] * k_i;
			tmp2 = avx_add_mod_cr(tmp0, tmp1, conv_base->rns_b->avx_k[j]);
			rop[j] = avx_add_mod_cr(rop[j], tmp2, conv_base->rns_b->avx_k[j]);
		}
	}
}
```
##### Initialisation
***85 instructions***

```C
int i, j;
int size = conv_base->rns_a->size;

int r = 63;
int q = 7;

__m256i sigma = cox_sigma;

for (i = 0; i < size / 4; i++)
{
	rop[i] = _mm256_set1_epi64x(0);
}
__m256i xhi, trunk, k_i;
__m256i tmp0, tmp1, tmp2;
```

##### Operations
***2347 instructions***
```C
for (i = 0; i < size; i++)
{
	xhi = _mm256_set1_epi64x(mul_mod_cr(a[i], conv_base->rns_a->int_inv_Mi[i], conv_base->rns_a->k[i]));

	trunk = xhi & cox_mask;
	sigma = _mm256_add_epi64(sigma, _mm256_srli_epi64(trunk, r - q));
	k_i = sigma & cox_mask2;

	sigma = _mm256_sub_epi64(sigma, k_i);
	k_i = _mm256_srli_epi64(k_i, q);

	for (j = 0; j < size / 4; j++)
	{
		tmp0 = avx_mul_mod_cr(xhi, conv_base->avx_Mi_modPi[i][j], conv_base->rns_b->avx_k[j]);
		tmp1 = conv_base->avx_invM_modPi[j] * k_i;
		tmp2 = avx_add_mod_cr(tmp0, tmp1, conv_base->rns_b->avx_k[j]);
		rop[j] = avx_add_mod_cr(rop[j], tmp2, conv_base->rns_b->avx_k[j]);
	}
}
```
