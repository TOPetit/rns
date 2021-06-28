#include "structs_data.h"

// RNS arithmetic functions
void add_rns(int64_t *rop, struct rns_base_t *base, int64_t *pa, int64_t *pb);
void sub_rns(int64_t *rop, struct rns_base_t *base, int64_t *pa, int64_t *pb);
void mul_rns(int64_t *rop, struct rns_base_t *base, int64_t *pa, int64_t *pb, int64_t n);
void mult_mod_rns(int64_t *rop, int64_t *pa, int64_t *pab, int64_t *pb, int64_t *pbb, struct mod_mul_t *mult, int64_t *tmp[3]);

void init_rns(struct rns_base_t *base);
void from_int_to_rns(int64_t *rop, struct rns_base_t *base, mpz_t op);
