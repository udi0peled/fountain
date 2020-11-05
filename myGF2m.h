#ifndef __FOUNTAIN_MY_SSL__
#define __FOUNTAIN_MY_SSL__

#include <stdint.h>

typedef int BIGNUM;
typedef int BN_CTX;
typedef unsigned long BN_ULONG;

int RAND_bytes(unsigned char *buf, int num);

BN_CTX *BN_CTX_new(void);
BN_CTX *BN_CTX_secure_new(void);
BIGNUM *BN_CTX_get(BN_CTX *ctx);
void    BN_CTX_start(BN_CTX *ctx);
void    BN_CTX_end(BN_CTX *ctx);
void    BN_CTX_free(BN_CTX *ctx);

int     BN_set_word(BIGNUM *a, BN_ULONG w);
int     BN_bn2binpad(const BIGNUM *a, unsigned char *to, int tolen);
BIGNUM *BN_bin2bn(const unsigned char *s, int len, BIGNUM *ret);
char   *BN_bn2dec(const BIGNUM *a);
int     BN_num_bits(const BIGNUM *a);
# define BN_num_bytes(a) ((BN_num_bits(a)+7)/8)


int BN_GF2m_mod_mul_arr(BIGNUM *r, const BIGNUM *a, const BIGNUM *b, const int p[], BN_CTX *ctx);
int BN_GF2m_mod_inv_arr(BIGNUM *r, const BIGNUM *xx, const int p[], BN_CTX *ctx);
int BN_GF2m_add(BIGNUM *r, const BIGNUM *a, const BIGNUM *b);
#  define BN_GF2m_sub(r, a, b) BN_GF2m_add(r, a, b);

#endif