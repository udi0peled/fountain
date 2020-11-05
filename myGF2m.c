#include "my_ssl.h"
#include <stdlib.h>

int RAND_bytes(unsigned char *buf, int num) {
    return 0;
}

BN_CTX *BN_CTX_new(void) {
    return NULL;
}

BN_CTX *BN_CTX_secure_new(void) {
    return NULL;
}

BIGNUM *BN_CTX_get(BN_CTX *ctx) {
    return NULL;
}

void    BN_CTX_start(BN_CTX *ctx) {
    return;
}

void    BN_CTX_end(BN_CTX *ctx) {
    return ;
}

void    BN_CTX_free(BN_CTX *ctx) {
    return;
}

int     BN_set_word(BIGNUM *a, BN_ULONG w) {
    return 0;
}

int     BN_bn2binpad(const BIGNUM *a, unsigned char *to, int tolen) {
    return 0;
}

BIGNUM *BN_bin2bn(const unsigned char *s, int len, BIGNUM *ret) {
    return NULL;
}

char   *BN_bn2dec(const BIGNUM *a) {
    return NULL;
}

int     BN_num_bits(const BIGNUM *a) {
    return 0;
}

# define BN_num_bytes(a) ((BN_num_bits(a)+7)/8)


int BN_GF2m_mod_mul_arr(BIGNUM *r, const BIGNUM *a, const BIGNUM *b, const int p[], BN_CTX *ctx) {
    return 0;
}

int BN_GF2m_mod_inv_arr(BIGNUM *r, const BIGNUM *xx, const int p[], BN_CTX *ctx) {
    return 0;
}

int BN_GF2m_add(BIGNUM *r, const BIGNUM *a, const BIGNUM *b) {
    return 0;
}