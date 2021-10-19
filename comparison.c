#include "galois_16bit_log_table.h"
#include "galois_16bit_inv_log_table.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <time.h>
#include <inttypes.h>
#include <openssl/bn.h>

void printHexBytes(const char * prefix, const char *src, unsigned len, const char * suffix, int print_len) {
  if (len == 0) {
    printf("%s <0 len char array> %s", prefix, suffix);
    return;
  }

  if (print_len) printf("[%u]", len);
  printf("%s", prefix);
  unsigned int i;
  for (i = 0; i < len-1; ++i) {
    printf("%02x",src[i] & 0xff);
  }
  printf("%02x%s",src[i] & 0xff, suffix);
}


void usage_error() {
    printf("usage: <program> <number of repititions>\n");
    exit(1);
}

/****************************************
 * 
 *   Galois function on 16-bit field
 * 
 * *********************************** */ 

int galois_mult(int x, int y)
{
  if (x == 0 || y == 0) return 0;
  
  int sum_j = galois_16bit_log_table[x] + galois_16bit_log_table[y];
  
  return galois_16bit_inv_log_table[sum_j];
}

int galois_div(int a, int b)
{
    int sum_j;

    if (b == 0) return 0;
    if (a == 0) return 0;

    sum_j = galois_16bit_log_table[a] - galois_16bit_log_table[b];
    if (sum_j < 0) sum_j += (int) (1 << 16)-1;
    return galois_16bit_inv_log_table[sum_j];
}

typedef unsigned short field_el;


void time_access_vs_addition(uint64_t num) {
    
    printf("\nAccess vs Addition...\n");

    clock_t start, diff;
    double time_ms;

    field_el a, b;
    a = rand();
    b = rand();

    printf("Accessing table %lu times...\n", num);
    start = clock();

    for (uint32_t i = 0; i < num; ++i) {
        a = galois_16bit_log_table[a] + b;    
    }

    diff = clock() - start;
    time_ms = ((double) diff * 1000/ CLOCKS_PER_SEC);
    printf("Done. Time: %.3f ms\n", time_ms);


    printf("Adding %lu times...\n", num);
    start = clock();
    for (uint32_t i = 0; i < num; ++i) {
        b = galois_16bit_log_table[b] + a;
    }

    diff = clock() - start;
    time_ms = ((double) diff * 1000/ CLOCKS_PER_SEC);
    printf("Done. Time: %.3f ms\n", time_ms);

    if (b == 0) printf("hello\n");
    if (a == 0) printf("world\n");
}

void time_galois_16bit_log_table(uint64_t num) {
    
    printf("\nGalois 16bit table operations...\n");

    clock_t start, diff;
    double time_ms;

    field_el a, b;
    a = rand();
    b = rand();

    printf("Multiplying %lu times...\n", num);
    start = clock();

    for (uint32_t i = 0; i < num; ++i) {
        b = galois_mult(a, b);
    }

    diff = clock() - start;
    time_ms = ((double) diff * 1000/ CLOCKS_PER_SEC);
    printf("Done. Time: %.3f ms\n", time_ms);


    printf("Dividing %lu times...\n", num);
    start = clock();
    for (uint32_t i = 0; i < num; ++i) {
        b = galois_div(b, a);
    }

    diff = clock() - start;
    time_ms = ((double) diff * 1000/ CLOCKS_PER_SEC);
    printf("Done. Time: %.3f ms\n", time_ms);

    
    if (b == 0) printf("hello\n");
    if (a == 0) printf("world\n");
}

void time_openssl_gf2m(uint64_t bits, uint64_t num) {

    int irr_poly_16[]   = {16, 5 ,3, 1, 0, -1};
    int irr_poly_32[]   = {32, 7 ,3, 2, 0, -1};
    int irr_poly_64[]   = {64, 4 ,3, 1, 0, -1};
    int irr_poly_128[]  = {128, 7, 2, 1, 0, -1};
    int irr_poly_256[]  = {256, 10, 5, 2, 0, -1};

    printf("\nOpenSSL GF2^%lu operations...\n", bits);

    clock_t start, diff;
    double time_ms;

    int field[6];

    if (bits == 16) {
      memcpy(field, irr_poly_16, sizeof(field));
    } else if (bits == 32) {
        memcpy(field, irr_poly_32, sizeof(field));
    } else if (bits == 64) {
        memcpy(field, irr_poly_64, sizeof(field));
    } else if (bits == 128) {
        memcpy(field, irr_poly_128, sizeof(field));
    } else if (bits == 256) {
        memcpy(field, irr_poly_256, sizeof(field));
    } else {
      exit(1);
    }

    BN_CTX *bn_ctx = BN_CTX_new();
    BIGNUM *a = BN_CTX_get(bn_ctx);
    BIGNUM *b = BN_CTX_get(bn_ctx);

    BN_rand(a, bits, 1, 0);
    BN_rand(b, bits, 1, 0);

    printf("Multiplying %lu times...\n", num);
    start = clock();

    for (uint32_t i = 0; i < num; ++i) {
       BN_GF2m_mod_mul_arr(b, a, b, field, bn_ctx);
    }

    diff = clock() - start;
    time_ms = ((double) diff * 1000/ CLOCKS_PER_SEC);
    printf("Done. Time: %.3f ms\n", time_ms);

    if (BN_is_zero(a)) printf("Hello\n");
    if (BN_is_zero(b)) printf("World\n");
    // printf("Dividing %lu times...\n", num);
    // start = clock();
    // for (uint32_t i = 0; i < num; ++i) {
    //     BN_GF2m_mod_div_arr(b, b, a, field, bn_ctx);
    // }

    // diff = clock() - start;
    // time_ms = ((double) diff * 1000/ CLOCKS_PER_SEC);
    // printf("Done. Time: %.3f ms\n", time_ms);

    BN_CTX_free(bn_ctx);
}

int main(int argc, char* argv[]) {

    srand(time(NULL));

    if (argc < 2) usage_error();

    uint64_t num_reps = strtoul(argv[1], NULL, 10);

    time_access_vs_addition(num_reps);
    time_galois_16bit_log_table(num_reps);
    time_openssl_gf2m(16, num_reps);
    time_openssl_gf2m(32, num_reps);
    time_openssl_gf2m(64, num_reps);
    time_openssl_gf2m(128, num_reps);
    time_openssl_gf2m(256, num_reps);

}