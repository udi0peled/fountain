#include <assert.h>
#include <string.h>
#include <openssl/rand.h>
#include <time.h>
#include <stdio.h>

#define PRIME_FIELD 4294967311;
#define CHUNK_BYTES 4;


typedef struct
{ 
    uint64_t base_data_len;
    uint64_t *base_data;
    uint64_t *lagrange;
} fountain_ctx;

fountain_ctx *fountain_ctx_new(uint64_t chunk_bytelen){
    fountain_ctx *ctx = malloc(sizeof(fountain_ctx));
    ctx->base_data_len = 0;
    ctx->base_data = NULL;
    ctx->lagrange = NULL;
    return ctx;
}

void fountain_ctx_free(fountain_ctx *ctx) {
    free(ctx->base_data);
    free(ctx->base_data);
    free(ctx);
}

// C function for extended Euclidean Algorithm  
uint64_t gcd_extended(uint64_t a, uint64_t b, uint64_t *x, uint64_t *y);  
  
// Function to find modulo inverse of a  
void mod_inverse(uint64_t a, uint64_t m)  
{  
    uint64_t x, y;  
    uint64_t g = gcd_extended(a, m, &x, &y);  
    if (g != 1)  
        printf("Inverse doesn't exist\n");  
    else
    {  
        // m is added to handle negative x  
        uint64_t res = (x%m + m) % m;  
        printf("Modular multiplicative inverse is %ld\n",res);  
    }  
}  
  
// C function for extended Euclidean Algorithm  
uint64_t gcd_extended(uint64_t a, uint64_t b, uint64_t *x, uint64_t *y)  
{  
    // Base Case  
    if (a == 0)  
    {  
        *x = 0, *y = 1;  
        return b;  
    }  
  
    uint64_t x1, y1; // To store results of recursive call  
    uint64_t gcd = gcd_extended(b%a, a, &x1, &y1);  
  
    // Update x and y using results of recursive  
    // call  
    *x = y1 - (b * x1 / a);  
    *y = x1;  
  
    return gcd;  
}  
  
void set_base_data(fountain_ctx *ctx, uint32_t *base_data, uint64_t base_data_len) {
    uint64_t temp;

    ctx->base_data_len = base_data_len;
    ctx->base_data = calloc(base_data_len, sizeof(uint64_t));
    ctx->lagrange = calloc(base_data_len, sizeof(uint64_t));

    for (uint64_t i = 0; i < base_data_len; ++i) {
        ctx->base_data[i] = base_data[i] % PRIME_FIELD;
        
        // Set lagrange multiplier denominator (inverse once at end)
        ctx->lagrange[i] = 1;
        for (uint64_t j = 0; j < base_data_len; ++j) {
            if (j == i) continue;
            temp = (i - j) % PRIME_FIELD;
            (ctx->lagrange[i] * temp) % PRIME_FIELD;
        }
    }
}
/*
uint64_t compute_chunk_at(uint64_t chunk_index, fountain_ctx *ctx)
{

    uint64_t chunk_val = 0;
    uint64_t lagrange;
    uint64_t temp;
    
    for (uint64_t i = 0; i < ctx->base_data_len; ++i) {        
        lagrange = ctx->base_data[i];

        for (uint64_t j = 0; j < ctx->base_data_len; ++j) {
            if (j == i) continue;
            temp = (chunk_index - i) % PRIME_FIELD;
            lagrange = (lagrange * temp) % PRIME_FIELD;
        }
        
        chunk_val = (chunk_val + lagrange) % PRIME_FIELD;

    }

    assert((uint64_t) BN_num_bytes(chunk_val) <= ctx->chunk_bytelen);
    BN_bn2binpad(chunk_val, chunk_bytes, ctx->chunk_bytelen);

    BN_free(chunk_ind);
    BN_free(chunk_val);
    BN_free(temp);
    BN_free(lagrange);
}
*/
int main(int argc, char* argv[]) {
    clock_t start, diff;
    double time_ms;

    uint64_t chunk_byte_len;
    uint64_t base_length;
    uint64_t num_chunks;

    if (argc < 2) {
        printf("usage: ./main <chunk_bytes>\n");
        return 0;
    }

    num_chunks  = strtoul(argv[1], NULL, 10);
/*
    // Setup the fountain context and print values
    fountain_ctx *ctx = fountain_ctx_new(chunk_byte_len);
    
    /*
    // Set random base
    uint8_t **base = calloc(base_length, sizeof(uint8_t *));
    for (uint64_t i = 0; i < base_length; ++i)
    {
        base[i] = malloc(ctx->chunk_bytelen);
        RAND_bytes(base[i], ctx->chunk_bytelen);
        // printf("base[%ld]=", i);
        // printBIGNUM("", temp, " ");
    }
    //printf("\n");
    
    // Set base data
    set_base_data(ctx, base, base_length);

    // Test base chunks are correct
    uint8_t *chunk_bytes = malloc(ctx->chunk_bytelen);
    for (uint64_t i = 0; i < base_length; ++i) {
        compute_chunk_at(i, chunk_bytes, ctx);
        assert(memcmp(base[i], chunk_bytes, ctx->chunk_bytelen) == 0);
    }
*/

    start = clock();
    // Compute all needed chunks (also after base)
    uint64_t a = 2;
    uint64_t b = 3;
    for (uint64_t i = 0; i < num_chunks; ++i) {
        a = (b*b) % PRIME_FIELD;
        for (uint64_t j = 0; j < num_chunks; ++j) {
            a = (a * b) % PRIME_FIELD;
        }
        a = (a*a) % PRIME_FIELD;
    }

    diff = clock() - start;
    time_ms = ((double) diff * 1000/ CLOCKS_PER_SEC);
    
    printf("Time computing %ld chunks: %.3f ms\n", num_chunks, time_ms);
}