#include <reed_solomon.h>
#include <assert.h>
#include <string.h>
#include <openssl/rand.h>
#include <time.h>

typedef struct
{ 
    BIGNUM *field;
    BN_CTX *bn_ctx;
    uint64_t chunk_bytelen;
    uint64_t base_data_len;
    BIGNUM **base_data;
    BIGNUM **lagrange;

} reed_solomon_ctx;

void printHexBytes(const char * prefix, const uint8_t *src, unsigned len, const char * suffix, int print_len) {
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

void printBIGNUM(const char * prefix, const BIGNUM *bn, const char * suffix) {
  char *bn_str = BN_bn2dec(bn);
  printf("%s%s%s", prefix, bn_str, suffix);
  free(bn_str);
}

reed_solomon_ctx *reed_solomon_field_ctx_new(uint64_t needed_chunk_bytelen){
    BN_CTX *bn_ctx = BN_CTX_new();

    int irr_poly_2048[] = {2048, 19, 14, 13, 0, -1};
    int irr_poly_4096[] = {4096, 27, 15, 1 , 0, -1};
    int irr_poly_8192[] = {8192, 9 , 5 , 2 , 0, -1};

    BIGNUM *field = BN_new();
    if (needed_chunk_bytelen * 8 > 8192) {
        return NULL;
    } else if (needed_chunk_bytelen * 8 > 4096) {
        BN_GF2m_arr2poly(irr_poly_8192, field);
    } else if (needed_chunk_bytelen * 8 > 2048) {
        BN_GF2m_arr2poly(irr_poly_4096, field);
    } else {
        BN_GF2m_arr2poly(irr_poly_2048, field);
    }
    
    reed_solomon_ctx *rs_ctx = malloc(sizeof(reed_solomon_ctx));
    rs_ctx->bn_ctx = bn_ctx;
    rs_ctx->field = field;
    rs_ctx->chunk_bytelen = (BN_num_bits(field) - 1)/8;
    rs_ctx->base_data_len = 0;
    rs_ctx->base_data = NULL;
    rs_ctx->lagrange = NULL;

    return rs_ctx;
}

void reed_solomon_ctx_free(reed_solomon_ctx *rs_ctx) {
    for (uint64_t i = 0; i < rs_ctx->base_data_len; ++i) {
        BN_free(rs_ctx->base_data[i]);
        BN_free(rs_ctx->lagrange[i]);
    }
    BN_free(rs_ctx->field);
    BN_CTX_free(rs_ctx->bn_ctx);
    free(rs_ctx);
}
void set_base_data(reed_solomon_ctx *rs_ctx, uint8_t **base_data, uint64_t base_data_len) {
    BIGNUM *temp = BN_new();
    BIGNUM *temp_i = BN_new();

    rs_ctx->base_data_len = base_data_len;
    rs_ctx->base_data = calloc(base_data_len, sizeof(BIGNUM*));
    rs_ctx->lagrange = calloc(base_data_len, sizeof(BIGNUM*));

    for (uint64_t i = 0; i < base_data_len; ++i) {
        rs_ctx->base_data[i] = BN_bin2bn(base_data[i], rs_ctx->chunk_bytelen, NULL);
        
        // Set lagrange multiplier denominator (inverse once at end)
        rs_ctx->lagrange[i] = BN_new();
        BN_set_word(rs_ctx->lagrange[i], 1);
        for (uint64_t j = 0; j < base_data_len; ++j) {
            if (j == i) continue;
            BN_set_word(temp_i, i);
            BN_set_word(temp, j);
            BN_GF2m_sub(temp, temp_i, temp);
            BN_GF2m_mod_mul(rs_ctx->lagrange[i], rs_ctx->lagrange[i], temp, rs_ctx->field, rs_ctx->bn_ctx);
        }
        BN_GF2m_mod_inv(rs_ctx->lagrange[i], rs_ctx->lagrange[i], rs_ctx->field, rs_ctx->bn_ctx);
        // printBIGNUM("lagrange ", rs_ctx->lagrange[i], "\n");
        // printBIGNUM("base ", rs_ctx->base_data[i], "\n");

        assert((uint64_t) BN_num_bytes(rs_ctx->base_data[i]) <= rs_ctx->chunk_bytelen);
    }
}

void compute_chunk_at(uint64_t chunk_index, uint8_t *chunk_bytes, reed_solomon_ctx *rs_ctx)
{
    if (!chunk_bytes) return;

    BIGNUM *chunk_ind = BN_new();
    BIGNUM *chunk_val = BN_new();
    BIGNUM *lagrange = BN_new();
    BIGNUM *temp = BN_new();
    
    BN_set_word(chunk_ind, chunk_index);
    BN_set_word(chunk_val, 0);
    
    for (uint64_t i = 0; i < rs_ctx->base_data_len; ++i) {
        
        BN_GF2m_mod_mul(lagrange, rs_ctx->base_data[i], rs_ctx->lagrange[i], rs_ctx->field, rs_ctx->bn_ctx);

        for (uint64_t j = 0; j < rs_ctx->base_data_len; ++j) {
            if (j == i) continue;
            BN_set_word(temp, j);
            BN_GF2m_sub(temp, chunk_ind, temp);
            BN_GF2m_mod_mul(lagrange, lagrange, temp, rs_ctx->field, rs_ctx->bn_ctx);
        }
        
        BN_GF2m_add(chunk_val, chunk_val, lagrange);
        // printBIGNUM("lagr ", rs_ctx->lagrange[i], "\n");
        // printBIGNUM("base ", rs_ctx->base_data[i], "\n");
        // printBIGNUM("sum  ", chunk_val, "\n");
    }
    //printf("\n");

    assert((uint64_t) BN_num_bytes(chunk_val) <= rs_ctx->chunk_bytelen);
    BN_bn2binpad(chunk_val, chunk_bytes, rs_ctx->chunk_bytelen);

    BN_free(chunk_ind);
    BN_free(chunk_val);
    BN_free(temp);
    BN_free(lagrange);
}

// void printBinary(BIGNUM* num, int temp[])
// {
//     BN_GF2m_poly2arr(num, temp, sizeof(temp));
//     int j = 0;
//     while (temp[j] != -1) {
//         printf("%d ",temp[j]);
//         ++j;
//     }
//     printf("\n");
// }

int main(int argc, char* argv[]) {
    clock_t start, diff;
    double time_ms;

    uint64_t chunk_byte_len;
    uint64_t base_length;
    uint64_t num_chunks;

    if (argc < 4) {
        printf("usage: ./main <chunk_bytes> <base_length> <num_chunks>\n");
        return 0;
    }

    chunk_byte_len = strtoul(argv[1], NULL, 10);
    base_length = strtoul(argv[2], NULL, 10);
    num_chunks  = strtoul(argv[3], NULL, 10);

    // Setup the fountain context and print values

    start = clock();
    reed_solomon_ctx *rs_ctx = reed_solomon_field_ctx_new(chunk_byte_len);
    diff = clock() - start;
    time_ms = ((double) diff * 1000/ CLOCKS_PER_SEC);

    printf("generation time: %.3f ms\n", time_ms);
    printf("chunk bytes = %ld\n", rs_ctx->chunk_bytelen);
    //printf("maximal number of chunks: "); printBIGNUM("", rs_ctx->field, "\n");
    BIGNUM *num_bytes = BN_dup(rs_ctx->field);
    BN_mul_word(num_bytes, rs_ctx->chunk_bytelen);
    printf("total data bytes (including redundancies): "); printBIGNUM("", num_bytes, "\n");

    // Set random base
    uint8_t **base = calloc(base_length, sizeof(uint8_t *));
    for (uint64_t i = 0; i < base_length; ++i)
    {
        base[i] = malloc(rs_ctx->chunk_bytelen);
        RAND_bytes(base[i], rs_ctx->chunk_bytelen);
        // printf("base[%ld]=", i);
        // printBIGNUM("", temp, " ");
    }
    //printf("\n");
    
    // Set base data
    set_base_data(rs_ctx, base, base_length);

    // Test base chunks are correct
    printf("Testing correct decoding of base data...\n");
    uint8_t *chunk_bytes = malloc(rs_ctx->chunk_bytelen);
    for (uint64_t i = 0; i < base_length; ++i) {
        compute_chunk_at(i, chunk_bytes, rs_ctx);
        assert(memcmp(base[i], chunk_bytes, rs_ctx->chunk_bytelen) == 0);
    }
    printf("Done\n");

    start = clock();
    // Compute all needed chunks (also after base)
    for (uint64_t i = 0; i < num_chunks; ++i) {
        compute_chunk_at(i, chunk_bytes, rs_ctx);
        //BN_bin2bn(chunk_bytes, rs_ctx->chunk_bytelen, temp);
    }
    //printf("\n");
    diff = clock() - start;
    time_ms = ((double) diff * 1000/ CLOCKS_PER_SEC);
    
    printf("Time computing %ld chunks (%ld base bytes): %.3f ms\n", num_chunks, base_length*chunk_byte_len, time_ms);
}