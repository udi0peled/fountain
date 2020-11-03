#include <reed_solomon.h>
#include <assert.h>
#include <string.h>
#include <openssl/rand.h>
#include <time.h>

typedef struct
{ 
    int field[6];
    BN_CTX *bn_ctx;
    uint64_t chunk_bytelen;
    
    uint64_t data_len;
    uint64_t next_data_i;
    BIGNUM **data_indices; 
    BIGNUM **data;
    
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

reed_solomon_ctx *reed_solomon_ctx_new(uint64_t needed_chunk_bytelen, uint64_t data_len){

    int irr_poly_2048[] = {2048, 19, 14, 13, 0, -1};
    int irr_poly_4096[] = {4096, 27, 15, 1 , 0, -1};
    //int irr_poly_7728[] = {7728, 21, 6, 3, 0, -1};
    int irr_poly_8192[] = {8192, 9 , 5 , 2 , 0, -1};
    int irr_poly_10K[]  = {10000, 19, 13 , 9, 0, -1};

    int *field;

    if (needed_chunk_bytelen * 8 > 10000) {
        return NULL;
    } else if (needed_chunk_bytelen * 8 > 8192) {
        field = irr_poly_10K;
    } else if (needed_chunk_bytelen * 8 > 4096) {
        field = irr_poly_8192;
    } else if (needed_chunk_bytelen * 8 > 2048) {
        field = irr_poly_4096;
    } else {
        field = irr_poly_2048;
    }
    
    reed_solomon_ctx *rs_ctx = malloc(sizeof(reed_solomon_ctx));
    if (!rs_ctx) return NULL;
    
    memcpy(rs_ctx->field, field, sizeof(rs_ctx->field));

    rs_ctx->bn_ctx = BN_CTX_new();
    rs_ctx->chunk_bytelen = field[0]/8;
    rs_ctx->data_len     = data_len;
    rs_ctx->next_data_i  = 0;
    rs_ctx->data_indices = calloc(data_len, sizeof(BIGNUM*));
    rs_ctx->data         = calloc(data_len, sizeof(BIGNUM*));
    rs_ctx->lagrange     = calloc(data_len, sizeof(BIGNUM*));
    for (uint64_t i = 0; i < data_len; ++i)
    {
        rs_ctx->data[i]         = BN_CTX_get(rs_ctx->bn_ctx);
        rs_ctx->data_indices[i] = BN_CTX_get(rs_ctx->bn_ctx);
        rs_ctx->lagrange[i]     = BN_CTX_get(rs_ctx->bn_ctx);
    }

    return rs_ctx;
}

void reed_solomon_ctx_free(reed_solomon_ctx *rs_ctx) {
    for (uint64_t i = 0; i < rs_ctx->data_len; ++i) {
        BN_free(rs_ctx->data_indices[i]);
        BN_free(rs_ctx->data[i]);
        BN_free(rs_ctx->lagrange[i]);
    }
    BN_CTX_free(rs_ctx->bn_ctx);
    free(rs_ctx);
}

void set_data_at(reed_solomon_ctx *rs_ctx, uint8_t *data, uint64_t data_index) {
    uint64_t i = rs_ctx->next_data_i;

    if (i >= rs_ctx->data_len) return;

    BN_set_word(rs_ctx->data_indices[i], data_index); 
    BN_bin2bn(data, rs_ctx->chunk_bytelen, rs_ctx->data[i]);    

    BIGNUM *temp = BN_CTX_get(rs_ctx->bn_ctx);

    BN_set_word(rs_ctx->lagrange[i], 1);
    for (uint64_t j = 0; j < i; ++j) {
        BN_GF2m_sub(temp, rs_ctx->data_indices[i], rs_ctx->data_indices[j]);
        BN_GF2m_mod_mul_arr(rs_ctx->lagrange[i], rs_ctx->lagrange[i], temp, rs_ctx->field, rs_ctx->bn_ctx);
        BN_GF2m_mod_mul_arr(rs_ctx->lagrange[j], rs_ctx->lagrange[j], temp, rs_ctx->field, rs_ctx->bn_ctx);
    }
    BN_free(temp);

    rs_ctx->next_data_i += 1;
}

void set_base_data(reed_solomon_ctx *rs_ctx, uint8_t **data, uint64_t data_len) {
    for (uint64_t i = 0; i < data_len; ++i) {
        set_data_at(rs_ctx, data[i], i);
        assert((uint64_t) BN_num_bytes(rs_ctx->data[i]) <= rs_ctx->chunk_bytelen);
    }
    assert(rs_ctx->next_data_i == data_len);

    for (uint64_t i = 0; i < data_len; ++i) {
        BN_GF2m_mod_inv_arr(rs_ctx->lagrange[i], rs_ctx->lagrange[i], rs_ctx->field, rs_ctx->bn_ctx);
    }
}

void compute_chunk_at(reed_solomon_ctx *rs_ctx, uint64_t chunk_index, uint8_t *chunk_bytes)
{
    if (!chunk_bytes) return;
    if (rs_ctx->next_data_i != rs_ctx->data_len) return;

    BIGNUM *chunk_ind = BN_CTX_get(rs_ctx->bn_ctx);
    BIGNUM *chunk_val = BN_CTX_get(rs_ctx->bn_ctx);
    BIGNUM *lagrange = BN_CTX_get(rs_ctx->bn_ctx);
    BIGNUM *temp = BN_CTX_get(rs_ctx->bn_ctx);
    
    BN_set_word(chunk_ind, chunk_index);
    BN_set_word(chunk_val, 0);
    
    for (uint64_t i = 0; i < rs_ctx->data_len; ++i) {
        
        BN_GF2m_mod_mul_arr(lagrange, rs_ctx->data[i], rs_ctx->lagrange[i], rs_ctx->field, rs_ctx->bn_ctx);

        for (uint64_t j = 0; j < rs_ctx->data_len; ++j) {
            if (j == i) continue;
            BN_GF2m_sub(temp, chunk_ind, rs_ctx->data_indices[j]);
            BN_GF2m_mod_mul_arr(lagrange, lagrange, temp, rs_ctx->field, rs_ctx->bn_ctx);
        }
        
        BN_GF2m_add(chunk_val, chunk_val, lagrange);
        // printBIGNUM("lagr ", rs_ctx->lagrange[i], "\n");
        // printBIGNUM("base ", rs_ctx->data[i], "\n");
        // printBIGNUM("sum  ", chunk_val, "\n");
    }
    //printf("\n");

    //assert((uint64_t) BN_num_bytes(chunk_val) <= rs_ctx->chunk_bytelen);
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
    uint64_t num_base_chunks;

    if (argc < 3) {
        printf("usage: ./main <chunk_byte_length> <num_base_chunks>\n");
        return 0;
    }

    chunk_byte_len = strtoul(argv[1], NULL, 10);
    num_base_chunks = strtoul(argv[2], NULL, 10);

    // Setup the fountain context and print values
    reed_solomon_ctx *rs_ctx = reed_solomon_ctx_new(chunk_byte_len, num_base_chunks);
    printf("chunk bytes = %ld\n", rs_ctx->chunk_bytelen);

    // Set random base
    uint8_t **base = calloc(num_base_chunks, sizeof(uint8_t *));
    for (uint64_t i = 0; i < num_base_chunks; ++i)
    {
        base[i] = malloc(rs_ctx->chunk_bytelen);
        RAND_bytes(base[i], rs_ctx->chunk_bytelen);
    }
    
    // Set base data
    printf("set base data %ld chunks...\n", num_base_chunks);
    start = clock();

    set_base_data(rs_ctx, base, num_base_chunks);
    
    diff = clock() - start;
    time_ms = ((double) diff * 1000/ CLOCKS_PER_SEC);
    printf("Done. Time: %.3f ms\n", time_ms);

    // Test base chunks are correct
    printf("Testing correct decoding of %ld base data...\n", num_base_chunks);
    start = clock();

    uint8_t *chunk_bytes = malloc(rs_ctx->chunk_bytelen);
    for (uint64_t i = 0; i < num_base_chunks; ++i) {
        compute_chunk_at(rs_ctx, i, chunk_bytes);
        assert(memcmp(base[i], chunk_bytes, rs_ctx->chunk_bytelen) == 0);
    }

    diff = clock() - start;
    time_ms = ((double) diff * 1000/ CLOCKS_PER_SEC);
    printf("Done. Time: %.3f ms\n", time_ms);

    // Decode base from non base

    reed_solomon_ctx *decoder_ctx = reed_solomon_ctx_new(chunk_byte_len, num_base_chunks);

    // Set non-base data
    printf("Compute non-base data %ld chunks...\n", num_base_chunks);
    start = clock();

    uint64_t nonbase_shift = num_base_chunks;
    uint8_t **nonbase = calloc(num_base_chunks, sizeof(uint8_t*));
    for (uint64_t i = 0; i < num_base_chunks; ++i)
    {
        nonbase[i] = malloc(rs_ctx->chunk_bytelen);
        compute_chunk_at(rs_ctx, nonbase_shift + i, nonbase[i]);
    }

    diff = clock() - start;
    time_ms = ((double) diff * 1000/ CLOCKS_PER_SEC);
    printf("Done. Time: %.3f ms\n", time_ms);

    printf("Setting non-base data %ld chunks...\n", num_base_chunks);
    start = clock();

    for (uint64_t i = 0; i < num_base_chunks; ++i)
    {
        set_data_at(decoder_ctx, nonbase[i], nonbase_shift + i);
    }

    diff = clock() - start;
    time_ms = ((double) diff * 1000/ CLOCKS_PER_SEC);
    printf("Done. Time: %.3f ms\n", time_ms);


    printf("Decoding base from non-base data %ld chunks...\n", num_base_chunks);
    start = clock();

    uint8_t *decoded_base = malloc(chunk_byte_len);
    for (uint64_t i = 0; i < num_base_chunks; ++i)
    {
        compute_chunk_at(decoder_ctx, i, decoded_base);
        //assert(memcmp(decoded_base, base[i], chunk_byte_len));
    }
    
    assert(decoder_ctx->data_len == num_base_chunks);
    assert(decoder_ctx->data_len == decoder_ctx->next_data_i);

    diff = clock() - start;
    time_ms = ((double) diff * 1000/ CLOCKS_PER_SEC);
    printf("Done. Time: %.3f ms\n", time_ms);
}