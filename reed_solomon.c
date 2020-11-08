#include "reed_solomon.h"
#include <assert.h>
#include <string.h>
#include <openssl/rand.h>
#include <time.h>

typedef struct
{ 
    BN_CTX *bn_ctx;
    int *field;

    uint8_t base_size;
    uint64_t data_bytelen;
    uint8_t next_data_i;
    BIGNUM **data_indices; 
    BIGNUM **bn_data;
    
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

void readHexBytes(uint8_t* dest, uint64_t dest_len, const char* src, uint64_t src_len)
{
    uint64_t i = 0;
    uint64_t j = 0;
    unsigned int del;
    while (i < dest_len && j < src_len) {
        del = 0;
        sscanf(src + j,"%2hhx%n", dest + i, &del);
        j += del;
        i++;
    }
    //if (j%2 == 1) fprintf("Warning: odd hex length for \"%s\"\n", src);
}

void printBIGNUM(const char * prefix, const BIGNUM *bn, const char * suffix) {
  char *bn_str = BN_bn2dec(bn);
  printf("%s%s%s", prefix, bn_str, suffix);
  free(bn_str);
}

reed_solomon_ctx *reed_solomon_ctx_new(uint64_t base_data_bytelen, uint8_t base_size){
    BN_CTX *bn_ctx = BN_CTX_new();

    int irr_poly_64[]   = {64, 4 ,3, 1, 0, -1};
    int irr_poly_80[]   = {80, 9, 4, 2, 0, -1};
    int irr_poly_2048[] = {2048, 19, 14, 13, 0, -1};
    int irr_poly_4096[] = {4096, 27, 15, 1 , 0, -1};
    int irr_poly_8192[] = {8192, 9 , 5 , 2 , 0, -1};
    int irr_poly_10K[]  = {10000, 19, 13 , 9, 0, -1};

    int *field = malloc(6 * sizeof(int));

    if (base_data_bytelen * 8 > 10000) {
        memcpy(field, irr_poly_10K, sizeof(irr_poly_10K));
    } else if (base_data_bytelen * 8 == 8192) {
        memcpy(field, irr_poly_8192, sizeof(irr_poly_8192));
    } else if (base_data_bytelen * 8 == 4096) {
        memcpy(field, irr_poly_4096, sizeof(irr_poly_4096));
    } else if (base_data_bytelen * 8 == 2048) {
        memcpy(field, irr_poly_2048, sizeof(irr_poly_2048));
    } else if (base_data_bytelen * 8 == 80) {
        memcpy(field, irr_poly_80, sizeof(irr_poly_80));
    } else if (base_data_bytelen * 8 == 64) {
        memcpy(field, irr_poly_64, sizeof(irr_poly_80));
    } else {
        BN_CTX_free(bn_ctx);
        return NULL;
    }

    BN_CTX_start(bn_ctx);

    reed_solomon_ctx *rs_ctx = malloc(sizeof(reed_solomon_ctx));
    rs_ctx->bn_ctx = bn_ctx;
    rs_ctx->field = field;
    rs_ctx->data_bytelen = base_data_bytelen;
    
    rs_ctx->base_size    = base_size;
    rs_ctx->next_data_i  = 0;
    rs_ctx->data_indices = calloc(base_size, sizeof(BIGNUM*));
    rs_ctx->bn_data      = calloc(base_size, sizeof(BIGNUM*));
    rs_ctx->lagrange     = calloc(base_size, sizeof(BIGNUM*));

    for (uint8_t i = 0; i < base_size; ++i)
    {
        rs_ctx->bn_data[i]      = BN_CTX_get(bn_ctx);
        rs_ctx->data_indices[i] = BN_CTX_get(bn_ctx);
        rs_ctx->lagrange[i]     = BN_CTX_get(bn_ctx);
    }

    return rs_ctx;
}

void reed_solomon_ctx_free(reed_solomon_ctx *rs_ctx) {
    BN_CTX_end(rs_ctx->bn_ctx);
    BN_CTX_free(rs_ctx->bn_ctx);
    free(rs_ctx->data_indices);
    free(rs_ctx->bn_data);
    free(rs_ctx->lagrange);
    free(rs_ctx->field);
    free(rs_ctx);
}

void invert_lagrange(reed_solomon_ctx *rs_ctx) {
    for (uint8_t i = 0; i < rs_ctx->base_size; ++i) {
        BN_GF2m_mod_inv_arr(rs_ctx->lagrange[i], rs_ctx->lagrange[i], rs_ctx->field, rs_ctx->bn_ctx);
    }
}

void set_data_at(reed_solomon_ctx *rs_ctx, const uint8_t *data, uint8_t data_index) {
    uint8_t i = rs_ctx->next_data_i;

    if (i >= rs_ctx->base_size) return;

    BN_set_word(rs_ctx->data_indices[i], data_index);
    BN_bin2bn(data, rs_ctx->data_bytelen, rs_ctx->bn_data[i]);    

    BIGNUM *temp = BN_CTX_get(rs_ctx->bn_ctx);

    BN_set_word(rs_ctx->lagrange[i], 1);
    for (uint8_t j = 0; j < i; ++j) {
        BN_GF2m_sub(temp, rs_ctx->data_indices[i], rs_ctx->data_indices[j]);
        BN_GF2m_mod_mul_arr(rs_ctx->lagrange[i], rs_ctx->lagrange[i], temp, rs_ctx->field, rs_ctx->bn_ctx);
        BN_GF2m_mod_mul_arr(rs_ctx->lagrange[j], rs_ctx->lagrange[j], temp, rs_ctx->field, rs_ctx->bn_ctx);
    }

    rs_ctx->next_data_i += 1;
    
    if (rs_ctx->next_data_i == rs_ctx->base_size) invert_lagrange(rs_ctx);
}

void set_base_data(reed_solomon_ctx *rs_ctx, uint8_t *data) {
    uint8_t *curr_data = data;
    for (uint8_t i = 0; i < rs_ctx->base_size; ++i) {
        set_data_at(rs_ctx, curr_data, i);
        curr_data += rs_ctx->data_bytelen;
        assert((uint64_t) BN_num_bytes(rs_ctx->bn_data[i]) <= rs_ctx->data_bytelen);
    }
    assert(rs_ctx->next_data_i == rs_ctx->base_size);
}

void compute_data_at(reed_solomon_ctx *rs_ctx, uint8_t chunk_index, uint8_t *comp_data)
{
    if (!comp_data) return;
    if (rs_ctx->next_data_i != rs_ctx->base_size) return;

    BIGNUM *chunk_ind = BN_CTX_get(rs_ctx->bn_ctx);
    BIGNUM *chunk_val = BN_CTX_get(rs_ctx->bn_ctx);
    BIGNUM *lagrange  = BN_CTX_get(rs_ctx->bn_ctx);
    BIGNUM *temp      = BN_CTX_get(rs_ctx->bn_ctx);
    
    BN_set_word(chunk_ind, chunk_index);
    BN_set_word(chunk_val, 0);
    
    for (uint8_t i = 0; i < rs_ctx->base_size; ++i) {
        
        BN_GF2m_mod_mul_arr(lagrange, rs_ctx->bn_data[i], rs_ctx->lagrange[i], rs_ctx->field, rs_ctx->bn_ctx);

        for (uint8_t j = 0; j < rs_ctx->base_size; ++j) {
            if (j == i) continue;
            BN_GF2m_sub(temp, chunk_ind, rs_ctx->data_indices[j]);
            BN_GF2m_mod_mul_arr(lagrange, lagrange, temp, rs_ctx->field, rs_ctx->bn_ctx);
        }
        
        BN_GF2m_add(chunk_val, chunk_val, lagrange);
        // printBIGNUM("lagr ", rs_ctx->lagrange[i], "\n");
        // printBIGNUM("base ", rs_ctx->bn_data[i], "\n");
        // printBIGNUM("sum  ", chunk_val, "\n");
    }
    //printf("\n");

    //assert((uint64_t) BN_num_bytes(chunk_val) <= rs_ctx->data_bytelen);
    BN_bn2binpad(chunk_val, comp_data, rs_ctx->data_bytelen);
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



void time_GF2m(uint64_t reps, uint64_t data_bytelen) {
    clock_t start, diff;
    double time_ms;

    reed_solomon_ctx *rs_ctx = reed_solomon_ctx_new(data_bytelen, 1);
    if (!rs_ctx) {
        printf("Initializaion Error, aborting\n");
        exit(1);
    }
    assert(rs_ctx->data_bytelen == data_bytelen);
    printf("data byte length = %ld\n", rs_ctx->data_bytelen);

    BIGNUM *field = BN_CTX_get(rs_ctx->bn_ctx);
    BN_GF2m_arr2poly(rs_ctx->field, field);
    for (int j = 0; rs_ctx->field[j] != -1; ++j) printf("%d ", rs_ctx->field[j]);
    printf("\n");
    
    BIGNUM *a = BN_CTX_get(rs_ctx->bn_ctx);
    BIGNUM *b = BN_CTX_get(rs_ctx->bn_ctx); 
    BIGNUM *a_inv = BN_CTX_get(rs_ctx->bn_ctx); 

    uint8_t *bytes = malloc(data_bytelen);
    
    RAND_bytes(bytes, data_bytelen);
    BN_bin2bn(bytes, data_bytelen, a);

    RAND_bytes(bytes, data_bytelen);
    BN_bin2bn(bytes, data_bytelen, b);

    start = clock();

    for (uint64_t i = 0; i < reps; ++i) BN_GF2m_mod_mul_arr(a, a, a, rs_ctx->field, rs_ctx->bn_ctx);
    
    diff = clock() - start;
    time_ms = ((double) diff * 1000/ CLOCKS_PER_SEC);
    printf("Done. (a*a)x%lu, Time: %.3f ms\n", reps, time_ms);

    start = clock();

    for (uint64_t i = 0; i < reps; ++i){
        BN_GF2m_mod_inv(a_inv, a, field, rs_ctx->bn_ctx);
        BN_GF2m_add(a, a, b);
    }

    diff = clock() - start;
    time_ms = ((double) diff * 1000/ CLOCKS_PER_SEC);
    printf("Done. (a_inv + b)x%lu, Time: %.3f ms\n", reps, time_ms);

    reed_solomon_ctx_free(rs_ctx);
}

void test(uint64_t data_bytelen, uint8_t base_size) {
    clock_t start, diff;
    double time_ms;

    // Setup the fountain context and print values
    reed_solomon_ctx *rs_ctx = reed_solomon_ctx_new(data_bytelen, base_size);
    if (!rs_ctx) {
        printf("Initializaion Error, aborting\n");
        exit(1);
    }
    assert(rs_ctx->data_bytelen == data_bytelen);
    printf("data byte length = %ld\n", rs_ctx->data_bytelen);

    // Sample random base
    uint8_t *base = calloc(base_size, data_bytelen);
    RAND_bytes(base, base_size * rs_ctx->data_bytelen);
    
    // Set base data
    printf("set base data %d chunks...\n", base_size);
    start = clock();

    set_base_data(rs_ctx, base);

    diff = clock() - start;
    time_ms = ((double) diff * 1000/ CLOCKS_PER_SEC);
    printf("Done. Time: %.3f ms\n", time_ms);

    // Test base chunks are correct
    printf("Testing correct decoding of %d base data...\n", base_size);
    start = clock();

    uint8_t *curr_data = malloc(rs_ctx->data_bytelen);
    for (uint8_t i = 0; i < base_size; ++i) {
        compute_data_at(rs_ctx, i, curr_data);
        assert(memcmp(base + i*data_bytelen, curr_data, rs_ctx->data_bytelen) == 0);
    }
    free(curr_data);

    diff = clock() - start;
    time_ms = ((double) diff * 1000/ CLOCKS_PER_SEC);
    printf("Done. Time: %.3f ms\n", time_ms);

    // Decode base from non base

    reed_solomon_ctx *decoder_ctx = reed_solomon_ctx_new(data_bytelen, base_size);
    if (!decoder_ctx) {
        printf("Initializaion Error, aborting\n");
        exit(1);
    }

    // Set non-base data
    printf("Compute non-base data %d chunks...\n", base_size);
    start = clock();

    uint8_t nonbase_shift = base_size;
    uint8_t *nonbase = calloc(base_size, data_bytelen);
    for (uint8_t i = 0; i < base_size; ++i)
    {
        compute_data_at(rs_ctx, nonbase_shift + i, nonbase + i*data_bytelen);
    }

    diff = clock() - start;
    time_ms = ((double) diff * 1000/ CLOCKS_PER_SEC);
    printf("Done. Time: %.3f ms\n", time_ms);
    
    reed_solomon_ctx_free(rs_ctx);

    printf("Setting non-base data %d chunks...\n", base_size);
    start = clock();

    for (uint8_t i = 0; i < base_size; ++i)
    {
        set_data_at(decoder_ctx, nonbase + i*data_bytelen, nonbase_shift + i);
    }

    diff = clock() - start;
    time_ms = ((double) diff * 1000/ CLOCKS_PER_SEC);
    printf("Done. Time: %.3f ms\n", time_ms);


    printf("Decoding base from non-base data %d chunks...\n", base_size);
    start = clock();

    uint8_t *decoded_base = malloc(data_bytelen);
    for (uint8_t i = 0; i < base_size; ++i)
    {
        compute_data_at(decoder_ctx, i, decoded_base);
        assert(memcmp(decoded_base, base + i * data_bytelen, data_bytelen) == 0);
    }
    free(decoded_base);
    
    assert(decoder_ctx->base_size == base_size);
    assert(decoder_ctx->base_size == decoder_ctx->next_data_i);

    diff = clock() - start;
    time_ms = ((double) diff * 1000/ CLOCKS_PER_SEC);
    printf("Done. Time: %.3f ms\n", time_ms);

    free(base);
    free(nonbase);
    reed_solomon_ctx_free(decoder_ctx);
}

void print_chunk(const reed_solomon_ctx *rs_ctx, const uint8_t *data, uint8_t index) {
    printHexBytes("", &rs_ctx->base_size, 1, "", 0);
    printHexBytes("", &index, 1, "", 0);
    printHexBytes("", data, rs_ctx->data_bytelen, " ", 0);
}

// typedef struct bignum_st {
//     BN_ULONG *d;                /* Pointer to an array of 'BN_BITS2' bit
//                                  * chunks. */
//     int top;                    /* Index of last used d +1. */
//     /* The next are internal book keeping for bn_expand. */
//     int dmax;                   /* Size of the d array. */
//     int neg;                    /* one if the number is negative */
//     int flags;
// } bignum_st;

void encode_random(uint64_t data_bytelen, uint8_t base_size, uint8_t chunk_amount) {
    
    // Setup the fountain context
    reed_solomon_ctx *rs_ctx = reed_solomon_ctx_new(data_bytelen, base_size);
    if (!rs_ctx) {
        printf("Initializaion Error, aborting\n");
        exit(1);
    }

    // bignum_st *s = (bignum_st *) BN_new();
    // BN_GF2m_arr2poly(rs_ctx->field, (BIGNUM*) s);
    // printHexBytes("s = ", s->d, s->top*8, "\n", 1);

    // Sample random base
    uint8_t *base = calloc(base_size, data_bytelen);
    RAND_bytes(base, base_size*data_bytelen);
    set_base_data(rs_ctx, base);

    uint8_t *curr_data = malloc(rs_ctx->data_bytelen);
    for (uint8_t i = 0; i < chunk_amount; ++i) {
        compute_data_at(rs_ctx, i, curr_data);
        print_chunk(rs_ctx, curr_data, i);
    }
    printf("\n");

    free(base);
    free(curr_data);
    reed_solomon_ctx_free(rs_ctx);
}

void decode(uint64_t data_bytelen, uint8_t base_size, const uint8_t *chunks) {
    reed_solomon_ctx *decoder_ctx = reed_solomon_ctx_new(data_bytelen, base_size);
    if (!decoder_ctx) {
        printf("Initializaion Error, aborting\n");
        exit(1);
    }

    uint8_t chunk_index;
    for (uint8_t i = 0; i < base_size; ++i)
    {
        if (base_size != *chunks) {
            printf("Wrong base_size encoded in chunk %d\n", i);
            exit(1);
        }
        ++chunks;
        chunk_index = *chunks;
        ++chunks;
        set_data_at(decoder_ctx, chunks, chunk_index);
        chunks += data_bytelen;
    }

    uint8_t *decoded_base = malloc(data_bytelen);
    for (uint8_t i = 0; i < base_size; ++i)
    {
        compute_data_at(decoder_ctx, i, decoded_base);
        print_chunk(decoder_ctx, decoded_base, i);
    }
    printf("\n");

    free(decoded_base);
    reed_solomon_ctx_free(decoder_ctx);
}

void usage_error() {
    printf("usage: ./main test <data_bytelen> <base_size>\n");
    printf("usage: ./main encode <data_bytelen> <base_size> <amount index encoded chunks>\n");
    printf("usage: ./main decode <data_bytelen> <base_size> <spaced list of indexed chunks to decode>\n");
    exit(1);
}


int main(int argc, char* argv[]) {
    uint64_t data_bytelen;
    uint8_t base_size;

    if (argc < 4) usage_error();

    data_bytelen = strtoul(argv[2], NULL, 10);
    base_size = strtoul(argv[3], NULL, 10);

    if (strcmp(argv[1], "test") == 0) {

        time_GF2m(1000, data_bytelen);
        test(data_bytelen, base_size);

    } else if ((strcmp(argv[1], "encode") == 0) || (strcmp(argv[1], "enc") == 0)) {
        if (argc < 5) usage_error();

        uint64_t chunk_amount = strtoul(argv[4], NULL, 10);
        encode_random(data_bytelen, base_size, chunk_amount);
        
    } else if ((strcmp(argv[1], "decode") == 0) || (strcmp(argv[1], "dec") == 0)) {
        if (argc < 4+base_size) usage_error();

        uint8_t *chunks = calloc(base_size, 2+data_bytelen);
        uint8_t *curr_chunk = chunks;
        for (uint8_t i = 0; i < base_size; ++i) {
            readHexBytes(curr_chunk, data_bytelen+2, argv[4+i], strlen(argv[4+i]));
            curr_chunk += 2+data_bytelen;
        }
        decode(data_bytelen, base_size, chunks);

    } else usage_error();
}