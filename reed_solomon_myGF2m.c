#include "reed_solomon_myGF2m.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <time.h> 

typedef struct
{ 
    uint8_t base_size;
    uint8_t next_data_i;
    GF2m_el *data_indices; 
    GF2m_el *data;
    GF2m_el *lagrange;
    GF2m_extended_el field;
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

int test_myGF2m() {
    printf("Field: %d", GF2m_BITLEN);
    for (int j = 0; j < GF2m_FIELD_WEIGHT-1; ++j) printf(" %d", GF2m_FIELD[j]);
    printf("\n");

    GF2m_extended_el field;
    GF2m_get_field(field);
    printHexBytes("field =  ", (uint8_t*) field, GF2m_BYTELEN+1, "\n", 1);
    printHexBytes("field = ", (uint8_t*) field, sizeof(field), "\n", 1);

    GF2m_el a, b, c, a_inv;

    srand((unsigned int) time(NULL));
    GF2m_rand(a);
    GF2m_rand(b);

    printHexBytes("a = ", (uint8_t*) a, GF2m_BYTELEN, "\n", 1);
    printHexBytes("b = ", (uint8_t*) b, GF2m_BYTELEN, "\n", 1);

    GF2m_mul(c, a, a);
    printHexBytes("a*a = ", (uint8_t*) c, GF2m_BYTELEN, "\n", 1);
    
    GF2m_mul(c, c, a);
    printHexBytes("a*a*a = ", (uint8_t*) c, GF2m_BYTELEN, "\n", 1);

    GF2m_mul(c, a, c);
    printHexBytes("a*a*a*a = ", (uint8_t*) c, GF2m_BYTELEN, "\n", 1);

    GF2m_mul(c, a, b);
    printHexBytes("a*b = ", (uint8_t*) c, GF2m_BYTELEN, "\n", 1);

    GF2m_inv(a_inv, a, field);
    printHexBytes("a inv = ", (uint8_t*) a_inv, GF2m_BYTELEN, "\n", 1);

    GF2m_mul(a, a, a_inv);
    printHexBytes("a*(a inv) = ", (uint8_t*) a, GF2m_BYTELEN, "\n", 1);

    GF2m_mul(a, c, a_inv);
    assert(memcmp(a,b,GF2m_BYTELEN) == 0);
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

reed_solomon_ctx *reed_solomon_ctx_new(uint8_t base_size) {
    reed_solomon_ctx *rs_ctx = malloc(sizeof(reed_solomon_ctx));
    GF2m_get_field(rs_ctx->field);
    rs_ctx->base_size    = base_size;
    rs_ctx->next_data_i  = 0;
    rs_ctx->data_indices = calloc(base_size, sizeof(GF2m_el));
    rs_ctx->data         = calloc(base_size, sizeof(GF2m_el));
    rs_ctx->lagrange     = calloc(base_size, sizeof(GF2m_el));
    return rs_ctx;
}

void reed_solomon_ctx_free(reed_solomon_ctx *rs_ctx) {
    free(rs_ctx->data_indices);
    free(rs_ctx->data);
    free(rs_ctx->lagrange);
    free(rs_ctx);
}

void invert_lagrange(reed_solomon_ctx *rs_ctx) {
    for (uint8_t i = 0; i < rs_ctx->base_size; ++i) {
        GF2m_inv(rs_ctx->lagrange[i], rs_ctx->lagrange[i], rs_ctx->field);
    }
}

void set_data_at(reed_solomon_ctx *rs_ctx, const GF2m_el data, uint8_t data_index) {
    
    // TODO: ignore if data_index exists in indices
    uint8_t i = rs_ctx->next_data_i;
    
    if (i >= rs_ctx->base_size) return;

    GF2m_from_bytes(rs_ctx->data_indices[i], &data_index, 1);
    GF2m_from_bytes(rs_ctx->data[i], (uint8_t*) data, sizeof(GF2m_el));

    GF2m_el temp;
    uint8_t j = 1;
    GF2m_from_bytes(rs_ctx->lagrange[i], &j, 1);

    for (j = 0; j < i; ++j) {
        GF2m_add(temp, rs_ctx->data_indices[i], rs_ctx->data_indices[j]);
        GF2m_mul(rs_ctx->lagrange[i], rs_ctx->lagrange[i], temp);
        GF2m_mul(rs_ctx->lagrange[j], rs_ctx->lagrange[j], temp);
    }
    rs_ctx->next_data_i += 1;
    if (rs_ctx->next_data_i == rs_ctx->base_size) invert_lagrange(rs_ctx);
}

void set_base_data(reed_solomon_ctx *rs_ctx, const GF2m_el *data) {
    for (uint8_t i = 0; i < rs_ctx->base_size; ++i) {
        set_data_at(rs_ctx, data[i], i);
    }
    assert(rs_ctx->next_data_i == rs_ctx->base_size);   
}

void compute_data_at(reed_solomon_ctx *rs_ctx, uint8_t chunk_index, GF2m_el comp_data) {
    // TODO, if index is base, return quick

    if (rs_ctx->next_data_i != rs_ctx->base_size) return;

    GF2m_el chunk_ind;
    GF2m_el lagrange ;
    GF2m_el temp     ;
    
    GF2m_from_bytes(chunk_ind, &chunk_index, 1);
    GF2m_from_bytes(comp_data, NULL, 0);
    
    for (uint8_t i = 0; i < rs_ctx->base_size; ++i) {
        
        GF2m_mul(lagrange, rs_ctx->data[i], rs_ctx->lagrange[i]);

        for (uint8_t j = 0; j < rs_ctx->base_size; ++j) {
            if (j == i) continue;
            GF2m_sub(temp, chunk_ind, rs_ctx->data_indices[j]);
            GF2m_mul(lagrange, lagrange, temp);
        }
        
        GF2m_add(comp_data, comp_data, lagrange);
        // printBIGNUM("lagr ", rs_ctx->lagrange[i], "\n");
        // printBIGNUM("base ", rs_ctx->data[i], "\n");
        // printBIGNUM("sum  ", chunk_val, "\n");
    }
    //printf("\n");

    //assert((uint64_t) BN_num_bytes(chunk_val) <= rs_ctx->data_bytelen);
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

void test(uint8_t base_size) {
    clock_t start, diff;
    double time_ms;

    // Setup the fountain context and print values
    reed_solomon_ctx *rs_ctx = reed_solomon_ctx_new(base_size);
    if (!rs_ctx) {
        printf("Initializaion Error, aborting\n");
        exit(1);
    }
    printf("data byte length = %ld\n", sizeof(GF2m_el));

    // Sample random base
    GF2m_el *base = calloc(base_size, sizeof(GF2m_el));
    for (uint64_t i = 0; i < base_size; ++i) {
        GF2m_rand(base[i]);
    }
    
    // Set base data
    printf("set base data %d chunks...\n", rs_ctx->base_size);
    start = clock();

    set_base_data(rs_ctx, base);

    diff = clock() - start;
    time_ms = ((double) diff * 1000/ CLOCKS_PER_SEC);
    printf("Done. Time: %.3f ms\n", time_ms);

    // Test base chunks are correct
    printf("Testing correct decoding of %d base data...\n", base_size);
    start = clock();

    GF2m_el curr_data;
    for (uint8_t i = 0; i < base_size; ++i) {
        compute_data_at(rs_ctx, i, curr_data);
        assert(memcmp(&base[i], curr_data, sizeof(GF2m_el)) == 0);
    }

    diff = clock() - start;
    time_ms = ((double) diff * 1000/ CLOCKS_PER_SEC);
    printf("Done. Time: %.3f ms\n", time_ms);

    // Decode base from non base

    reed_solomon_ctx *decoder_ctx = reed_solomon_ctx_new(base_size);
    if (!decoder_ctx) {
        printf("Initializaion Error, aborting\n");
        exit(1);
    }

    // Set non-base data
    printf("Compute non-base data %d chunks...\n", base_size);
    start = clock();

    uint8_t nonbase_shift = base_size;
    GF2m_el *nonbase = calloc(base_size, sizeof(GF2m_el));
    for (uint8_t i = 0; i < base_size; ++i)
    {
        compute_data_at(rs_ctx, nonbase_shift + i, nonbase[i]);
    }

    diff = clock() - start;
    time_ms = ((double) diff * 1000/ CLOCKS_PER_SEC);
    printf("Done. Time: %.3f ms\n", time_ms);
    
    reed_solomon_ctx_free(rs_ctx);

    printf("Setting non-base data %d chunks...\n", base_size);
    start = clock();

    for (uint8_t i = 0; i < base_size; ++i)
    {
        set_data_at(decoder_ctx, nonbase[i], nonbase_shift + i);
    }

    diff = clock() - start;
    time_ms = ((double) diff * 1000/ CLOCKS_PER_SEC);
    printf("Done. Time: %.3f ms\n", time_ms);


    printf("Decoding base from non-base data %d chunks...\n", base_size);
    start = clock();

    GF2m_el decoded_base;
    for (uint8_t i = 0; i < base_size; ++i)
    {
        compute_data_at(decoder_ctx, i, decoded_base);
        assert(memcmp(decoded_base, base[i], sizeof(GF2m_el)) == 0);
    }
    assert(decoder_ctx->base_size == base_size);
    assert(decoder_ctx->base_size == decoder_ctx->next_data_i);

    diff = clock() - start;
    time_ms = ((double) diff * 1000/ CLOCKS_PER_SEC);
    printf("Done. Time: %.3f ms\n", time_ms);

    free(base);
    free(nonbase);
    reed_solomon_ctx_free(decoder_ctx);
}

void print_chunk(const reed_solomon_ctx *rs_ctx, const GF2m_el data, uint8_t index) {
    printHexBytes("", &rs_ctx->base_size, 1, "", 0);
    printHexBytes("", &index, 1, "", 0);
    printHexBytes("", (uint8_t*) data, sizeof(GF2m_el), " ", 0);
}

void encode_random(uint8_t base_size, uint8_t chunk_amount) {
    
    // Setup the fountain context
    reed_solomon_ctx *rs_ctx = reed_solomon_ctx_new(base_size);
    if (!rs_ctx) {
        printf("Initializaion Error, aborting\n");
        exit(1);
    }

    // Sample random base
    GF2m_el *base = calloc(base_size, sizeof(GF2m_el));
    for (uint64_t i = 0; i < base_size; ++i) GF2m_rand(base[i]);
    
    set_base_data(rs_ctx, base);

    GF2m_el curr_data;
    for (uint8_t i = 0; i < chunk_amount; ++i) {
        compute_data_at(rs_ctx, i, curr_data);
        print_chunk(rs_ctx, curr_data, i);
    }
    printf("\n");

    free(base);
    reed_solomon_ctx_free(rs_ctx);
}

void decode(uint8_t base_size, const uint8_t *chunks) {
    reed_solomon_ctx *decoder_ctx = reed_solomon_ctx_new(base_size);
    if (!decoder_ctx) {
        printf("Initializaion Error, aborting\n");
        exit(1);
    }

    GF2m_el chunk_el;
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
        GF2m_from_bytes(chunk_el, chunks, sizeof(GF2m_el));
        set_data_at(decoder_ctx, chunk_el, chunk_index);
        chunks += sizeof(GF2m_el);
    }

    GF2m_el decoded_base;
    for (uint8_t i = 0; i < base_size; ++i)
    {
        compute_data_at(decoder_ctx, i, decoded_base);
        print_chunk(decoder_ctx, decoded_base, i);
    }
    printf("\n");

    reed_solomon_ctx_free(decoder_ctx);
}

void usage_error() {
    printf("usage: ./main test   <base_size>\n");
    printf("usage: ./main encode <base_size> <amount index encoded chunks>\n");
    printf("usage: ./main decode <base_size> <spaced list of indexed chunks to decode>\n");
    exit(1);
}

int main(int argc, char* argv[]) {
    uint8_t base_size;
    
    if (argc < 3) usage_error();

    base_size = strtoul(argv[2], NULL, 10);

    if (strcmp(argv[1], "test") == 0) {

        //test_myGF2m();
        test(base_size);

    } else if ((strcmp(argv[1], "encode") == 0) || (strcmp(argv[1], "enc") == 0)) {
        if (argc < 4) usage_error();

        uint64_t chunk_amount = strtoul(argv[3], NULL, 10);
        encode_random(base_size, chunk_amount);
        
    } else if ((strcmp(argv[1], "decode") == 0) || (strcmp(argv[1], "dec") == 0)) {
        if (argc < 3+base_size) usage_error();

        uint8_t *chunks = calloc(base_size, 2 + sizeof(GF2m_el));
        uint8_t *curr_chunk = chunks;
        for (uint8_t i = 0; i < base_size; ++i) {
            readHexBytes(curr_chunk, 2+sizeof(GF2m_el), argv[3+i], strlen(argv[4+i]));
            curr_chunk += 2+sizeof(GF2m_el);
        }
        decode(base_size, chunks);

    } else usage_error();
}
