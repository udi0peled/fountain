#include "reed_solomon_galois.h"
//#include "galois.h"
#include "galois_16bit_log_table.h"
#include "galois_16bit_inv_log_table.h"

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <time.h>

typedef unsigned short field_el;
#define GALOIS_BYTES 2
//#define galois_region_mult galois_w16_region_multiply
//#define galois_multiply galois_logtable_multiply

typedef struct reed_solomon_ctx
{ 
    uint32_t num_data;
    uint32_t data_bytelen; 
    uint32_t next_data_ind;
    field_el *lagrange_w;    // auxiliary for quicker interpolation
    field_el *data_pos;      // [num_data] positions
    char     **data_val;     // [num_data] regions (each data_bytelen size), each representing multiple field_elements;
    
} reed_solomon_ctx;

// Galois function on GALOIS_BITS=16

int galois_mult(int x, int y)
{
  if (x == 0 || y == 0) return 0;
  
  int sum_j = galois_16bit_log_table[x] + galois_16bit_log_table[y];
  
  return galois_16bit_inv_log_table[sum_j];
}


void galois_region_mult(char *region, int multby, int nbytes, char *r2, int add)
{
  unsigned short *ur1, *ur2, *cp;
  int prod;
  int i, log1, j, log2;
  unsigned long l, *lp2, *lptop;
  unsigned short *lp;
  int sol;

  ur1 = (unsigned short *) region;
  ur2 = (r2 == NULL) ? ur1 : (unsigned short *) r2;
  nbytes /= 2;

  if (multby == 0) {
    if (!add) {
      lp2 = (unsigned long *) ur2;
      ur2 += nbytes;
      lptop = (unsigned long *) ur2;
      while (lp2 < lptop) { *lp2 = 0; lp2++; }
    }
    return;
  }
    
  log1 = galois_16bit_log_table[multby];

  if (r2 == NULL || !add) {
    for (i = 0; i < nbytes; i++) {
      if (ur1[i] == 0) {
        ur2[i] = 0;
      } else {
        prod = galois_16bit_log_table[ur1[i]] + log1;
        ur2[i] = galois_16bit_inv_log_table[prod];
      }
    }
  } else {
    sol = sizeof(long)/2;
    lp2 = &l;
    lp = (unsigned short *) lp2;
    for (i = 0; i < nbytes; i += sol) {
      cp = ur2+i;
      lp2 = (unsigned long *) cp;
      for (j = 0; j < sol; j++) {
        if (ur1[i+j] == 0) {
          lp[j] = 0;
        } else {
          log2 = galois_16bit_log_table[ur1[i+j]];
          prod = log2 + log1;
          lp[j] = galois_16bit_inv_log_table[prod];
        }
      }
      *lp2 = (*lp2) ^ l;
    }
  }
  return; 
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

int galois_inv(int y)
{
  if (y == 0) return 0;
  return galois_div(1, y);
}

// Printing functions -- For debug

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

// Encoding/Decoding functions

reed_solomon_ctx *reed_solomon_ctx_new(uint32_t num_data, uint32_t data_bytelen) {

    if (sizeof(galois_16bit_log_table) != 65536 * sizeof(galois_16bit_log_table[0])) {
        fprintf(stderr, "reed_solomon_ctx_new -- worng log table size\n");
        return NULL;
    }

    if (sizeof(galois_16bit_inv_log_table) != 196608 * sizeof(galois_16bit_inv_log_table[0])) {
        fprintf(stderr, "reed_solomon_ctx_new -- worng inv log table size\n");
        return NULL;
    }

    if (data_bytelen % GALOIS_BYTES != 0) {
        fprintf(stderr, "reed_solomon_ctx_new -- need bytelen which is multiple of %u\n", GALOIS_BYTES);
        return NULL;
    }

    reed_solomon_ctx *ctx = malloc(sizeof(reed_solomon_ctx));
    if (!ctx) {
        fprintf(stderr, "reed_solomon_ctx_new -- couldn't allocate encoder\n");
        return NULL;
    }

    ctx->num_data      = num_data;
    ctx->data_bytelen  = data_bytelen;
    ctx->next_data_ind = 0;
    ctx->data_pos      = (field_el *) calloc(num_data, sizeof(field_el));
    ctx->data_val      = (char **)    calloc(num_data, sizeof(char*));
    ctx->lagrange_w    = (field_el *) calloc(num_data, sizeof(field_el));

    if (!ctx->data_pos || !ctx->data_val || !ctx->lagrange_w) {
        free(ctx->data_pos);
        free(ctx->data_val);
        free(ctx->lagrange_w);
        free(ctx);

        fprintf(stderr, "reed_solomon_ctx_new -- couldn't allocate encoder members\n");
        return NULL;
    }

    for (uint32_t i = 0; i < num_data; ++i) {
        ctx->data_val[i] = (char *) calloc(data_bytelen, sizeof(char));

        if (!ctx->data_val[i]) {
            for (uint32_t j = 0; j < i; ++j) free(ctx->data_val[j]);
            free(ctx->data_pos);
            free(ctx->data_val);
            free(ctx->lagrange_w);
            free(ctx);
            return NULL;
        }
    }

    return ctx;
}

void reed_solomon_ctx_free(reed_solomon_ctx *ctx) {
    for (uint32_t j = 0; j < ctx->num_data; ++j) free(ctx->data_val[j]);
    free(ctx->data_pos);
    free(ctx->data_val);
    free(ctx->lagrange_w);
    free(ctx);
}

void set_data_at(reed_solomon_ctx *ctx, uint32_t data_index, const char *data) {

    if (!ctx || !data) return;
    if (ctx->next_data_ind >= ctx->num_data){
        fprintf(stderr, "set_data_at -- trying to set data pos %u, after full data elements needed, ignoring\n", data_index);
        return ;
    }

    field_el new_data_pos = (field_el) data_index;
    uint32_t ind = ctx->next_data_ind;
    for (uint32_t j = 0; j < ind; ++j) {
        if (new_data_pos == ctx->data_pos[j]) {
            fprintf(stderr, "set_data_at -- trying to set data pos %u twice, ignoring\n", data_index);
            return;
        }
    }

    ctx->data_pos[ind] = new_data_pos;
    memcpy(ctx->data_val[ind], data, ctx->data_bytelen);  

    ctx->lagrange_w[ind] = 1;
    field_el temp;

    // Multiply new and i'th lagrange_w with (new_pos - pos_i) for all i < ind
    for (uint32_t j = 0; j < ind; ++j) {
        temp = new_data_pos ^ ctx->data_pos[j];
        if (temp == 0) {
            fprintf(stderr, "set_data_at -- trying to set data pos %u twice, ignoring\n", data_index);
            return ;
        }
        ctx->lagrange_w[j]   = galois_mult(ctx->lagrange_w[j],   temp);
        ctx->lagrange_w[ind] = galois_mult(ctx->lagrange_w[ind], temp);
    }

    ctx->next_data_ind += 1;
}

// Computes lagrange_w interpolation at x=data_index by formula \prod(x-x_i)\cdot \sum y_i/(w_i\cdot (x-x_i))
// If x=x_i, just return value at x_i pos
void compute_data_at(reed_solomon_ctx *ctx, uint32_t data_index, char *computed_data)
{
    if (!ctx || !computed_data) return;
    if (ctx->next_data_ind != ctx->num_data) {
        fprintf(stderr, "compute_data_at -- trying to compute data pos %u, while not all input data set, ignoring\n", data_index);
        return;
    }

    field_el data_pos = (field_el) data_index;

    // Check if x_i as one of the set data_pos, return its data_val
    for (uint32_t i = 0; i < ctx->num_data; ++i) {
        if (ctx->data_pos[i] == data_pos)
        {
            memcpy(computed_data, ctx->data_val[i], ctx->data_bytelen);
            return;
        }
    }

    field_el temp;
    field_el lagrange_prod = 1;
    for (uint32_t i = 0; i < ctx->num_data; ++i) {
        temp = data_pos ^ ctx->data_pos[i];
        lagrange_prod = galois_mult(lagrange_prod, temp);
    }

    memset(computed_data, 0x00, ctx->data_bytelen);
    for (uint32_t i = 0; i < ctx->num_data; ++i) {
        temp = data_pos ^ ctx->data_pos[i];
        temp = galois_mult(temp, ctx->lagrange_w[i]);
        temp = galois_div(lagrange_prod, temp);
        galois_region_mult(ctx->data_val[i], (int) temp, ctx->data_bytelen, computed_data, 1);
    }
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

void time_basic(uint32_t num) {
    clock_t start, diff;
    double time_ms;

    field_el a, b;
    a = 3;
    b = 5;

    printf("Multiplying %u times...\n", num);
    start = clock();

    for (uint32_t i = 0; i < num; ++i) {
        b = galois_mult(a, b);
    }

    diff = clock() - start;
    time_ms = ((double) diff * 1000/ CLOCKS_PER_SEC);
    printf("Done. Time: %.3f ms\n", time_ms);


    printf("Inverting (+1) %u times...\n", num);
    start = clock();
    for (uint32_t i = 0; i < num; ++i) {
        a = galois_inv(a);
        a += 1;
    }

    diff = clock() - start;
    time_ms = ((double) diff * 1000/ CLOCKS_PER_SEC);
    printf("Done. Time: %.3f ms\n", time_ms);

}

void test_basic(uint32_t num) {

    // Test log table:
    uint32_t log_table_size = 65535;

    field_el val; 
    for (field_el i = 1; i < log_table_size; ++i) {
        val = galois_16bit_inv_log_table[i];
        if (galois_16bit_log_table[val] != i) {
            printf("error at exp(%u) = %u\n", i, val);
            printf("error at %u = log(%u)\n", galois_16bit_log_table[val], val);
            exit(1);
        }

        val = galois_16bit_log_table[i];
        if (galois_16bit_inv_log_table[val] != i) {
            printf("error at log(%u) = %u\n", i, val);
            printf("error at %u = exp(%u)\n", galois_16bit_inv_log_table[val], val);
            exit(1);
        }
    }

    unsigned int rand_seed = (unsigned int) time(NULL);
    printf("Seeding randomness with %u\n", rand_seed);
    srand(rand_seed);

    field_el a,b,c;
    a = rand();
    b = rand();
    for (uint32_t i = 0; i < num; ++i) {
        c = galois_mult(a, b);
        if ((a != 0) && (galois_div(c, a) != b)) {
            fprintf(stderr, "error at %u*%u == %u\n", a, b, c);
            fprintf(stderr, "error at %u/%u = %u \n", c, a, galois_div(c, a));
            exit(1);
        }
        b = c ^ a;
    }
}
void test_rs(uint32_t num_data, uint32_t data_bytelen) {
    clock_t start, diff;
    double time_ms;

    unsigned int rand_seed = (unsigned int) time(NULL);
    printf("Seeding randomness with %u\n", rand_seed);
    srand(rand_seed);

    // Setup the encoder and print values
    reed_solomon_ctx *rs_enc = reed_solomon_ctx_new(num_data, data_bytelen);
    if (!rs_enc) {
        printf("Encoder initializaion Error, aborting\n");
        exit(1);
    }
    printf("Number of data values = %u, each of byte length = %u\n", rs_enc->num_data, rs_enc->data_bytelen);
    assert(data_bytelen == rs_enc->data_bytelen);
    assert(num_data == rs_enc->num_data);

    // Sample random base
    char **base = (char**) calloc(num_data, sizeof(char *));
    for (uint32_t i = 0; i < num_data; ++i) {
        base[i] = (char *) calloc(data_bytelen, sizeof(char));
        for (uint32_t j = 0; j < data_bytelen; ++j) base[i][j] = rand();

        // printf("base[%u] = %u = ", i, base[i]);
        // printHexBytes("", base[i], data_bytelen, "\n", 0);
    }
    
    // Set base data
    printf("Setting base data %u chunks...\n", num_data);
    start = clock();

    for (uint32_t i = 0; i < num_data; ++i) set_data_at(rs_enc, i, base[i]);

    diff = clock() - start;
    time_ms = ((double) diff * 1000/ CLOCKS_PER_SEC);
    printf("Done. Time: %.3f ms\n", time_ms);

    // TODO: test duplicate indices

    // Test base chunks are correct
    printf("Testing correct computing of %u base data values...\n", num_data);
    start = clock();

    char *curr_data = malloc(rs_enc->data_bytelen);
    for (uint32_t i = 0; i < num_data; ++i) {
        compute_data_at(rs_enc, i, curr_data);
        assert(memcmp(base[i], curr_data, rs_enc->data_bytelen) == 0);
    }
    free(curr_data);

    diff = clock() - start;
    time_ms = ((double) diff * 1000/ CLOCKS_PER_SEC);
    printf("Done. Time: %.3f ms\n", time_ms);

    // Decode base from non base

    reed_solomon_ctx *rs_dec = reed_solomon_ctx_new(num_data, data_bytelen);
    if (!rs_dec) {
        printf("Decoder initializaion Error, aborting\n");
        exit(1);
    }
    assert(data_bytelen == rs_dec->data_bytelen);
    assert(num_data == rs_dec->num_data);

    // Set non-base data for decoder

    uint32_t nonbase_shift = num_data-3;
    printf("Computing non-base data %u values [%u, %u)...\n", num_data, nonbase_shift, nonbase_shift+num_data);
    start = clock();

    char **nonbase = (char **) calloc(num_data, sizeof(char *));
    for (uint32_t i = 0; i < num_data; ++i)
    {
        nonbase[i] = (char *) calloc(data_bytelen, sizeof(char));
        compute_data_at(rs_enc, nonbase_shift + i, nonbase[i]);
        
        // printf("nonbase[%u] = %u = ", nonbase_shift + i, nonbase[i]);
        // printHexBytes("", nonbase[i], data_bytelen, "\n", 0);
    }

    diff = clock() - start;
    time_ms = ((double) diff * 1000/ CLOCKS_PER_SEC);
    printf("Done. Time: %.3f ms\n", time_ms);
    
    reed_solomon_ctx_free(rs_enc);

    printf("Setting non-base data %u values [%u, %u) for decoder...\n", num_data, nonbase_shift, nonbase_shift+num_data);
    start = clock();

    for (uint32_t i = 0; i < num_data; ++i)
    {
        set_data_at(rs_dec, nonbase_shift + i, nonbase[i]);
    }

    diff = clock() - start;
    time_ms = ((double) diff * 1000/ CLOCKS_PER_SEC);
    printf("Done. Time: %.3f ms\n", time_ms);


    printf("Decoding (and verifying) base from non-base data %u values...\n", num_data);
    start = clock();

    char *decoded_base = malloc(data_bytelen);
    for (uint32_t i = 0; i < num_data; ++i)
    {
        compute_data_at(rs_dec, i, decoded_base);
        assert(memcmp(decoded_base, base[i], data_bytelen) == 0);
        // printf("decoded_base[%u] = %u = ", i, decoded_base);
        // printHexBytes("", decoded_base, data_bytelen, "\n", 0);
    }

    diff = clock() - start;
    time_ms = ((double) diff * 1000/ CLOCKS_PER_SEC);
    printf("Done. Time: %.3f ms\n", time_ms);

    // Check changing value doesnt decode correct
    *rs_dec->data_val[0] += 1;
    compute_data_at(rs_dec, 0, decoded_base);
    assert(memcmp(decoded_base, base[0], data_bytelen) != 0);

    free(base);
    free(nonbase);
    free(decoded_base);
    reed_solomon_ctx_free(rs_dec);
}

// void print_chunk(const reed_solomon_ctx *rs_enc, const uint8_t *data, uint8_t index) {
//     printHexBytes("", &rs_enc->base_size, 1, "", 0);
//     printHexBytes("", &index, 1, "", 0);
//     printHexBytes("", data, rs_enc->data_bytelen, " ", 0);
// }

// typedef struct bignum_st {
//     BN_ULONG *d;                /* Pointer to an array of 'BN_BITS2' bit
//                                  * chunks. */
//     int top;                    /* Index of last used d +1. */
//     /* The next are internal book keeping for bn_expand. */
//     int dmax;                   /* Size of the d array. */
//     int neg;                    /* one if the number is negative */
//     int flags;
// } bignum_st;

void encode_random(uint32_t data_bytelen, uint8_t base_size, uint8_t chunk_amount) {
    data_bytelen += 0;
    base_size += 0;
    chunk_amount += 0;
    return;
    // // Setup the fountain context
    // reed_solomon_ctx *rs_enc = reed_solomon_ctx_new(data_bytelen, base_size);
    // if (!rs_enc) {
    //     printf("Initializaion Error, aborting\n");
    //     exit(1);
    // }

    // // bignum_st *s = (bignum_st *) BN_new();
    // // BN_GF2m_arr2poly(rs_enc->field, (BIGNUM*) s);
    // // printHexBytes("s = ", s->d, s->top*8, "\n", 1);

    // // Sample random base
    // uint8_t *base = calloc(base_size, data_bytelen);
    // RAND_bytes(base, base_size*data_bytelen);
    // set_base_data(rs_enc, base);

    // uint8_t *curr_data = malloc(rs_enc->data_bytelen);
    // for (uint8_t i = 0; i < chunk_amount; ++i) {
    //     compute_data_at(rs_enc, i, curr_data);
    //     print_chunk(rs_enc, curr_data, i);
    // }
    // printf("\n");

    // free(base);
    // free(curr_data);
    // reed_solomon_encoder_free(rs_enc);
}

void decode(uint32_t data_bytelen, uint8_t base_size, const uint8_t *chunks) {
    data_bytelen += 0;
    base_size += 0;
    chunks += 0;
    return;
    // reed_solomon_ctx *decoder_ctx = reed_solomon_ctx_new(data_bytelen, base_size);
    // if (!decoder_ctx) {
    //     printf("Initializaion Error, aborting\n");
    //     exit(1);
    // }

    // uint8_t chunk_index;
    // for (uint8_t i = 0; i < base_size; ++i)
    // {
    //     if (base_size != *chunks) {
    //         printf("Wrong base_size encoded in chunk %d\n", i);
    //         exit(1);
    //     }
    //     ++chunks;
    //     chunk_index = *chunks;
    //     ++chunks;
    //     set_data_at(decoder_ctx, chunks, chunk_index);
    //     chunks += data_bytelen;
    // }

    // uint8_t *decoded_base = malloc(data_bytelen);
    // for (uint8_t i = 0; i < base_size; ++i)
    // {
    //     compute_data_at(decoder_ctx, i, decoded_base);
    //     print_chunk(decoder_ctx, decoded_base, i);
    // }
    // printf("\n");

    // free(decoded_base);
    // reed_solomon_encoder_free(decoder_ctx);
}

void usage_error() {
    printf("usage: ./main test <data_bytelen> <base_size>\n");
    printf("usage: ./main encode <data_bytelen> <base_size> <amount index encoded chunks>\n");
    printf("usage: ./main decode <data_bytelen> <base_size> <spaced list of indexed chunks to decode>\n");
    exit(1);
}

int main(int argc, char* argv[]) {
    uint32_t data_bytelen;
    uint32_t base_size;

    if (argc < 4) usage_error();

    data_bytelen = strtoul(argv[2], NULL, 10);
    base_size = strtoul(argv[3], NULL, 10);

    // data_bytelen = 1024;
    // base_size = 100;

    if (strcmp(argv[1], "test") == 0) {
        
        //galois_create_log_tables(GALOIS_BITS);
        time_basic(base_size);
        time_basic(base_size*base_size);
        test_basic(base_size);
        test_rs(base_size, data_bytelen);

    } else if ((strcmp(argv[1], "encode") == 0) || (strcmp(argv[1], "enc") == 0)) {
        if (argc < 5) usage_error();

        uint32_t chunk_amount = strtoul(argv[4], NULL, 10);
        encode_random(data_bytelen, base_size, chunk_amount);
        
    } else if ((strcmp(argv[1], "decode") == 0) || (strcmp(argv[1], "dec") == 0)) {
        if ((unsigned) argc < 4+base_size) usage_error();

        uint8_t *chunks = calloc(base_size, 2+data_bytelen);
        uint8_t *curr_chunk = chunks;
        for (uint8_t i = 0; i < base_size; ++i) {
            readHexBytes(curr_chunk, data_bytelen+2, argv[4+i], strlen(argv[4+i]));
            curr_chunk += 2+data_bytelen;
        }
        decode(data_bytelen, base_size, chunks);

    } else usage_error();
}