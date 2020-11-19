#include "reed_solomon_galois.h"
//#include "galois.h"
#include "galois_16bit_log_table.h"
#include "galois_16bit_exp_table.h"
#include "galois_16bit_inv_log_table.h"

#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <openssl/rand.h>
#include <time.h>

typedef unsigned short field_el;
typedef field_el exp_el; // 0 to 2^16-2, EXP_SIZE=2^16-1 represents 0 field element
#define EXP_SIZE 65535
#define GALOIS_BYTES 2
#define GALOIS_BITS (8*GALOIS_BYTES)

#define galois_add_field(a,b) ((a)^(b))
#define galois_mult_exp_nonzero(a,b) (((a) + (b)) % EXP_SIZE)
#define glaois_inv_exp_nonzero(a) (a == 0 ? 0 : EXP_SIZE-a)
#define galois_field_to_exp(a) (galois_16bit_log_table[a])
#define galois_exp_to_field(e) (galois_16bit_exp_table[e])

//#define galois_region_mult galois_w16_region_multiply
//#define galois_multiply galois_logtable_multiply

typedef struct reed_solomon_encoder
{ 
    uint64_t num_data;
    uint64_t num_elements;
    uint64_t data_bytelen; 
    uint64_t next_data_ind;
    field_el *data_pos;       // [num_data] positions
    exp_el  **data_val_exp;   // [num_data] regions (each data_bytelen size), each representing multiple field_elements;
    exp_el  *lagrange_w_exp;  // auxiliary for quicker interpolation
    exp_el  *temp_val_exp;   // same as above, multiples of data_val
    
} reed_solomon_encoder;

// Galois function on GALOIS_BITS=16

// Adds a (exp) to region in place, Assumes a is non-zero, but not elements of regions
void galois_add_exp_region_nonzero(exp_el *res_region, exp_el *in_region, exp_el a_exp, uint64_t num_el) {
    for (uint64_t i = 0; i < num_el; ++i) {
        // If in_region is zero, keep zero
        if (in_region[i] == EXP_SIZE) {
            res_region[i] = EXP_SIZE;
        } else {
            res_region[i] = galois_mult_exp_nonzero(a_exp, in_region[i]);
        }
    }
}

exp_el galois_add_exp(exp_el a_exp, exp_el b_exp) {
    if ((a_exp == EXP_SIZE) || (b_exp == EXP_SIZE)) return EXP_SIZE;
    return galois_mult_exp_nonzero(a_exp, b_exp);
}

void galois_region_field_to_exp(field_el *region, uint64_t num_el) {
    for (uint64_t i = 0; i < num_el; ++i) {
        region[i] = galois_field_to_exp(region[i]);
    }
}

void galois_region_exp_to_field(field_el *region, uint64_t num_el) {
    for (uint64_t i = 0; i < num_el; ++i) {
        region[i] = galois_exp_to_field(region[i]);
    }
}

// Result is a_region
void galois_add_regions(field_el *a_region, field_el *b_region, uint64_t num_el) {
    for (uint64_t i = 0; i < num_el; ++i) {
        a_region[i] = galois_add_field(a_region[i], b_region[i]);
    }
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

reed_solomon_encoder *reed_solomon_encoder_new(uint64_t num_data, uint64_t data_bytelen) {

    // TODO: Add other checks of field_el exp_el sizes
    
    if (sizeof(galois_16bit_log_table) != (EXP_SIZE+1) * sizeof(galois_16bit_log_table[0])) {
        fprintf(stderr, "reed_solomon_encoder_new -- worng log table size\n");
        return NULL;
    }

    if (sizeof(galois_16bit_exp_table) != (EXP_SIZE+1) * sizeof(galois_16bit_exp_table[0])) {
        fprintf(stderr, "reed_solomon_encoder_new -- worng exp table size\n");
        return NULL;
    }

    if (sizeof(galois_16bit_inv_log_table) != 196608 * sizeof(galois_16bit_inv_log_table[0])) {
        fprintf(stderr, "reed_solomon_encoder_new -- worng inv log table size\n");
        return NULL;
    }

    if (data_bytelen % sizeof(field_el) != 0) {
        fprintf(stderr, "reed_solomon_encoder_new -- need bytelen which is multiple of %u\n", sizeof(field_el));
        return NULL;
    }

    reed_solomon_encoder *rs_enc = malloc(sizeof(reed_solomon_encoder));
    if (!rs_enc) {
        fprintf(stderr, "reed_solomon_encoder_new -- couldn't allocate encoder\n");
        return NULL;
    }

    rs_enc->num_data       = num_data;
    rs_enc->data_bytelen   = data_bytelen;
    rs_enc->num_elements   = data_bytelen / sizeof(exp_el);
    rs_enc->next_data_ind  = 0;
    rs_enc->data_pos       = (field_el *) calloc(num_data, sizeof(field_el));
    rs_enc->lagrange_w_exp = (exp_el *)  calloc(num_data, sizeof(field_el));
    rs_enc->data_val_exp   = (exp_el **) calloc(num_data, sizeof(exp_el *));
    rs_enc->temp_val_exp   = (exp_el *) calloc(rs_enc->num_elements, sizeof(exp_el));

    if (!rs_enc->data_pos || !rs_enc->data_val_exp || !rs_enc->temp_val_exp || !rs_enc->lagrange_w_exp) {
        free(rs_enc->data_pos);
        free(rs_enc->data_val_exp);
        free(rs_enc->temp_val_exp);
        free(rs_enc->lagrange_w_exp);
        free(rs_enc);

        fprintf(stderr, "reed_solomon_encoder_new -- couldn't allocate encoder members\n");
        return NULL;
    }

    for (uint64_t i = 0; i < num_data; ++i) {
        rs_enc->data_val_exp[i] = calloc(rs_enc->num_elements, sizeof(exp_el));

        if (!rs_enc->data_val_exp[i]) {
            for (uint64_t j = 0; j <= i; ++j) {
                free(rs_enc->data_val_exp[j]);
            }
            free(rs_enc->data_pos);
            free(rs_enc->data_val_exp);
            free(rs_enc->temp_val_exp);
            free(rs_enc->lagrange_w_exp);
            free(rs_enc);
            return NULL;
        }
    }

    return rs_enc;
}

void reed_solomon_encoder_free(reed_solomon_encoder *rs_enc) {
    for (uint64_t j = 0; j < rs_enc->num_data; ++j) {
        free(rs_enc->data_val_exp[j]);
    }
    free(rs_enc->data_pos);
    free(rs_enc->data_val_exp);
    free(rs_enc->temp_val_exp);
    free(rs_enc->lagrange_w_exp);
    free(rs_enc);
}

void set_data_at(reed_solomon_encoder *rs_enc, uint64_t data_index, const char *data) {

    if (!rs_enc || !data) return;
    if (rs_enc->next_data_ind >= rs_enc->num_data){
        fprintf(stderr, "set_data_at -- trying to set data pos %lu, after full data elements needed, ignoring\n", data_index);
        return ;
    }

    field_el new_data_pos = (field_el) data_index;
    uint64_t ind = rs_enc->next_data_ind;
    for (uint64_t j = 0; j < ind; ++j) {
        if (new_data_pos == rs_enc->data_pos[j]) {
            fprintf(stderr, "set_data_at -- trying to set data pos %lu twice, ignoring\n", data_index);
            return;
        }
    }

    rs_enc->data_pos[ind] = new_data_pos;
    memcpy(rs_enc->data_val_exp[ind], data, rs_enc->data_bytelen);
    galois_region_exp_to_field((field_el *) rs_enc->data_val_exp[ind], rs_enc->num_elements);

    rs_enc->lagrange_w_exp[ind] = 0;
    exp_el temp;

    // Multiply new and i'th lagrange_w with (new_pos - pos_i) for all i < ind
    for (uint64_t j = 0; j < ind; ++j) {
        temp = galois_field_to_exp(galois_add_field(new_data_pos, rs_enc->data_pos[j]));
        if (temp == 0) {
            fprintf(stderr, "set_data_at -- trying to set data pos %lu twice, ignoring\n", data_index);
            return ;
        }
        rs_enc->lagrange_w_exp[j] = galois_add_exp(rs_enc->lagrange_w_exp[j], temp);
        rs_enc->lagrange_w_exp[ind] = galois_add_exp(rs_enc->lagrange_w_exp[ind], temp);
    }

    rs_enc->next_data_ind += 1;
}

// Computes lagrange_w interpolation at x=data_index by formula \prod(x-x_i)\cdot \sum y_i/(w_i\cdot (x-x_i))
// If x=x_i, just return value at x_i pos
void compute_data_at(reed_solomon_encoder *rs_enc, uint64_t data_index, char *computed_data)
{
    if (!rs_enc || !computed_data) return;
    if (rs_enc->next_data_ind != rs_enc->num_data) {
        fprintf(stderr, "compute_data_at -- trying to compute data pos %lu, while not all input data set, ignoring\n", data_index);
        return;
    }

    field_el data_pos = (field_el) data_index;

    // Check if x_i as one of the set data_pos, return its data_val
    for (uint64_t i = 0; i < rs_enc->num_data; ++i) {
        if (rs_enc->data_pos[i] == data_pos)
        {
            memcpy(computed_data, rs_enc->data_val_exp[i], rs_enc->data_bytelen);
            galois_region_exp_to_field((exp_el *) computed_data, rs_enc->num_elements);
            return;
        }
    }

    exp_el temp;
    exp_el lagrange_prod_exp = 0;

    memset(computed_data, 0, rs_enc->data_bytelen);
    for (uint64_t i = 0; i < rs_enc->num_data; ++i) {
        temp = galois_field_to_exp(galois_add_field(data_pos, rs_enc->data_pos[i]));
        lagrange_prod_exp = galois_add_exp(lagrange_prod_exp, temp);

        temp = galois_add_exp(temp, rs_enc->lagrange_w_exp[i]);
        temp = glaois_inv_exp_nonzero(temp);
        temp = galois_add_exp(temp, lagrange_prod_exp);

        galois_add_exp_region_nonzero(rs_enc->temp_val_exp, rs_enc->data_val_exp[i], temp, rs_enc->num_elements);
        galois_region_exp_to_field(rs_enc->temp_val_exp, rs_enc->num_elements);
        galois_add_regions(computed_data, rs_enc->temp_val_exp, rs_enc->num_elements);
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

// void time_basic(uint64_t num) {
//     clock_t start, diff;
//     double time_ms;

//     field_el a, b;
//     a = 3;
//     b = 5;

//     printf("Multiplying %lu times...\n", num);
//     start = clock();

//     for (uint64_t i = 0; i < num; ++i) {
//         b = galois_mult(a, b);
//     }

//     diff = clock() - start;
//     time_ms = ((double) diff * 1000/ CLOCKS_PER_SEC);
//     printf("Done. Time: %.3f ms\n", time_ms);


//     printf("Inverting (+1) %lu times...\n", num);
//     start = clock();
//     for (uint64_t i = 0; i < num; ++i) {
//         a = galois_inv(a);
//         a += 1;
//     }

//     diff = clock() - start;
//     time_ms = ((double) diff * 1000/ CLOCKS_PER_SEC);
//     printf("Done. Time: %.3f ms\n", time_ms);
// }

// void test_basic(uint64_t num) {

//     // Test log table:
//     uint64_t log_table_size = 65535;

//     field_el val; 
//     for (field_el i = 1; i < log_table_size; ++i) {
//         val = galois_16bit_inv_log_table[i];
//         if (galois_16bit_log_table[val] != i) {
//             printf("error at exp(%u) = %u\n", i, val);
//             printf("error at %u = log(%u)\n", galois_16bit_log_table[val], val);
//             exit(1);
//         }

//         val = galois_16bit_log_table[i];
//         if (galois_16bit_inv_log_table[val] != i) {
//             printf("error at log(%u) = %u\n", i, val);
//             printf("error at %u = exp(%u)\n", galois_16bit_inv_log_table[val], val);
//             exit(1);
//         }
//     }

//     unsigned int rand_seed = (unsigned int) time(NULL);
//     printf("Seeding randomness with %u\n", rand_seed);
//     srand(rand_seed);

//     field_el a,b,c;
//     a = rand();
//     b = rand();
//     for (uint64_t i = 0; i < num; ++i) {
//         c = galois_mult(a, b);
//         if ((a != 0) && (galois_div(c, a) != b)) {
//             fprintf(stderr, "error at %u*%u == %u\n", a, b, c);
//             fprintf(stderr, "error at %u/%u = %u \n", c, a, galois_div(c, a));
//             exit(1);
//         }
//         b = c ^ a;
//     }
// }

void test_rs(uint64_t num_data, uint64_t data_bytelen) {
    clock_t start, diff;
    double time_ms;

    unsigned int rand_seed = (unsigned int) time(NULL);
    printf("Seeding randomness with %u\n", rand_seed);
    srand(rand_seed);

    // Setup the encoder and print values
    reed_solomon_encoder *rs_enc = reed_solomon_encoder_new(num_data, data_bytelen);
    if (!rs_enc) {
        printf("Encoder initializaion Error, aborting\n");
        exit(1);
    }
    printf("Number of data values = %lu, each of byte length = %lu\n", rs_enc->num_data, rs_enc->data_bytelen);
    assert(data_bytelen == rs_enc->data_bytelen);
    assert(num_data == rs_enc->num_data);

    // Sample random base
    char **base = (char**) calloc(num_data, sizeof(char *));
    for (uint64_t i = 0; i < num_data; ++i) {
        base[i] = (char *) calloc(data_bytelen, sizeof(char));
        for (uint64_t j = 0; j < data_bytelen; ++j) base[i][j] = rand();

        // printf("base[%lu] = %u = ", i, base[i]);
        // printHexBytes("", base[i], data_bytelen, "\n", 0);
    }
    
    // Set base data
    printf("Setting base data %lu chunks...\n", num_data);
    start = clock();

    for (uint64_t i = 0; i < num_data; ++i) set_data_at(rs_enc, i, base[i]);

    diff = clock() - start;
    time_ms = ((double) diff * 1000/ CLOCKS_PER_SEC);
    printf("Done. Time: %.3f ms\n", time_ms);

    // TODO: test duplicate indices

    // Test base chunks are correct
    printf("Testing correct computing of %lu base data values...\n", num_data);
    start = clock();

    char *curr_data = malloc(rs_enc->data_bytelen);
    for (uint64_t i = 0; i < num_data; ++i) {
        compute_data_at(rs_enc, i, curr_data);
        assert(memcmp(base[i], curr_data, rs_enc->data_bytelen) == 0);
    }
    free(curr_data);

    diff = clock() - start;
    time_ms = ((double) diff * 1000/ CLOCKS_PER_SEC);
    printf("Done. Time: %.3f ms\n", time_ms);

    // Decode base from non base

    reed_solomon_encoder *rs_dec = reed_solomon_encoder_new(num_data, data_bytelen);
    if (!rs_dec) {
        printf("Decoder initializaion Error, aborting\n");
        exit(1);
    }
    assert(data_bytelen == rs_dec->data_bytelen);
    assert(num_data == rs_dec->num_data);

    // Set non-base data for decoder

    uint64_t nonbase_shift = num_data-3;
    printf("Computing non-base data %lu values [%lu, %lu)...\n", num_data, nonbase_shift, nonbase_shift+num_data);
    start = clock();

    char **nonbase = (char **) calloc(num_data, sizeof(char *));
    for (uint64_t i = 0; i < num_data; ++i)
    {
        nonbase[i] = (char *) calloc(data_bytelen, sizeof(char));
        compute_data_at(rs_enc, nonbase_shift + i, nonbase[i]);
        
        // printf("nonbase[%u] = %u = ", nonbase_shift + i, nonbase[i]);
        // printHexBytes("", nonbase[i], data_bytelen, "\n", 0);
    }

    diff = clock() - start;
    time_ms = ((double) diff * 1000/ CLOCKS_PER_SEC);
    printf("Done. Time: %.3f ms\n", time_ms);
    
    reed_solomon_encoder_free(rs_enc);

    printf("Setting non-base data %lu values [%lu, %lu) for decoder...\n", num_data, nonbase_shift, nonbase_shift+num_data);
    start = clock();

    for (uint64_t i = 0; i < num_data; ++i)
    {
        set_data_at(rs_dec, nonbase_shift + i, nonbase[i]);
    }

    diff = clock() - start;
    time_ms = ((double) diff * 1000/ CLOCKS_PER_SEC);
    printf("Done. Time: %.3f ms\n", time_ms);


    printf("Decoding (and verifying) base from non-base data %lu values...\n", num_data);
    start = clock();

    char *decoded_base = malloc(data_bytelen);
    for (uint64_t i = 0; i < num_data; ++i)
    {
        compute_data_at(rs_dec, i, decoded_base);
        assert(memcmp(decoded_base, base[i], data_bytelen) == 0);
        // printf("decoded_base[%lu] = %u = ", i, decoded_base);
        // printHexBytes("", decoded_base, data_bytelen, "\n", 0);
    }

    diff = clock() - start;
    time_ms = ((double) diff * 1000/ CLOCKS_PER_SEC);
    printf("Done. Time: %.3f ms\n", time_ms);

    // Check changing value doesnt decode correct
    *rs_dec->data_val_exp[0] += 1;
    compute_data_at(rs_dec, 0, decoded_base);
    assert(memcmp(decoded_base, base[0], data_bytelen) != 0);

    free(base);
    free(nonbase);
    free(decoded_base);
    reed_solomon_encoder_free(rs_dec);
}

// void print_chunk(const reed_solomon_encoder *rs_enc, const uint8_t *data, uint8_t index) {
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

void encode_random(uint64_t data_bytelen, uint8_t base_size, uint8_t chunk_amount) {
    data_bytelen += 0;
    base_size += 0;
    chunk_amount += 0;
    return;
    // // Setup the fountain context
    // reed_solomon_encoder *rs_enc = reed_solomon_encoder_new(data_bytelen, base_size);
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

void decode(uint64_t data_bytelen, uint8_t base_size, const uint8_t *chunks) {
    data_bytelen += 0;
    base_size += 0;
    chunks += 0;
    return;
    // reed_solomon_encoder *decoder_ctx = reed_solomon_encoder_new(data_bytelen, base_size);
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
    uint64_t data_bytelen;
    uint16_t base_size;

    if (argc < 4) usage_error();

    data_bytelen = strtoul(argv[2], NULL, 10);
    base_size = strtoul(argv[3], NULL, 10);

    if (strcmp(argv[1], "test") == 0) {
        
        //galois_create_log_tables(GALOIS_BITS);
        //time_basic(base_size);
        //time_basic(base_size*base_size);
        //test_basic(base_size);
        test_rs(base_size, data_bytelen);

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