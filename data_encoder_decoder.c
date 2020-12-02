#include "data_encoder_decoder.h"
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

typedef struct data_encoder_ctx
{
    reed_solomon_ctx *rs_ctx;
    uint32_t next_index;
} data_encoder_ctx;

typedef struct data_decoder_ctx
{
    reed_solomon_ctx *rs_ctx;
    uint32_t chunk_bytelen;
} data_decoder_ctx;


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
        fprintf(stderr, "reed_solomon_ctx_new -- wrong log table size\n");
        return NULL;
    }

    if (sizeof(galois_16bit_inv_log_table) != 196608 * sizeof(galois_16bit_inv_log_table[0])) {
        fprintf(stderr, "reed_solomon_ctx_new -- wrong inv log table size\n");
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
    if (!ctx) return;
    for (uint32_t j = 0; j < ctx->num_data; ++j) free(ctx->data_val[j]);
    free(ctx->data_pos);
    free(ctx->data_val);
    free(ctx->lagrange_w);
    free(ctx);
}

// 0 - success, 1 - existing index, 2 - already enough data, -1 - other
int reed_solomon_set_data_at(reed_solomon_ctx *ctx, uint32_t data_index, const char *data) {

    if (!ctx || !data) return -1;
    if (ctx->next_data_ind >= ctx->num_data){
        //fprintf(stderr, "reed_solomon_set_data_at -- trying to set data pos %u, after full data elements needed, ignoring\n", data_index);
        return 2;
    }

    field_el new_data_pos = (field_el) data_index;
    uint32_t ind = ctx->next_data_ind;
    for (uint32_t j = 0; j < ind; ++j) {
        if (new_data_pos == ctx->data_pos[j]) {
            //fprintf(stderr, "reed_solomon_set_data_at -- trying to set data pos %u twice, ignoring\n", data_index);
            return 1;
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
            //fprintf(stderr, "reed_solomon_set_data_at -- trying to set data pos %u twice, ignoring\n", data_index);
            return 1;
        }
        ctx->lagrange_w[j]   = galois_mult(ctx->lagrange_w[j],   temp);
        ctx->lagrange_w[ind] = galois_mult(ctx->lagrange_w[ind], temp);
    }

    ctx->next_data_ind += 1;

    return 0;
}

// Computes lagrange_w interpolation at x=data_index by formula \prod(x-x_i)\cdot \sum y_i/(w_i\cdot (x-x_i))
// If x=x_i, just return value at x_i pos
void reed_solomon_compute_data_at(reed_solomon_ctx *ctx, uint32_t data_index, char *computed_data)
{
    if (!ctx || !computed_data) return;
    if (ctx->next_data_ind != ctx->num_data) {
        fprintf(stderr, "reed_solomon_compute_data_at -- trying to compute data pos %u, while not all input data set, ignoring\n", data_index);
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

/**************************** 
 * 
 *  Checksum function
 * 
 ***************************/ 

#define CRC_BYTES 4

uint32_t crc32(const char *s,size_t n) {
	uint32_t crc=0xFFFFFFFF;
	
	for(size_t i=0;i<n;i++) {
		char ch=s[i];
		for(size_t j=0;j<8;j++) {
			uint32_t b=(ch^crc)&1;
			crc>>=1;
			if(b) crc=crc^0xEDB88320;
			ch>>=1;
		}
	}
	return ~crc;
}

/**************************** 
 * 
 * Data Encoder Functions
 * 
 ***************************/ 

#define CHUNK_INDEX_BYTES 2

data_encoder_ctx *data_encoder_new(const char *data, uint32_t data_bytelen, uint32_t chunk_bytelen)
{
    uint32_t header_bytelen = CRC_BYTES + 2*CHUNK_INDEX_BYTES;
    uint32_t raw_chunk_bytelen = chunk_bytelen - header_bytelen;

    if ((!data) || (data_bytelen == 0) || (chunk_bytelen <= header_bytelen) || (chunk_bytelen % 2 != 0)) {
        printf("data_encoder_new -- invalid parameters\n");
        return NULL;
    }

    // Pad data length to be multiple of raw_chunk_bytelen (to avoid shorter last raw_chunk)
    // Pad with 4 bytes of length at the end, and all zeros before
    uint32_t zero_padding_bytelen = raw_chunk_bytelen - ((data_bytelen + sizeof(data_bytelen)) % raw_chunk_bytelen);
    uint32_t padded_bytelen = data_bytelen + zero_padding_bytelen + sizeof(data_bytelen);
    
    assert(padded_bytelen % raw_chunk_bytelen == 0);

    char *padded_data = malloc(padded_bytelen);
    if (!padded_data) {
        printf("data_encoder_new -- memory allocation failure\n");
        return NULL;
    }

    memcpy(padded_data, data, data_bytelen);
    memset(padded_data + data_bytelen, 0x00, zero_padding_bytelen);
    memcpy(padded_data + data_bytelen + zero_padding_bytelen, &data_bytelen, sizeof(data_bytelen));

    // Initialize Reed Solomon ctx with relevant num of chunks and chunk size
    data_encoder_ctx *enc = malloc(sizeof(data_encoder_ctx));
    if (!enc) {
        free(padded_data);
        printf("data_encoder_new -- memory allocation failure\n");
        return NULL;
    }

    uint32_t num_raw_chunks = padded_bytelen/raw_chunk_bytelen;
    enc->rs_ctx = reed_solomon_ctx_new(num_raw_chunks, raw_chunk_bytelen);
    if (!enc->rs_ctx) {
        free(padded_data);
        free(enc);
        printf("data_encoder_new -- reed_solomon_ctx allocation_failure\n");
        return NULL;
    }

    enc->next_index = 0;

    char *padded_data_iter = padded_data;
    for (uint32_t i = 0; i < num_raw_chunks; ++i) {
        int res = reed_solomon_set_data_at(enc->rs_ctx, i, padded_data_iter);
        if (res != 0) {
            printf("data_encoder_new -- error %d setting data index %d\n", res, i);
            free(padded_data);
            reed_solomon_ctx_free(enc->rs_ctx);
            free(enc);
            return NULL;
        }
        padded_data_iter += raw_chunk_bytelen;
    }
    assert(padded_data_iter == padded_data + padded_bytelen);

    free(padded_data);
    return enc;
}

void data_encoder_free(data_encoder_ctx *enc)
{
    if (!enc) return;
    reed_solomon_ctx_free(enc->rs_ctx);
    free(enc);
}

uint32_t data_encoder_initial_num_chunks(const data_encoder_ctx *enc)
{
    if (!enc) return 0;
    return enc->rs_ctx->num_data;
}

void data_encoder_next_chunk(data_encoder_ctx *enc, char *chunk)
{
    if (!enc) return;

    char *chunk_iter = chunk;
    reed_solomon_compute_data_at(enc->rs_ctx, enc->next_index, chunk_iter);
    chunk_iter += enc->rs_ctx->data_bytelen;
    
    memcpy(chunk_iter, &enc->next_index, CHUNK_INDEX_BYTES);
    chunk_iter += CHUNK_INDEX_BYTES;
    
    memcpy(chunk_iter, &enc->rs_ctx->num_data, CHUNK_INDEX_BYTES);
    chunk_iter += CHUNK_INDEX_BYTES;

    uint32_t crc = crc32(chunk, enc->rs_ctx->data_bytelen + 2*CHUNK_INDEX_BYTES);
    memcpy(chunk_iter, &crc, CRC_BYTES);
    chunk_iter += CRC_BYTES;

    enc->next_index++;
}

/**************************** 
 * 
 * Data Decoder Functions
 * 
 ***************************/ 

data_decoder_ctx_t *data_decoder_new(uint32_t chunk_bytelen)
{
    data_decoder_ctx *dec = malloc(sizeof(data_decoder_ctx));
    if (!dec) {
        printf("data_decoder_new -- memory allocation failure\n");
        return NULL;
    }

    dec->chunk_bytelen = chunk_bytelen;
    dec->rs_ctx = NULL; // will be initialized upon receiving first chunk
    return dec;
}

void data_decoder_free(data_decoder_ctx *dec)
{
    reed_solomon_ctx_free(dec->rs_ctx);
    free(dec);
}

static void parse_chunk(const char *chunk, uint32_t chunk_bytelen, uint32_t *index, uint32_t *total, uint32_t *crc)
{
    *index = 0;
    *total = 0;
    *crc = 0;

    char *chunk_iter = (char *) chunk + chunk_bytelen - CRC_BYTES - 2*CHUNK_INDEX_BYTES;
    
    memcpy(index, chunk_iter, CHUNK_INDEX_BYTES);
    chunk_iter += CHUNK_INDEX_BYTES;

    memcpy(total, chunk_iter, CHUNK_INDEX_BYTES);
    chunk_iter += CHUNK_INDEX_BYTES;

    memcpy(crc, chunk_iter, CRC_BYTES);
    chunk_iter += CHUNK_INDEX_BYTES;
}

int data_decoder_received_chunk(data_decoder_ctx *dec, const char *chunk)
{
    if (!dec) return -1;

    uint32_t chunk_index;
    uint32_t num_chunks;
    uint32_t crc;
    
    parse_chunk(chunk, dec->chunk_bytelen, &chunk_index, &num_chunks, &crc);

    // Compute crc here and compare
    uint32_t comp_crc = crc32(chunk, dec->chunk_bytelen - CRC_BYTES);
    if (crc != comp_crc) {
        // Wrong crc, ignore and continue
        return 1;
    }

    if (!dec->rs_ctx)
    {
        dec->rs_ctx = reed_solomon_ctx_new(num_chunks, dec->chunk_bytelen - CRC_BYTES - 2*CHUNK_INDEX_BYTES);
        
        if (!dec->rs_ctx) {
            printf("data_decoder_received_chunk -- memory allocation failure\n");
            return -1;
        }
    }

    if (num_chunks != dec->rs_ctx->num_data) {
        printf("data_decoder_received_chunk -- got different total number of chunks\n");
        return 2;
    }

    reed_solomon_set_data_at(dec->rs_ctx, chunk_index, chunk);

    return 0;
}

int data_decoder_is_finished(const data_decoder_ctx *dec)
{
    if (!dec) return -1;
    if (!dec->rs_ctx) return 0;
    if (dec->rs_ctx->next_data_ind >= dec->rs_ctx->num_data) return 1;
    return 0;
}

char *data_decoder_reconstruct_data(const data_decoder_ctx *dec, uint32_t* data_bytelen)
{
    if (data_decoder_is_finished(dec) != 1) return NULL;
    
    char *padded_data = malloc(dec->rs_ctx->num_data * dec->rs_ctx->data_bytelen);

    if (!padded_data) {
        printf("data_decoder_reconstruct_data -- memory allocation failure\n");
        return NULL;
    }

    char *padded_data_iter = padded_data;
    
    for (uint32_t i = 0; i < dec->rs_ctx->num_data; ++i)
    {
        reed_solomon_compute_data_at(dec->rs_ctx, i, padded_data_iter);
        padded_data_iter += dec->rs_ctx->data_bytelen;
    }

    // Decode 4 last bytes of padded data as real data length
    memcpy(data_bytelen, padded_data_iter - sizeof(uint32_t), sizeof(uint32_t));

    return realloc(padded_data, *data_bytelen);
}

/******************************************************************
 * 
 *  Testing Library - After this point can be deleted if used as library
 * 
 ******************************************************************/

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

    for (uint32_t i = 0; i < num_data; ++i) reed_solomon_set_data_at(rs_enc, i, base[i]);

    diff = clock() - start;
    time_ms = ((double) diff * 1000/ CLOCKS_PER_SEC);
    printf("Done. Time: %.3f ms\n", time_ms);

    // TODO: test duplicate indices

    // Test base chunks are correct
    printf("Testing correct computing of %u base data values...\n", num_data);
    start = clock();

    char *curr_data = malloc(rs_enc->data_bytelen);
    for (uint32_t i = 0; i < num_data; ++i) {
        reed_solomon_compute_data_at(rs_enc, i, curr_data);
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
        reed_solomon_compute_data_at(rs_enc, nonbase_shift + i, nonbase[i]);
        
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
        reed_solomon_set_data_at(rs_dec, nonbase_shift + i, nonbase[i]);
    }

    diff = clock() - start;
    time_ms = ((double) diff * 1000/ CLOCKS_PER_SEC);
    printf("Done. Time: %.3f ms\n", time_ms);


    printf("Decoding (and verifying) base from non-base data %u values...\n", num_data);
    start = clock();

    char *decoded_base = malloc(data_bytelen);
    for (uint32_t i = 0; i < num_data; ++i)
    {
        reed_solomon_compute_data_at(rs_dec, i, decoded_base);
        assert(memcmp(decoded_base, base[i], data_bytelen) == 0);
        // printf("decoded_base[%u] = %u = ", i, decoded_base);
        // printHexBytes("", decoded_base, data_bytelen, "\n", 0);
    }

    diff = clock() - start;
    time_ms = ((double) diff * 1000/ CLOCKS_PER_SEC);
    printf("Done. Time: %.3f ms\n", time_ms);

    // Check changing value doesnt decode correct
    *rs_dec->data_val[0] += 1;
    reed_solomon_compute_data_at(rs_dec, 0, decoded_base);
    assert(memcmp(decoded_base, base[0], data_bytelen) != 0);

    free(base);
    free(nonbase);
    free(decoded_base);
    reed_solomon_ctx_free(rs_dec);
}

void print_chunk(const char *chunk, uint32_t chunk_bytelen) {
    uint16_t total_num = 0;
    uint16_t chunk_index = 0;
    uint32_t crc = 0;

    memcpy(&chunk_index, chunk + chunk_bytelen - CRC_BYTES - 2*CHUNK_INDEX_BYTES, CHUNK_INDEX_BYTES);
    memcpy(&total_num, chunk + chunk_bytelen - CRC_BYTES - CHUNK_INDEX_BYTES, CHUNK_INDEX_BYTES);
    memcpy(&crc, chunk + chunk_bytelen - CRC_BYTES, CRC_BYTES);

    printf("chunk[%u/%u]", chunk_index, total_num);
    printHexBytes(", crc: ", (uint8_t *) &crc, sizeof(uint32_t), "\n", 0);

    uint8_t * chunk_iter = (uint8_t *) chunk;
    printHexBytes("", chunk_iter, chunk_bytelen - CRC_BYTES - 2*CHUNK_INDEX_BYTES, " ", 0);
    chunk_iter +=  chunk_bytelen - CRC_BYTES - 2*CHUNK_INDEX_BYTES;
    
    printHexBytes("", chunk_iter, CHUNK_INDEX_BYTES, " ", 0);
    chunk_iter += CHUNK_INDEX_BYTES;

    printHexBytes("", chunk_iter, CHUNK_INDEX_BYTES , " ", 0);
    chunk_iter += CHUNK_INDEX_BYTES;

    printHexBytes("", chunk_iter, CRC_BYTES, "\n", 0);
    chunk_iter += CRC_BYTES;
}

void encode_decode_random(uint32_t data_bytelen, uint32_t chunk_bytelen) {
    char *data = malloc(data_bytelen);
    time_t time_seed = time(NULL);
    printf("setting time random seed: %ld", time_seed);
    
    srand(time_seed);
    for (uint32_t i = 0; i < data_bytelen; ++i) {
        data[i] = rand() % 256;
    }

    //printHexBytes("Randomized data\n", (uint8_t*) data, data_bytelen, "\n", 1);

    data_encoder_ctx *enc = data_encoder_new(data, data_bytelen, chunk_bytelen);
    uint32_t num_chunks = data_encoder_initial_num_chunks(enc);
    uint32_t computed_num_chunks = num_chunks * 2;

    printf("Number of chunks: %d\n", num_chunks);
    printf("computing chunks number: %d\n", computed_num_chunks);

    char *chunks = calloc(computed_num_chunks, chunk_bytelen);
    char *chunk_iter = chunks;
    for (uint32_t i = 0; i < computed_num_chunks; ++i) {
        data_encoder_next_chunk(enc, chunk_iter);
        //print_chunk(chunk_iter, chunk_bytelen);
        chunk_iter += chunk_bytelen;
    }

    data_encoder_free(enc);

    // Decode 
    data_decoder_ctx_t *dec = data_decoder_new(chunk_bytelen);
    
    int finished_decoding = 0;
    
    printf("\n\nDecoding\n");
    
    uint32_t i = num_chunks;
    uint32_t count = 0;
    while (finished_decoding != 1) {
        //print_chunk(chunk_iter, chunk_bytelen);
        data_decoder_received_chunk(dec, chunks + i * chunk_bytelen);
        finished_decoding = data_decoder_is_finished(dec);
        chunk_iter += chunk_bytelen; // * (rand() % computed_num_chunks);
        i = (i + rand()) % computed_num_chunks;
        //printf("%d ", i);
        count++;
        // if (i >= 10) return;
    }
    printf("\t\t%d\n", count);

    uint32_t recon_bytelen;
    char *comp_data = data_decoder_reconstruct_data(dec, &recon_bytelen);
    data_decoder_free(dec);
    
    if (recon_bytelen != data_bytelen) {
        printf("Reconstructed data length: %u (not %u)\n", recon_bytelen, data_bytelen);
    }

    //printHexBytes("Reconstructed data\n", (uint8_t*) data, recon_bytelen, "\n", 1);
    assert(memcmp(data, comp_data, data_bytelen) == 0);
    
    free(chunks);
    free(data);
    free(comp_data);
}

void usage_error() {
    printf("usage: ./main test <data_bytelen> <base_size>\n");
    printf("usage: ./main encdec <data_bytelen> <chunk_bytelen>\n");
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
        
        time_basic(base_size);
        time_basic(base_size*base_size);
        test_basic(base_size);
        test_rs(base_size, data_bytelen);

    } else if ((strcmp(argv[1], "encdec") == 0) || (strcmp(argv[1], "enc") == 0)) {

        encode_decode_random(data_bytelen, base_size);
        
    } else usage_error();
}