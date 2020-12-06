#include "data_encoder_decoder.h"
#include "galois_16bit_log_table.h"
#include "galois_16bit_inv_log_table.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// This allows disabling printing by commenting out
#define __PRINT_ENABLED__
#ifdef __PRINT_ENABLED__ 
#define debug_printf(...) printf(__VA_ARGS__)
#else
#define debug_printf(...) 
#endif

typedef unsigned short field_el;
#define GALOIS_BYTES 2

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

/****************************************
 * 
 *   Reed Solomon Encoding/Decoding
 * 
 * *********************************** */ 

reed_solomon_ctx *reed_solomon_ctx_new(uint32_t num_data, uint32_t data_bytelen) {

    if (sizeof(galois_16bit_log_table) != 65536 * sizeof(galois_16bit_log_table[0])) {
        debug_printf("reed_solomon_ctx_new -- wrong log table size\n");
        return NULL;
    }

    if (sizeof(galois_16bit_inv_log_table) != 196608 * sizeof(galois_16bit_inv_log_table[0])) {
        debug_printf("reed_solomon_ctx_new -- wrong inv log table size\n");
        return NULL;
    }

    if (data_bytelen % GALOIS_BYTES != 0) {
        debug_printf("reed_solomon_ctx_new -- need bytelen which is multiple of %u\n", GALOIS_BYTES);
        return NULL;
    }

    reed_solomon_ctx *ctx = malloc(sizeof(reed_solomon_ctx));
    if (!ctx) {
        debug_printf("reed_solomon_ctx_new -- couldn't allocate encoder\n");
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

        debug_printf("reed_solomon_ctx_new -- couldn't allocate encoder members\n");
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
        //debug_printf("reed_solomon_set_data_at -- trying to set data pos %u, after full data elements needed, ignoring\n", data_index);
        return 2;
    }

    field_el new_data_pos = (field_el) data_index;
    uint32_t ind = ctx->next_data_ind;
    for (uint32_t j = 0; j < ind; ++j) {
        if (new_data_pos == ctx->data_pos[j]) {
            //debug_printf("reed_solomon_set_data_at -- trying to set data pos %u twice, ignoring\n", data_index);
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
            //debug_printf("reed_solomon_set_data_at -- trying to set data pos %u twice, ignoring\n", data_index);
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
        debug_printf("reed_solomon_compute_data_at -- trying to compute data pos %u, while not all input data set, ignoring\n", data_index);
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
        debug_printf("data_encoder_new -- invalid parameters\n");
        return NULL;
    }

    // Pad data length to be multiple of raw_chunk_bytelen (to avoid shorter last raw_chunk)
    // Pad with 4 bytes of length at the end, and all zeros before
    uint32_t zero_padding_bytelen = raw_chunk_bytelen - ((data_bytelen + sizeof(data_bytelen)) % raw_chunk_bytelen);
    uint32_t padded_bytelen = data_bytelen + zero_padding_bytelen + sizeof(data_bytelen);
    
    char *padded_data = malloc(padded_bytelen);
    if (!padded_data) {
        debug_printf("data_encoder_new -- memory allocation failure\n");
        return NULL;
    }

    memcpy(padded_data, data, data_bytelen);
    memset(padded_data + data_bytelen, 0x00, zero_padding_bytelen);
    memcpy(padded_data + data_bytelen + zero_padding_bytelen, &data_bytelen, sizeof(data_bytelen));

    // Initialize Reed Solomon ctx with relevant num of chunks and chunk size
    data_encoder_ctx *enc = malloc(sizeof(data_encoder_ctx));
    if (!enc) {
        free(padded_data);
        debug_printf("data_encoder_new -- memory allocation failure\n");
        return NULL;
    }

    uint32_t num_raw_chunks = padded_bytelen/raw_chunk_bytelen;
    enc->rs_ctx = reed_solomon_ctx_new(num_raw_chunks, raw_chunk_bytelen);
    if (!enc->rs_ctx) {
        free(padded_data);
        free(enc);
        debug_printf("data_encoder_new -- reed_solomon_ctx allocation_failure\n");
        return NULL;
    }

    enc->next_index = 0;

    char *padded_data_iter = padded_data;
    for (uint32_t i = 0; i < num_raw_chunks; ++i) {
        int res = reed_solomon_set_data_at(enc->rs_ctx, i, padded_data_iter);
        if (res != 0) {
            debug_printf("data_encoder_new -- error %d setting data index %d\n", res, i);
            free(padded_data);
            reed_solomon_ctx_free(enc->rs_ctx);
            free(enc);
            return NULL;
        }
        padded_data_iter += raw_chunk_bytelen;
    }
    
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
        debug_printf("data_decoder_new -- memory allocation failure\n");
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
            debug_printf("data_decoder_received_chunk -- memory allocation failure\n");
            return -1;
        }
    }

    if (num_chunks != dec->rs_ctx->num_data) {
        debug_printf("data_decoder_received_chunk -- got different total number of chunks\n");
        return 2;
    }

    reed_solomon_set_data_at(dec->rs_ctx, chunk_index, chunk);

    return 0;
}


uint32_t data_decoder_num_received(const data_decoder_ctx *dec)
{
    if (!dec) return -1;
    if (!dec->rs_ctx) return 0;
    return dec->rs_ctx->next_data_ind;
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
        debug_printf("data_decoder_reconstruct_data -- memory allocation failure\n");
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

/****************************************
 * 
 *  binary <-> alphanumeric conversion
 * 
 ****************************************/

char alphanumerics[] = {'0','1','2','3','4','5','6','7','8','9','A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z','+','-','*','/','%'};

char shifted_alphanumerics_invers[] = {40, 255, 255, 255, 255, 38, 36, 255, 37, 255, 39, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 255, 255, 255, 255, 255, 255, 255, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35};
#define ALPHA_SHIFT 37
#define ALPHA_LEN (sizeof(alphanumerics))
void bytes_to_alphanumeric(const char *bytes, char *alpha, uint32_t bytelen)
{
    if (bytelen % 2 != 0) {
        debug_printf("Odd bytelen [%u] while converting bytes to alphanumeric\n", bytelen);
        return;
    }

    // Take every 2 bytes and convert to 3 alphnumeric
    unsigned char *bytes_end = (unsigned char *) bytes + bytelen;
    unsigned char *uc_bytes = (unsigned char *) bytes;
    uint16_t curr_2bytes = 0;
    while (uc_bytes < bytes_end) {
        // Convert 2 bytes to 16bit integer (little endian)
        curr_2bytes = (*uc_bytes) + (*(uc_bytes+1))*256;
        uc_bytes += 2;

        // Convert 2 bytes to 3 alphanumeric (base 41) values (little endian)
        for (int i = 0; i < 3; ++i) {
            *alpha = alphanumerics[curr_2bytes % ALPHA_LEN];
            alpha += 1;
            curr_2bytes /= ALPHA_LEN;
        }
    }
}

void alphanumeric_to_bytes(const char *alpha, char *bytes, uint32_t bytelen)
{
    if (bytelen % 2 != 0) {
        debug_printf("Odd bytelen [%u] while converting alphanumeric to bytes\n", bytelen);
        return;
    }

    uint16_t factor[3] = {1, sizeof(alphanumerics), sizeof(alphanumerics)*sizeof(alphanumerics)};

    // Take every 3 alphanumeric chars and convert to 2 bytes
    unsigned char *bytes_end = (unsigned char *) bytes + bytelen;
    unsigned char *uc_bytes = (unsigned char *) bytes;
    unsigned char *uc_alpha = (unsigned char *) alpha;
    uint16_t curr_3alpha ;
    while (uc_bytes < bytes_end) {
        curr_3alpha = 0;
        for (int i = 0; i < 3; ++i) {
            if (*uc_alpha < ALPHA_SHIFT) curr_3alpha += 255;
            else curr_3alpha += shifted_alphanumerics_invers[(*uc_alpha)-ALPHA_SHIFT] * factor[i];
            uc_alpha++;
        }
        *uc_bytes = curr_3alpha % 256;
        curr_3alpha /= 256;
        *(uc_bytes+1) = curr_3alpha % 256;
        uc_bytes += 2;
    }
}