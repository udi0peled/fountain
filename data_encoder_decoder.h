#ifndef __FOUNTAIN_DATA_ENCODER_DECODER_H__
#define __FOUNTAIN_DATA_ENCODER_DECODER_H__

#include <stdint.h>

typedef struct data_encoder_ctx data_encoder_ctx_t;
typedef struct data_decoder_ctx data_decoder_ctx_t;

/**************************** 
 * 
 * Data Encoder Functions
 * 
 ***************************/ 

// chunk_bytelen must be an even number
data_encoder_ctx_t *data_encoder_new(const char *data, uint32_t data_bytelen, uint32_t chunk_bytelen);
void data_encoder_free(data_encoder_ctx_t *enc);

// Minimal number of chunks needed to reconstruct the data later with the decoder
uint32_t data_encoder_initial_num_chunks(const data_encoder_ctx_t *enc);

// Repeatedly generte new chunks (of size chunk_bytelen) from data
// Might (and should) be called MORE then initial num of chunks (up to 65535)
void data_encoder_next_chunk(data_encoder_ctx_t *enc, char *chunk);

/**************************** 
 * 
 * Data Decoder Functions
 * 
 ***************************/ 

data_decoder_ctx_t *data_decoder_new(uint32_t chunk_bytelen);
void data_decoder_free(data_decoder_ctx_t *dec);

// Receive chunks, must be of chunk_bytelen size
// Returns 0 if ok, 1 if wrong checksum crc (continue without problem), 2 if differenct total number then previous, -1 other error
int data_decoder_received_chunk(data_decoder_ctx_t *dec, const char *chunk);

// Boolean if enough chunks were received to reconstruct the original data
int data_decoder_is_finished(const data_decoder_ctx_t *dec);

// Call only if is_finished is true, returns data and sets its byte length
// The returned value should be freed externally
char *data_decoder_reconstruct_data(const data_decoder_ctx_t *dec, uint32_t* data_bytelen);

#endif