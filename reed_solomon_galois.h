#ifndef __FOUNTAIN_REED_SOLOMON_H__
#define __FOUNTAIN_REED_SOLOMON_H__

#include <stdint.h>

typedef struct reed_solomon_ctx reed_solomon_ctx_t;
//typedef struct reed_solomon_decoder_ctx reed_solomon_decoder_ctx_t;

reed_solomon_ctx_t *reed_solomon_ctx_new(uint32_t num_data, uint32_t data_bytelen);
void reed_solomon_ctx_free(reed_solomon_ctx_t *rs_enc);

// Assumes data has rs_enc->data_bytelen bytes
void set_data_at(reed_solomon_ctx_t *rs_enc, uint32_t data_index, const char *data);

// Assumed computed_data has rs_enc->data_bytelen bytes allocated for result
void compute_data_at(reed_solomon_ctx_t *rs_enc, uint32_t data_index, char *computed_data);

void test_rs(uint32_t num_data, uint32_t data_bytelen);

//ctx *data_encoder_new(const char *data, uint32_t data_bytelen, uint32_t chunk_bytelen);
//uint32_t data_encoder_num_chunks(const ctx*)
//void data_encoder_at_index(const ctx *encoder, uint32_t index, char *chunk)
//void data_encoder_free(ctx *encoder)

//ctx *data_decoder_new(uint32_t chunk_bytelen);
//void data_decoder_received_chunk(ctx *decoder, const char *chunk)
//int data_decoder_is_finished(const *decoder)
//uint32_t data_decoder_get_data(const ctx *decoder, char **data)
//void data_encoder_free(ctx *decoder)

#endif
