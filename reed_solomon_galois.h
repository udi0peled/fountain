#ifndef __FOUNTAIN_REED_SOLOMON_H__
#define __FOUNTAIN_REED_SOLOMON_H__

#include <stdint.h>

typedef struct reed_solomon_ctx reed_solomon_ctx_t;
//typedef struct reed_solomon_decoder_ctx reed_solomon_decoder_ctx_t;

reed_solomon_ctx_t *reed_solomon_encoder_new(uint32_t num_data, uint32_t data_bytelen);
void reed_solomon_encoder_free(reed_solomon_ctx_t *rs_enc);

// Assumes data has rs_enc->data_bytelen bytes
void set_data_at(reed_solomon_ctx_t *rs_enc, uint32_t data_index, const char *data);

// Assumed computed_data has rs_enc->data_bytelen bytes allocated for result
void compute_data_at(reed_solomon_ctx_t *rs_enc, uint32_t data_index, char *computed_data);

void test_rs(uint32_t num_data, uint32_t data_bytelen);

#endif
