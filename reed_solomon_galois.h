#ifndef __FOUNTAIN_REED_SOLOMON_H__
#define __FOUNTAIN_REED_SOLOMON_H__

#include <stdint.h>

typedef struct reed_solomon_encoder reed_solomon_encoder_t;
//typedef struct reed_solomon_decoder_ctx reed_solomon_decoder_ctx_t;

// Assumes data has rs_enc->data_bytelen bytes
void set_data_at(reed_solomon_encoder_t *rs_enc, uint64_t data_index, const char *data);

// Assumed computed_data has rs_enc->data_bytelen bytes allocated for result
void compute_data_at(reed_solomon_encoder_t *rs_enc, uint64_t data_index, char *computed_data);

#endif
