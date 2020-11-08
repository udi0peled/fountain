#ifndef __FOUNTAIN_MY_GF2m_H__
#define __FOUNTAIN_MY_GF2m_H__

#include <stdint.h>

#define __GALOIS_FIELD_8192__

#ifdef __GALOIS_FIELD_9984__
#define GF2m_QWORDLEN 156
#define GF2m_FIELD_ARR {27, 10, 7, 0}

#elif defined __GALOIS_FIELD_8192__
#define GF2m_QWORDLEN 128
#define GF2m_FIELD_ARR {9 , 5 , 2 , 0}

#elif defined __GALOIS_FIELD_4096__
#define GF2m_QWORDLEN 64
#define GF2m_FIELD_ARR {27, 15, 1 , 0}

#elif defined __GALOIS_FIELD_2048__
#define GF2m_QWORDLEN 32
#define GF2m_FIELD_ARR {19, 14, 13, 0}

#else
#define GF2m_QWORDLEN 1
#define GF2m_FIELD_ARR {4, 3, 1, 0}

#endif


#define GF2m_BYTELEN (GF2m_QWORDLEN*8)
#define GF2m_BITLEN (GF2m_BYTELEN*8)
#define GF2m_FIELD_WEIGHT 5

// Powers in irreducible poly, must decrease, and first (GF2m_BITLEN) doesn't appear
extern const int GF2m_FIELD[GF2m_FIELD_WEIGHT-1];

typedef uint64_t GF2m_el[GF2m_QWORDLEN];
typedef uint64_t GF2m_extended_el[GF2m_QWORDLEN+1];

// Copies at most GF2m bytes needed from_bytes (if less, pad with 0x00)
void GF2m_from_bytes(GF2m_el el, const uint8_t *from_bytes, uint64_t from_len);
void GF2m_to_bytes(const GF2m_el el, uint8_t *to_bytes, uint64_t to_len);
void GF2m_copy(GF2m_el to, const GF2m_el from);
void GF2m_get_field(GF2m_extended_el field);
void GF2m_rand(GF2m_el el);

void GF2m_mul(GF2m_el r, const GF2m_el a, const GF2m_el b);
void GF2m_inv(GF2m_el r, const GF2m_el a, const GF2m_extended_el mod);
void GF2m_add(GF2m_el r, const GF2m_el a, const GF2m_el b);
#define GF2m_sub(r, a, b) GF2m_add(r, a, b);

#endif