#include "myGF2m.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// Should be of length GF2m_FIELD_WEIGHT-1, heighest order is GF2m_BITLEN (=64*GF2m_QWORDLEN)
const int GF2m_FIELD[GF2m_FIELD_WEIGHT-1] = GF2m_FIELD_ARR;

void GF2m_from_bytes(GF2m_el el, const uint8_t *from_bytes, uint64_t from_len) {
    uint8_t *p = (uint8_t*) el;
    uint64_t copy_len = (from_len < GF2m_BYTELEN ? from_len : GF2m_BYTELEN);
    memcpy(p, from_bytes, copy_len);
    memset(p + copy_len, 0x00, GF2m_BYTELEN-copy_len);
}

void GF2m_to_bytes(const GF2m_el el, uint8_t *to_bytes, uint64_t to_len) {
    uint8_t *p = (uint8_t*) el;
    uint64_t copy_len = (to_len < GF2m_BYTELEN ? to_len : GF2m_BYTELEN);
    memcpy(to_bytes, p, copy_len);
    memset(to_bytes + copy_len, 0x00, to_len-copy_len);
}

void GF2m_copy(GF2m_el to, const GF2m_el from) {
    memcpy(to, from, GF2m_BYTELEN);
}

void GF2m_rand(GF2m_el el) {
    uint8_t *el_bytes = (uint8_t*) el;
    for (uint64_t i = 0; i < GF2m_BYTELEN; ++i) el_bytes[i] = (uint8_t) rand();
}

#define QWORD_BITS 64

void GF2m_get_field(GF2m_extended_el field) {
    uint64_t i;
    uint64_t bit;

    memset(field, 0x00, sizeof(GF2m_extended_el));

    field[GF2m_QWORDLEN] = 1;
    for (i = 0; i < GF2m_FIELD_WEIGHT-1; ++i) {
        bit = GF2m_FIELD[i];
        field[bit / QWORD_BITS] |= 1L << (bit % QWORD_BITS);
    }
}

typedef uint64_t GF2m_full_product[GF2m_QWORDLEN+GF2m_QWORDLEN+2];

// Changes a and puts mod result in r
void GF2m_mod(GF2m_el a_mod, GF2m_full_product a) {
    int j, k;
    int n, dN, d0, d1;
    
    uint64_t zz, *z;

    z = a;

    /* start reduction */
    dN = GF2m_BITLEN / QWORD_BITS;
    for (j = 2*GF2m_QWORDLEN-1; j > dN;) {
        zz = z[j];
        if (z[j] == 0) {
            j--;
            continue;
        }
        z[j] = 0;

        for (k = 0; GF2m_FIELD[k] != 0; k++) {
            /* reducing component t^p[k] */
            n = GF2m_BITLEN - GF2m_FIELD[k];
            d0 = n % QWORD_BITS;
            d1 = QWORD_BITS - d0;
            n /= QWORD_BITS;
            z[j - n] ^= (zz >> d0);
            if (d0)
                z[j - n - 1] ^= (zz << d1);
        }

        /* reducing component t^0 */
        n = dN;
        d0 = GF2m_BITLEN % QWORD_BITS;
        d1 = QWORD_BITS - d0;
        z[j - n] ^= (zz >> d0);
        if (d0)
            z[j - n - 1] ^= (zz << d1);
    }

    /* final round of reduction */
    while (j == dN) {

        d0 = GF2m_BITLEN % QWORD_BITS;
        zz = z[dN] >> d0;
        if (zz == 0)
            break;
        d1 = QWORD_BITS - d0;

        /* clear up the top d1 bits */
        if (d0)
            z[dN] = (z[dN] << d1) >> d1;
        else
            z[dN] = 0;
        z[0] ^= zz;             /* reduction t^0 component */

        for (k = 0; GF2m_FIELD[k] != 0; k++) {
            uint64_t tmp_ulong;

            /* reducing component t^p[k] */
            n = GF2m_FIELD[k] / QWORD_BITS;
            d0 = GF2m_FIELD[k] % QWORD_BITS;
            d1 = QWORD_BITS - d0;
            z[n] ^= (zz << d0);
            if (d0 && (tmp_ulong = zz >> d1))
                z[n + 1] ^= tmp_ulong;
        }

    }

    memcpy(a_mod, a, GF2m_BYTELEN);
    //memset(a, 0x00, sizeof(a));
}

static void bn_GF2m_mul_1x1(uint64_t *r1, uint64_t *r0, const uint64_t a, const uint64_t b)
{
    register uint64_t h, l, s;
    uint64_t tab[16], top3b = a >> 61;
    register uint64_t a1, a2, a4, a8;

    a1 = a & (0x1FFFFFFFFFFFFFFFULL);
    a2 = a1 << 1;
    a4 = a2 << 1;
    a8 = a4 << 1;

    tab[0] = 0;
    tab[1] = a1;
    tab[2] = a2;
    tab[3] = a1 ^ a2;
    tab[4] = a4;
    tab[5] = a1 ^ a4;
    tab[6] = a2 ^ a4;
    tab[7] = a1 ^ a2 ^ a4;
    tab[8] = a8;
    tab[9] = a1 ^ a8;
    tab[10] = a2 ^ a8;
    tab[11] = a1 ^ a2 ^ a8;
    tab[12] = a4 ^ a8;
    tab[13] = a1 ^ a4 ^ a8;
    tab[14] = a2 ^ a4 ^ a8;
    tab[15] = a1 ^ a2 ^ a4 ^ a8;

    s = tab[b & 0xF];
    l = s;
    s = tab[b >> 4 & 0xF];
    l ^= s << 4;
    h = s >> 60;
    s = tab[b >> 8 & 0xF];
    l ^= s << 8;
    h ^= s >> 56;
    s = tab[b >> 12 & 0xF];
    l ^= s << 12;
    h ^= s >> 52;
    s = tab[b >> 16 & 0xF];
    l ^= s << 16;
    h ^= s >> 48;
    s = tab[b >> 20 & 0xF];
    l ^= s << 20;
    h ^= s >> 44;
    s = tab[b >> 24 & 0xF];
    l ^= s << 24;
    h ^= s >> 40;
    s = tab[b >> 28 & 0xF];
    l ^= s << 28;
    h ^= s >> 36;
    s = tab[b >> 32 & 0xF];
    l ^= s << 32;
    h ^= s >> 32;
    s = tab[b >> 36 & 0xF];
    l ^= s << 36;
    h ^= s >> 28;
    s = tab[b >> 40 & 0xF];
    l ^= s << 40;
    h ^= s >> 24;
    s = tab[b >> 44 & 0xF];
    l ^= s << 44;
    h ^= s >> 20;
    s = tab[b >> 48 & 0xF];
    l ^= s << 48;
    h ^= s >> 16;
    s = tab[b >> 52 & 0xF];
    l ^= s << 52;
    h ^= s >> 12;
    s = tab[b >> 56 & 0xF];
    l ^= s << 56;
    h ^= s >> 8;
    s = tab[b >> 60];
    l ^= s << 60;
    h ^= s >> 4;

    /* compensate for the top three bits of a */

    if (top3b & 01) {
        l ^= b << 61;
        h ^= b >> 3;
    }
    if (top3b & 02) {
        l ^= b << 62;
        h ^= b >> 2;
    }
    if (top3b & 04) {
        l ^= b << 63;
        h ^= b >> 1;
    }

    *r1 = h;
    *r0 = l;
}

static void bn_GF2m_mul_2x2(uint64_t *r, const uint64_t a1, const uint64_t a0, const uint64_t b1, const uint64_t b0)
{
    uint64_t m1, m0;
    /* r[3] = h1, r[2] = h0; r[1] = l1; r[0] = l0 */
    bn_GF2m_mul_1x1(r + 3, r + 2, a1, b1);
    bn_GF2m_mul_1x1(r + 1, r, a0, b0);
    bn_GF2m_mul_1x1(&m1, &m0, a0 ^ a1, b0 ^ b1);
    /* Correction on m1 ^= l1 ^ h1; m0 ^= l0 ^ h0; */
    r[2] ^= m1 ^ r[1] ^ r[3];   /* h0 ^= m1 ^ l1 ^ h1; */
    r[1] = r[3] ^ r[2] ^ r[0] ^ m1 ^ m0; /* l1 ^= l0 ^ h0 ^ m0; */
}

void GF2m_mul(GF2m_el r, const GF2m_el a, const GF2m_el b) {
    int i, j, k;

    uint64_t x0, x1, y0, y1, zz[4];

    GF2m_full_product s;
    memset(s, 0x00, sizeof(GF2m_full_product));

    for (j = 0; j < GF2m_QWORDLEN; j += 2) {
        y0 = b[j];
        y1 = ((j + 1) == GF2m_QWORDLEN) ? 0 : b[j + 1];
        for (i = 0; i < GF2m_QWORDLEN; i += 2) {
            x0 = a[i];
            x1 = ((i + 1) == GF2m_QWORDLEN) ? 0 : a[i + 1];
            bn_GF2m_mul_2x2(zz, x1, x0, y1, y0);
            for (k = 0; k < 4; k++) s[i + j + k] ^= zz[k];
        }
    }

    GF2m_mod(r, s);
}

#define BN_MASK2 (0xffffffffffffffffL)

static int BN_num_bits_word(uint64_t l)
{
    uint64_t x, mask;
    int bits = (l != 0);

#if QWORD_BITS > 32
    x = l >> 32;
    mask = (0 - x) & BN_MASK2;
    mask = (0 - (mask >> (QWORD_BITS - 1)));
    bits += 32 & mask;
    l ^= (x ^ l) & mask;
#endif

    x = l >> 16;
    mask = (0 - x) & BN_MASK2;
    mask = (0 - (mask >> (QWORD_BITS - 1)));
    bits += 16 & mask;
    l ^= (x ^ l) & mask;

    x = l >> 8;
    mask = (0 - x) & BN_MASK2;
    mask = (0 - (mask >> (QWORD_BITS - 1)));
    bits += 8 & mask;
    l ^= (x ^ l) & mask;

    x = l >> 4;
    mask = (0 - x) & BN_MASK2;
    mask = (0 - (mask >> (QWORD_BITS - 1)));
    bits += 4 & mask;
    l ^= (x ^ l) & mask;

    x = l >> 2;
    mask = (0 - x) & BN_MASK2;
    mask = (0 - (mask >> (QWORD_BITS - 1)));
    bits += 2 & mask;
    l ^= (x ^ l) & mask;

    x = l >> 1;
    mask = (0 - x) & BN_MASK2;
    mask = (0 - (mask >> (QWORD_BITS - 1)));
    bits += 1 & mask;

    return bits;
}

void GF2m_inv(GF2m_el r, const GF2m_el a, const GF2m_extended_el p) {
    GF2m_extended_el b_el, c_el, u_el, v_el;

    memset(b_el, 0x00, sizeof(GF2m_extended_el));
    memset(c_el, 0x00, sizeof(GF2m_extended_el));
    memset(u_el, 0x00, sizeof(GF2m_extended_el));
    memset(v_el, 0x00, sizeof(GF2m_extended_el));

    memcpy(u_el, a, sizeof(GF2m_el));
    memcpy(v_el, p, sizeof(GF2m_extended_el));

    uint64_t *b, *c, *u, *v, *tmp;
    
    b = b_el;
    c = c_el;
    u = u_el;
    v = v_el;

    int i;
    int ubits = GF2m_BITLEN;
    int vbits = GF2m_BITLEN + 1;
    int top   = GF2m_QWORDLEN + 1;

    uint64_t *udp, *bdp, *vdp, *cdp;
    udp = u;
    bdp = b;
    b[0] = 1;
    cdp = c;
    vdp = v;

    while (1) {
        while (ubits && !(udp[0] & 1)) {
            uint64_t u0, u1, b0, b1, mask;

            u0 = udp[0];
            b0 = bdp[0];
            mask = (uint64_t)0 - (b0 & 1);
            b0 ^= p[0] & mask;
            for (i = 0; i < top - 1; i++) {
                u1 = udp[i + 1];
                udp[i] = ((u0 >> 1) | (u1 << (QWORD_BITS - 1))) & BN_MASK2;
                u0 = u1;
                b1 = bdp[i + 1] ^ (p[i + 1] & mask);
                bdp[i] = ((b0 >> 1) | (b1 << (QWORD_BITS - 1))) & BN_MASK2;
                b0 = b1;
            }
            udp[i] = u0 >> 1;
            bdp[i] = b0 >> 1;
            ubits--;
        }

        if (ubits <= QWORD_BITS) {
            if (udp[0] == 0) /* poly was reducible */
                goto err;
            if (udp[0] == 1)
                break;
        }

        if (ubits < vbits) {
            i = ubits;
            ubits = vbits;
            vbits = i;

            tmp = u;
            u = v;
            v = tmp;
            
            tmp = b;
            b = c;
            c = tmp;

            udp = vdp;
            vdp = v;
            bdp = cdp;
            cdp = c;
        }
        for (i = 0; i < top; i++) {
            udp[i] ^= vdp[i];
            bdp[i] ^= cdp[i];
        }
        if (ubits == vbits) {
            uint64_t ul;
            int utop = (ubits - 1) / QWORD_BITS;

            while ((ul = udp[utop]) == 0 && utop)
                utop--;
            ubits = utop * QWORD_BITS + BN_num_bits_word(ul);
        }
    }
    
    GF2m_copy(r, b);
    return;

err:
    memset(r, 0x00, GF2m_BYTELEN);
}

void GF2m_add(GF2m_el r, const GF2m_el a, const GF2m_el b) {
    for (uint64_t i = 0; i < GF2m_QWORDLEN; ++i)  r[i] = a[i] ^ b[i];
}