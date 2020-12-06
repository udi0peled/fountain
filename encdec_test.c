#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <time.h>

#include "data_encoder_decoder.c"

// Printing functions -- For debug

void printHexBytes(const char * prefix, const char *src, unsigned len, const char * suffix, int print_len) {
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

void readHexBytes(char *dest, uint64_t dest_len, const char* src, uint64_t src_len)
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
    printHexBytes(", crc: ", (char *) &crc, sizeof(uint32_t), "\n", 0);

    uint8_t * chunk_iter = (uint8_t *) chunk;
    printHexBytes("", (char*) chunk_iter, chunk_bytelen - CRC_BYTES - 2*CHUNK_INDEX_BYTES, " ", 0);
    chunk_iter +=  chunk_bytelen - CRC_BYTES - 2*CHUNK_INDEX_BYTES;
    
    printHexBytes("", (char*) chunk_iter, CHUNK_INDEX_BYTES, " ", 0);
    chunk_iter += CHUNK_INDEX_BYTES;

    printHexBytes("", (char*) chunk_iter, CHUNK_INDEX_BYTES , " ", 0);
    chunk_iter += CHUNK_INDEX_BYTES;

    printHexBytes("", (char*) chunk_iter, CRC_BYTES, "\n", 0);
    chunk_iter += CRC_BYTES;
}

void encode_decode_random(uint32_t data_bytelen, uint32_t chunk_bytelen) {
    char *data = malloc(data_bytelen);
    time_t time_seed = time(NULL);
    printf("setting time random seed: %ld\n", time_seed);
    
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

char *arr[2] = {
    "abc",
    "def"
};

int main(int argc, char* argv[]) {
    
    char b[100];
    char a[3*sizeof(b)/2];
    char c[sizeof(b)];

    srand(time(NULL));
    for (uint32_t i = 0; i < sizeof(b); ++i) b[i] = rand() % 256;

    printHexBytes("b = ", b, sizeof(b), "\n", 1);
    bytes_to_alphanumeric(b, a, sizeof(b));
    printHexBytes("a = ", a, sizeof(a), "\n", 1);
    alphanumeric_to_bytes(a, c, sizeof(c));
    printHexBytes("c = ", c, sizeof(c), "\n", 1);
    
    assert(memcmp(b,c,sizeof(b)) == 0);
    
    if (argc < 4) usage_error();

    uint32_t data_bytelen = strtoul(argv[2], NULL, 10);
    uint32_t base_size = strtoul(argv[3], NULL, 10);

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