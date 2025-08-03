// This is the server side code for FSS which does the evaluation

#include "fss-server.h"

void initializeServer(Fss* fServer, Fss* fClient) {
    fServer->numKeys = fClient->numKeys;
    fServer->aes_keys = (AES_KEY*) malloc(sizeof(AES_KEY)*fClient->numKeys);
    memcpy(fServer->aes_keys, fClient->aes_keys, sizeof(AES_KEY)*fClient->numKeys);
    fServer->numBits = fClient->numBits;
    fServer->numParties = fClient->numParties;
    fServer->prime = fClient->prime;
}

#if SUMIT_MODIFICATION
#define IN_SIZE 48 // Value of 48 or 96 gives the best result
#include <immintrin.h>
void xor_with_key_avx(unsigned char* out, const unsigned char* key, size_t in_size) {
    size_t i = 0;
    __m256i key_vec = _mm256_loadu_si256((const __m256i*)key); // Load 32 bytes of key

    // Process in 32-byte chunks
    for (; i + 31 < in_size; i += 32) {
        __m256i out_vec = _mm256_loadu_si256((__m256i*)(out + i));
        __m256i key_repeated;
        // Repeat the 16-byte key twice to fill 32 bytes
        key_repeated = _mm256_set_m128i(_mm_loadu_si128((const __m128i*)key), _mm_loadu_si128((const __m128i*)key));
        __m256i res = _mm256_xor_si256(out_vec, key_repeated);
        _mm256_storeu_si256((__m256i*)(out + i), res);
    }
    // Handle remaining bytes
    for (; i < in_size; i++) {
        out[i] ^= key[i % 16];
    }
}

//Temporary function to check the time required for prf function
AES_KEY* prf1(unsigned char* out, unsigned char* key, uint64_t in_size, AES_KEY* aes_keys, uint32_t numKeys) {
    return aes_keys; // This is a dummy implementation for testing purposes
}
void prf2(unsigned char* out, unsigned char* key, uint64_t in_size, AES_KEY* aes_keys, uint32_t numKeys) {//Do not require to return, since it is not used
   
    // Do Matyas–Meyer–Oseas one-way compression function using different AES keys to get desired
    // output length
    uint32_t num_keys_required = in_size/16;

    for (int i = 0; i < num_keys_required; i++) {
        aesni_encrypt(key, out + (i*16), &aes_keys[i]);
    }

    xor_with_key_avx(out, key, in_size);

    return;
}
#endif
// evaluate whether x satisifies value in function stored in key k

#if SUMIT_MODIFICATION
mpz_class evaluateEq(Fss* f, ServerKeyEq *k, mpz_class x) {
#else
mpz_class evaluateEq(Fss* f, ServerKeyEq *k, uint64_t x) {
#endif
    // get num bits to be compared
#if SUMIT_MODIFICATION
    uint32_t n = f->numBits;
#else
    uint32_t n = f->numBits;
#endif

    // start at the correct LSB
    #if SUMIT_MODIFICATION
    int xi = mpz_tstbit(x.get_mpz_t(), (n-1));
    #else
    int xi = getBit(x, (64-n+1));
    #endif
    unsigned char s[16] __attribute__((aligned(16)));
    memcpy(s, k->s[xi], 16);
    unsigned char t = k->t[xi];
    
    unsigned char sArray[32];
    unsigned char temp[2];
    #if SUMIT_MODIFICATION
    unsigned char out[IN_SIZE] __attribute__((aligned(16)));
    #else
    unsigned char out[48];
    #endif

    for (uint32_t i = 1; i < n+1; i++) {
        if(i!=n) {
            #if SUMIT_MODIFICATION
            xi = mpz_tstbit(x.get_mpz_t(), (n-i-1));
            #else
            xi = getBit(x, (64-n+i+1));
            #endif
        } else {
            xi = 0;
        }
        #if SUMIT_MODIFICATION
        prf2(out, s, IN_SIZE, f->aes_keys, f->numKeys);//This function is taking too much time
        #else
        prf(out, s, 48, f->aes_keys, f->numKeys);
        #endif
        memcpy(sArray, out, 32);
        #if SUMIT_MODIFICATION
        temp[0] = out[32] & 0x01;
        temp[1] = out[33] & 0x01;
        #else
        temp[0] = out[32] % 2;
        temp[1] = out[33] % 2;
        #endif
        //printf("s: ");
        //printByteArray(s, 16);
        //printf("out: %d %d\n", out[32], out[33]);
        if (i == n) {
            break;
        }
        int xStart = 16 * xi;
        memcpy(s, (unsigned char*) (sArray + xStart), 16);
        #if SUMIT_MODIFICATION
        {
            __m256i s_vec = _mm256_loadu_si256((__m256i*)s);
            __m256i cw_vec = _mm256_loadu_si256((__m256i*)k->cw[t][i-1].cs[xi]);
            s_vec = _mm256_xor_si256(s_vec, cw_vec);
            _mm256_storeu_si256((__m256i*)s, s_vec);
        }
        #else
        for (uint32_t j = 0; j < 16; j++) {
            s[j] = s[j] ^ k->cw[t][i-1].cs[xi][j];
        }
        #endif
        //printf("After XOR: ");
        //printByteArray(s, 16);
        //printf("%d: t: %d %d, ct: %d, bit: %d\n", i, temp[0], temp[1], k->cw[t][i-1].ct[xi], xi);
        t = temp[xi] ^ k->cw[t][i-1].ct[xi];
    }

    mpz_class ans;
    unsigned char sIntArray[34];
    memcpy(sIntArray, sArray, 32);
    sIntArray[32] = temp[0];
    sIntArray[33] = temp[1];
    mpz_import(ans.get_mpz_t(), 34, 1, sizeof(sIntArray[0]), 0, 0, sIntArray);
    ans = ans * k->w;
    ans = ans % f->prime;
    return ans;
}
