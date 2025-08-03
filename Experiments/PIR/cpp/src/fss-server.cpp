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

mpz_class evaluateEq(Fss* f, ServerKeyEq *k, uint64_t x) {
    // get num bits to be compared
#if SUMIT_MODIFICATION
    uint32_t n = f->numBits;
#else
    uint32_t n = f->numBits;
#endif

    // start at the correct LSB
    #if SUMIT_MODIFICATION
    int xi = mpz_tstbit(mpz_class(x).get_mpz_t(), (n-1));
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
            xi = mpz_tstbit(mpz_class(x).get_mpz_t(), (n-i-1));
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

// Evaluate whether x < value in function stored in key k

uint64_t evaluateLt(Fss* f, ServerKeyLt *k, uint64_t x) {

    uint32_t n = f->numBits;

    int xi = getBit(x, (64-n+1));
    unsigned char s[16];
    memcpy(s, k->s[xi], 16);
    unsigned char t = k->t[xi];
    uint64_t v = k->v[xi];

    unsigned char sArray[32];
    unsigned char temp[2];
    unsigned char out[64];
    uint64_t temp_v;
    for (uint32_t i = 1; i < n; i++) {
        if(i!=n) {
            xi = getBit(x, (64-n+i+1));
        } else {
            xi = 0;
        }
        prf(out, s, 64, f->aes_keys, f->numKeys);
        memcpy(sArray, out, 32);
        temp[0] = out[32] % 2;
        temp[1] = out[33] % 2;

        temp_v = byteArr2Int64((unsigned char*) (out + 40 + (8*xi)));
        int xStart = 16 * xi;
        memcpy(s, (unsigned char*) (sArray + xStart), 16);
        for (uint32_t j = 0; j < 16; j++) {
            s[j] = s[j] ^ k->cw[t][i-1].cs[xi][j];
        }
        //printf("%d: t: %d %d, ct: %d, bit: %d\n", i, temp[0], temp[1], k->cw[t][i-1].ct[xi], xi);
        //printf("temp_v: %lld\n", temp_v);
        v = (v + temp_v);
        v = (v + k->cw[t][i-1].cv[xi]);
        t = temp[xi] ^ k->cw[t][i-1].ct[xi];
    }
    
    return v;
}

// This function is for multi-party (3 or more parties) FSS
// for equality functions
// The API interface is similar to the 2 party version.
// One main difference is the output of the evaluation function
// is XOR homomorphic, so for additive queries like SUM and COUNT,
// the client has to add it locally.

uint32_t evaluateEqMParty(Fss *f, MPKey* key, uint32_t x)
{
    uint32_t m = 4; // Assume messages are 4 bytes long 
    uint64_t n = f->numBits;
    uint32_t p = f->numParties;
    uint32_t p2 = (uint32_t)(pow(2, p-1));
    uint64_t mu = (uint64_t)ceil((pow(2, n/2.0) * pow(2,(p-1)/2.0)));

    // sigma is last n/2 bits
    uint32_t delta = x & ((1 << (n/2)) - 1);
    uint32_t gamma = (x & (((1 << (n+1)/2) - 1) << n/2)) >> n/2;

    unsigned char** sigma = key->sigma;
    uint32_t** cw = key->cw;
    uint32_t m_bytes = m*mu;

    uint32_t* y = (uint32_t*) malloc(m_bytes);
    unsigned char* temp_out = (unsigned char*) malloc(m_bytes);
    memset(y, 0, m_bytes);
    f->numKeys = mu;
    for (int i = 0; i < p2; i++) {
        unsigned char* s = (unsigned char*)sigma[gamma] + i*16;
        bool all_zero_bytes = true;
        for (int j = 0; j < 16; j++) {
            if (s[j] != 0) {
                all_zero_bytes = false;
                break;
            }
        }
        if (!all_zero_bytes) {
            prf(temp_out, s, m_bytes, f->aes_keys, f->numKeys);
            for (int k = 0; k < mu; k++) {
                unsigned char tempIntBytes[4];
                memcpy(tempIntBytes, &temp_out[4*k], 4);
                y[k] = y[k] ^ byteArr2Int32(tempIntBytes);
            }

            for (int j = 0; j < mu; j++) {
                y[j] = cw[i][j] ^ y[j];
            }
        }
    }

    uint32_t final_ans = y[delta];
    free(y);
    free(temp_out);
    return final_ans;
}
