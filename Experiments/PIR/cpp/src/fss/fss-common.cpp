#include "fss-common.h"

#ifdef AESNI
int aesni_set_encrypt_key(const unsigned char *userKey, int bits, AES_KEY *key) { // TODO: Dummy function, implement later
    #if 1 //Currently this is calling the standard AES_set_encrypt_key, but it should be replaced with an optimized version of aesni_set_encrypt_key
    return AES_set_encrypt_key(userKey, bits, key);
    #else
    #endif 
}
void aesni_encrypt(const unsigned char *in, unsigned char *out, const AES_KEY *key) { // TODO: Dummy function, implement later 
    #if 1 //Currently this is calling the standard AES_encrypt, but it should be replaced with an optimized version of aesni_encrypt
    AES_encrypt(in, out, key);
    #else
    #endif
    return;
}
#endif

AES_KEY* prf(unsigned char* out, unsigned char* key, uint64_t in_size, AES_KEY* aes_keys, uint32_t numKeys) {
#ifndef AESNI
    // check if there is aes-ni instruction
    uint32_t eax, ebx, ecx, edx;

    eax = ebx = ecx = edx = 0;
    __get_cpuid(1, &eax, &ebx, &ecx, &edx);
#endif
    
    AES_KEY* temp_keys = aes_keys;
    // Do Matyas–Meyer–Oseas one-way compression function using different AES keys to get desired
    // output length
    uint32_t num_keys_required = in_size/16;
    if (num_keys_required > numKeys) {
        free(temp_keys);
        temp_keys = (AES_KEY*) malloc(sizeof(AES_KEY)*num_keys_required); 
        for (int i = 0; i < num_keys_required; i++) {
            unsigned char rand_bytes[16];
            if (!RAND_bytes(rand_bytes, 16)) {
                printf("Random bytes failed.\n");
            }
#ifndef AESNI
            if ((ecx & bit_AES) > 0) {
                aesni_set_encrypt_key(rand_bytes, 128, &(temp_keys[i]));
            } else {
                AES_set_encrypt_key(rand_bytes, 128, &(temp_keys[i]));
            }
#else
            aesni_set_encrypt_key(rand_bytes, 128, &(temp_keys[i]));
#endif
        }
    }
    for (int i = 0; i < num_keys_required; i++) {
#ifndef AESNI
        if ((ecx & bit_AES) > 0) {
            aesni_encrypt(key, out + (i*16), &temp_keys[i]);
        } else {
            AES_encrypt(key, out + (i*16), &temp_keys[i]);
        }
#else
        aesni_encrypt(key, out + (i*16), &temp_keys[i]);
#endif
    }
    for (int i = 0; i < in_size; i++) {
        out[i] = out[i] ^ key[i%16];
    }
    return temp_keys;
}
