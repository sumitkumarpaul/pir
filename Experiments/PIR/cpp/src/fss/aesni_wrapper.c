#include <stddef.h>
/* This is a dummy structure for passing the compilation, it will not cause any problem since it is not used in our execution */
unsigned int OPENSSL_ia32cap_P[4] = {0, 0, 0, 0};


extern void aesni_encrypt(const unsigned char *in, unsigned char *out, const void *key);
extern void aesni_set_encrypt_key(const unsigned char *userKey, int bits, void *key);

void aesni_encrypt_wrapper(const unsigned char *in, unsigned char *out, const void *key) {
    aesni_encrypt(in, out, key);
}

void aesni_set_encrypt_key_wrapper(const unsigned char *userKey, int bits, void *key) {
    aesni_set_encrypt_key(userKey, bits, key);
}