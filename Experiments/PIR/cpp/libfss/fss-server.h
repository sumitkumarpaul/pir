#ifndef FSS_SERVER_H
#define FSS_SERVER_H

#include "fss-common.h"
#include <cmath>

// Initializes server with information from the client, namely aes_keys for PRF and numBits in input domain
void initializeServer(Fss* fServer, Fss* fClient);

#if SUMIT_MODIFICATION
mpz_class evaluateEq(Fss* f, ServerKeyEq *k, mpz_class x);
#else
// Runs point(delta) FSS given key on input x for 2 parties/providers
mpz_class evaluateEq(Fss* f, ServerKeyEq *k, uint64_t x);
#endif

#endif
