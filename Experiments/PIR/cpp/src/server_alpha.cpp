#include <assert.h>
#include <stdio.h>
#include "fss/fss-common.h"
#include "fss/fss-server.h"
#include "fss/fss-client.h"
#include <immintrin.h>  // Include AVX header
#include <omp.h>
#include <vector>
#include <gmp.h>
#include <gmpxx.h>
#include <random>

#include <iterator>
#include <cstring>
#include <unistd.h>
#include "pir_common.h"


static int sock_alpha_to_beta = -1, sock_alpha_to_gamma = -1;
static char net_buf[NET_BUF_SZ] = {0};

// The value of N will determine the bitlength during the client initialization
#define NUM_TAG_BITS 3072 // 16 bits can represent up to 65536, which is more than enough for N=50000
#define B 512 // Block size in bits, can be adjusted as needed
// And number of bits determine the evalution time drastically
static mpz_class sh[N][2]; // Database to store values, each entry is a pair, {Tag, Block-content}.

// Function declarations
static int InitSrv_alpha();
static int RecvInitParamsFromBeta();
static int FinSrv_alpha();
static void TestSrv_alpha();

// Function definitions
static int InitSrv_alpha(){
    int ret = -1;
    
    // Server_alpha only connects to other servers, it does not listen
    InitConnectingSocket(SERVER_BETA_IP, BETA_LISTENING_TO_ALPHA_PORT, &sock_alpha_to_beta);
    PrintLog(LOG_LEVEL_TRACE, __FILE__, __LINE__, "Established connection with Server Beta");

    #if 0
    InitConnectingSocket(SERVER_GAMMA_IP, GAMMA_LISTENING_TO_ALPHA_PORT, &sock_alpha_to_gamma);
    PrintLog(LOG_LEVEL_TRACE, __FILE__, __LINE__, "Established connection with Server Gamma");
    #endif

    ret = RecvInitParamsFromBeta();

    if (ret != 0) {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Server Alpha: Failed to receive initialization parameters from Server Beta");
        close(sock_alpha_to_beta);
        return -1;
    } else {
        PrintLog(LOG_LEVEL_TRACE, __FILE__, __LINE__, "Server Alpha: Successfully received initialization parameters from Server Beta");
    }

    PrintLog(LOG_LEVEL_INFO, __FILE__, __LINE__, "Server Alpha initialization complete");

    return 0;
}

static int RecvInitParamsFromBeta() {
    size_t received_sz = 0;
    int ret_recv = 0;

    // Receive p
    ret_recv = recvAll(sock_alpha_to_beta, net_buf, sizeof(net_buf), &received_sz);
    if (ret_recv != 0)
    {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Failed to receive p from Server Beta");
        return -1;
    }
    p = mpz_class(std::string(net_buf, received_sz));

    // Receive q
    ret_recv = recvAll(sock_alpha_to_beta, net_buf, sizeof(net_buf), &received_sz);
    if (ret_recv != 0)
    {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Failed to receive q from Server Beta");
        return -1;
    }
    q = mpz_class(std::string(net_buf, received_sz));

    // Receive g
    ret_recv = recvAll(sock_alpha_to_beta, net_buf, sizeof(net_buf), &received_sz);
    if (ret_recv != 0)
    {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Failed to receive g from Server Beta");
        return -1;
    }
    g = mpz_class(std::string(net_buf, received_sz));

    // Receive r
    ret_recv = recvAll(sock_alpha_to_beta, net_buf, sizeof(net_buf), &received_sz);
    if (ret_recv != 0)
    {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Failed to receive r from Server Beta");
        return -1;
    }
    r = mpz_class(std::string(net_buf, received_sz));

    // Receive pk_E
    ret_recv = recvAll(sock_alpha_to_beta, net_buf, sizeof(net_buf), &received_sz);
    if (ret_recv != 0)
    {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Failed to receive pk_E from Server Beta");
        return -1;
    }
    pk_E = mpz_class(std::string(net_buf, received_sz));

    // Receive FHEcryptoContext
    ret_recv = recvAll(sock_alpha_to_beta, net_buf, sizeof(net_buf), &received_sz);
    if (ret_recv != 0)
    {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Failed to receive FHEcryptoContext from Server Beta");
        return -1;
    }
    Serial::DeserializeFromString(FHEcryptoContext, std::string(net_buf, received_sz));

    // Receive pk_F
    ret_recv = recvAll(sock_alpha_to_beta, net_buf, sizeof(net_buf), &received_sz);
    if (ret_recv != 0)
    {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Failed to receive pk_F from Server Beta");
        return -1;
    }
    Serial::DeserializeFromString(pk_F, std::string(net_buf, received_sz));

    return 0;
}

static int FinSrv_alpha(){
    int ret = -1;

    // Close the sockets
    if (sock_alpha_to_beta != -1) {
        close(sock_alpha_to_beta);
        sock_alpha_to_beta = -1;
    }
    if (sock_alpha_to_gamma != -1) {
        close(sock_alpha_to_gamma);
        sock_alpha_to_gamma = -1;
    }

    PrintLog(LOG_LEVEL_SPECIAL, __FILE__, __LINE__, "Finalized Server Alpha");

    return ret;
}

static void TestSrv_alpha(){
    int ret;
    gmp_randclass rng(gmp_randinit_default);

    //Choose random message and random exponent
    mpz_class m1 = ElGamal_randomGroupElement();
    mpz_class m2 = ElGamal_randomGroupElement();
    mpz_class exp = rng.get_z_range(q);

    mpz_class m3 = (m1*m2)%p;
    mpz_class m4;
    mpz_powm(m4.get_mpz_t(), m1.get_mpz_t(), exp.get_mpz_t(), p.get_mpz_t());

    auto [c11, c12] = ElGamal_encrypt(m1, pk_E);
    auto [c21, c22] = ElGamal_encrypt(m2, pk_E);
    auto [c31, c32] = ElGamal_mult_ct({c11, c12}, {c21, c22});
    auto [c41, c42] = ElGamal_exp_ct({c11, c12}, exp, pk_E);

    std::vector<int64_t> PIRBlockVector;
    PIRBlockVector.reserve(NUM_FHE_BLOCKS_PER_PIR_BLOCK);
    for (int64_t i = 1; i <= NUM_FHE_BLOCKS_PER_PIR_BLOCK; ++i) {
        PIRBlockVector.push_back(i);
    }
    Plaintext PIRBlockPlaintext = FHEcryptoContext->MakePackedPlaintext(PIRBlockVector);

    // The encoded vectors are encrypted
    //auto FHE_c1 = FHEcryptoContext->Encrypt(pk_F, plaintext1);
    auto FHE_c1 = FHE_encSingleMsg(PIRBlockPlaintext);

    (void)sendAll(sock_alpha_to_beta, m1.get_str().c_str(), m1.get_str().size());
    (void)sendAll(sock_alpha_to_beta, m2.get_str().c_str(), m2.get_str().size());
    (void)sendAll(sock_alpha_to_beta, m3.get_str().c_str(), m3.get_str().size());
    (void)sendAll(sock_alpha_to_beta, m4.get_str().c_str(), m4.get_str().size());
    (void)sendAll(sock_alpha_to_beta, c11.get_str().c_str(), c11.get_str().size());
    (void)sendAll(sock_alpha_to_beta, c12.get_str().c_str(), c12.get_str().size());
    (void)sendAll(sock_alpha_to_beta, c21.get_str().c_str(), c21.get_str().size());
    (void)sendAll(sock_alpha_to_beta, c22.get_str().c_str(), c22.get_str().size());
    (void)sendAll(sock_alpha_to_beta, c31.get_str().c_str(), c31.get_str().size());
    (void)sendAll(sock_alpha_to_beta, c32.get_str().c_str(), c32.get_str().size());
    (void)sendAll(sock_alpha_to_beta, c41.get_str().c_str(), c41.get_str().size());
    (void)sendAll(sock_alpha_to_beta, c42.get_str().c_str(), c42.get_str().size());
    (void)sendAll(sock_alpha_to_beta, Serial::SerializeToString(FHE_c1).c_str(), Serial::SerializeToString(FHE_c1).size());

    return;
}

int main()
{
    InitSrv_alpha();

    TestSrv_alpha();

    FinSrv_alpha();

    return 0;
}