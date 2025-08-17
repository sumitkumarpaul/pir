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
    int valread = recv(sock_alpha_to_beta, net_buf, sizeof(net_buf), 0);

    if (valread <= 0) {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Server Alpha: Read failed, closing connection with Server Beta. recv() returned: " + std::to_string(valread));
        close(sock_alpha_to_beta);
        return -1;
    }

    std::string data(net_buf, valread);
    std::vector<std::string> params;
    size_t pos = 0;
    while ((pos = data.find("\n")) != std::string::npos) {
        params.push_back(data.substr(0, pos));
        data.erase(0, pos + 1);
    }

    // If last param is not empty, add it
    if (!data.empty()) params.push_back(data);

    // Expecting at least 7 parameters: p, q, g, r, pk_E, FHEcryptoContext, pk_F
    //if (params.size() < 7) {
    if (params.size() < 5) {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Server Alpha: Received insufficient parameters (" + std::to_string(params.size()) + ")");
        close(sock_alpha_to_beta);
        return -1;
    }

    // Extract and convert
    try {
        p = mpz_class(params[0]);
        q = mpz_class(params[1]);
        g = mpz_class(params[2]);
        r = mpz_class(params[3]);
        pk_E = mpz_class(params[4]);
        //Serial::DeserializeFromString(FHEcryptoContext, params[5]);
        //Serial::DeserializeFromString(pk_F, params[6]);
    } catch (...) {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Server Alpha: Error converting received message to expected format");
        close(sock_alpha_to_beta);
        return -1;
    }

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

    #if 0
    std::vector<int64_t> vectorOfInts1 = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12};
    Plaintext plaintext1               = FHEcryptoContext->MakePackedPlaintext(vectorOfInts1);

    // The encoded vectors are encrypted
    //auto FHE_c1 = FHEcryptoContext->Encrypt(pk_F, plaintext1);
    auto FHE_c1 = FHE_encSingleMsg(plaintext1);

    std::string msg = m1.get_str() + "\n" + m2.get_str() + "\n" + m3.get_str() + "\n" + m4.get_str() + "\n" + c11.get_str() + "\n" + c12.get_str() + "\n" + c21.get_str() + "\n" + c22.get_str() + "\n" + c31.get_str() + "\n" + c32.get_str() + "\n" + c41.get_str() + "\n" + c42.get_str() + "\n" + Serial::SerializeToString(FHE_c1) + "\n";
    #else
    std::string msg = m1.get_str() + "\n" + m2.get_str() + "\n" + m3.get_str() + "\n" + m4.get_str() + "\n" + c11.get_str() + "\n" + c12.get_str() + "\n" + c21.get_str() + "\n" + c22.get_str() + "\n" + c31.get_str() + "\n" + c32.get_str() + "\n" + c41.get_str() + "\n" + c42.get_str() + "\n";
    #endif
    ret = send(sock_alpha_to_beta, msg.c_str(), msg.size(), 0);

    if (ret != msg.size()) {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Failed to send the El-Gamal ciphertexts to Server Beta");
        return;
    }

    return;
}

int main()
{
    InitSrv_alpha();

    TestSrv_alpha();

    FinSrv_alpha();

    return 0;
}