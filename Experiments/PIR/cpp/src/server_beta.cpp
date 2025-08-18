#include <iostream>
#include <vector>
#include <random>
#include <cmath>    // std::sqrt, std::ceil
#include <algorithm> // std::swap
#include "pir_common.h"


//Global variables
static char net_buf[NET_BUF_SZ] = {0};
static int sock_beta_alpha_srv = -1, sock_beta_alpha_con = -1;
static int sock_beta_gamma_srv = -1, sock_beta_gamma_con = -1;


// Function declarations
static void Init_p_q_g_r(int p_bits = 3072, int q_bits = 256, int r_bits = 64);// Initializes p, q, g, GG(cyclic group) and r
static int shuffle();
static void TestSrv_beta();
static int shuffle();
static int FinSrv_beta();
static int InitSrv_beta();
static int OneTimeInitialization();
static int SendInitializedParamsToAllServers();

static void Init_p_q_g_r(int p_bits, int q_bits, int r_bits) {
    gmp_randclass rng(gmp_randinit_default);
    rng.seed(time(NULL));
    
    //Randomly choose q
    do {
        q = rng.get_z_bits(q_bits);
        mpz_nextprime(q.get_mpz_t(), q.get_mpz_t());
    } while (mpz_sizeinbase(q.get_mpz_t(), 2) != q_bits);

    //Accordingly choose p
    mpz_class temp;
    do {
        temp = rng.get_z_bits(p_bits - q_bits);
        p = temp * q + 1;
    } while (!mpz_probab_prime_p(p.get_mpz_t(), 25));

    //Then determine the generator g, which will define the cyclic group GG
    mpz_class h, exp;
    exp = (p - 1) / q;
    do {
        h = rng.get_z_range(p - 1) + 1;
        mpz_powm(g.get_mpz_t(), h.get_mpz_t(), exp.get_mpz_t(), p.get_mpz_t());
    } while (g == 1);
 
    //Randomly choose a prime number r
    do {
        r = rng.get_z_bits(r_bits);
        mpz_nextprime(r.get_mpz_t(), r.get_mpz_t());
    } while (mpz_sizeinbase(r.get_mpz_t(), 2) != r_bits);

    PrintLog(LOG_LEVEL_TRACE, __FILE__, __LINE__, "Initialized p,q,g and r");

    return;
}

static int SendInitializedParamsToAllServers(){
    int ret = 0;

    //Send parameters to Server Alpha
    (void)sendAll(sock_beta_alpha_con, p.get_str().c_str(), p.get_str().size());
    (void)sendAll(sock_beta_alpha_con, q.get_str().c_str(), q.get_str().size());
    (void)sendAll(sock_beta_alpha_con, g.get_str().c_str(), g.get_str().size());
    (void)sendAll(sock_beta_alpha_con, r.get_str().c_str(), r.get_str().size());
    (void)sendAll(sock_beta_alpha_con, pk_E.get_str().c_str(), pk_E.get_str().size());
    (void)sendAll(sock_beta_alpha_con, Serial::SerializeToString(FHEcryptoContext).c_str(), Serial::SerializeToString(FHEcryptoContext).size());
    (void)sendAll(sock_beta_alpha_con, Serial::SerializeToString(pk_F).c_str(), Serial::SerializeToString(pk_F).size());

    #if 0//Now the gamma server does not exist
    ret = send(sock_beta_to_gamma, msg.c_str(), msg.size(), 0);
    if (ret != msg.size()) {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Failed to send all parameters to Server Gamma");
        return -1;
    } else {
            PrintLog(LOG_LEVEL_TRACE, __FILE__, __LINE__, "Successfully sent all parameters to Server Gamma");
    }
    #endif

    return 0;
}

static int OneTimeInitialization(){
    int ret = 0;

    //Initialize p, q, g, GG(cyclic group) and r
    Init_p_q_g_r(P_BITS, Q_BITS, R_BITS);

    //Initialize El-Gamal key-pair
    std::tie(pk_E, sk_E) = ElGamal_keyGen(p, q, g);
    PrintLog(LOG_LEVEL_TRACE, __FILE__, __LINE__, "Generated El-Gamal key-pair");

    //Initialize FHE key-pair
    ret = FHE_keyGen();//TODO: Check allocation, call by reference etc.
    
    if (ret != 0) {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Failed to generate FHE key-pair");
        return -1;
    } else {
        PrintLog(LOG_LEVEL_TRACE, __FILE__, __LINE__, "Generated FHE key-pair");
    }

    //Initialize sockets for communication with server alpha
    ret = InitAcceptingSocket(BETA_LISTENING_TO_ALPHA_PORT, &sock_beta_alpha_srv, &sock_beta_alpha_con);

    if (ret != 0) {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Cannot establish communication with Server Alpha!!");
        return -1;
    }

    #if 0
    ret = InitAcceptingSocket(BETA_LISTENING_TO_GAMMA_PORT, &sock_beta_gamma_srv, &sock_beta_gamma_con);

    if (ret != 0) {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Cannot establish communication with Server Gamma!!");
        return -1;
    }
    #endif

    //TODO send {p, q, r, g, pk_E, pk_F} to other parties over network in serialized format
    ret = SendInitializedParamsToAllServers();

    //Other parties must store them in their own pir_common.cpp file

    return ret;
}

static int InitSrv_beta(){
    int ret = 0;

    // Initialize server beta
    ret = OneTimeInitialization();

    if (ret != 0) {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Failed to initialize Server Beta");
        return -1;
    } else {
        PrintLog(LOG_LEVEL_INFO, __FILE__, __LINE__, "Server Beta initialization complete");
    }

    return ret;
}

static int FinSrv_beta(){
    int ret = 0;

    // Close the sockets
    if (sock_beta_alpha_srv != -1) {
        close(sock_beta_alpha_srv);
        sock_beta_alpha_srv = -1;
    }
    if (sock_beta_alpha_con != -1) {
        close(sock_beta_alpha_con);
        sock_beta_alpha_con = -1;
    }

    PrintLog(LOG_LEVEL_SPECIAL, __FILE__, __LINE__, "Finalized Server Beta");

    return ret;
}

static int shuffle() {
    int ret = 0;
    // compute size M = N + ceil(sqrt(N))
    unsigned int sqrt_N = static_cast<unsigned int>(std::ceil(std::sqrt(static_cast<double>(N))));
    unsigned int M = N + sqrt_N;

    // build SS = {1, 2, ..., (N + sqrt_N))}
    std::vector<unsigned int> SS;
    SS.reserve(M);
    for (unsigned int i = 1; i <= M; ++i) SS.push_back(i);

    if (SS.empty()) return -1; // nothing to do for N == 0

    // RNG: mt19937 seeded from random_device
    std::random_device rd;
    std::mt19937 gen(rd());

    // repeatedly pick a random index in [0, SS.size()-1], print element,
    // then remove it by swapping with the last element and pop_back()
    for (unsigned int iter = 0; iter < M; ++iter) {
        std::uniform_int_distribution<std::size_t> dist(0, SS.size() - 1);
        std::size_t idx = dist(gen);
        unsigned int chosen = SS[idx];

        std::cout << chosen;
        if (iter + 1 < M) std::cout << ' ';

        // remove chosen element (order not preserved)
        std::swap(SS[idx], SS.back());
        SS.pop_back();
    }
    std::cout << '\n';

    return 0;
}

static void TestSrv_beta() {
    int ret = -1;

    size_t received_sz = 0;
    int ret_recv = 0;

    mpz_class m1_local, m2_local, m3_local, m4_local, c11_local, c12_local, c21_local, c22_local, c31_local, c32_local, c41_local, c42_local;
    Ciphertext<DCRTPoly> FHE_c1_local;

    // Receive m1
    ret_recv = recvAll(sock_beta_alpha_con, net_buf, sizeof(net_buf), &received_sz);
    if (ret_recv != 0) {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Server Beta: Failed to receive m1 from Server Alpha");
        close(sock_beta_alpha_con);
        return;
    }
    m1_local = mpz_class(std::string(net_buf, received_sz));

    // Receive m2
    ret_recv = recvAll(sock_beta_alpha_con, net_buf, sizeof(net_buf), &received_sz);
    if (ret_recv != 0) {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Server Beta: Failed to receive m2 from Server Alpha");
        close(sock_beta_alpha_con);
        return;
    }
    m2_local = mpz_class(std::string(net_buf, received_sz));

    // Receive m3
    ret_recv = recvAll(sock_beta_alpha_con, net_buf, sizeof(net_buf), &received_sz);
    if (ret_recv != 0) {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Server Beta: Failed to receive m3 from Server Alpha");
        close(sock_beta_alpha_con);
        return;
    }
    m3_local = mpz_class(std::string(net_buf, received_sz));

    // Receive m4
    ret_recv = recvAll(sock_beta_alpha_con, net_buf, sizeof(net_buf), &received_sz);
    if (ret_recv != 0) {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Server Beta: Failed to receive m4 from Server Alpha");
        close(sock_beta_alpha_con);
        return;
    }
    m4_local = mpz_class(std::string(net_buf, received_sz));

    // Receive c11
    ret_recv = recvAll(sock_beta_alpha_con, net_buf, sizeof(net_buf), &received_sz);
    if (ret_recv != 0) {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Server Beta: Failed to receive c11 from Server Alpha");
        close(sock_beta_alpha_con);
        return;
    }
    c11_local = mpz_class(std::string(net_buf, received_sz));

    // Receive c12
    ret_recv = recvAll(sock_beta_alpha_con, net_buf, sizeof(net_buf), &received_sz);
    if (ret_recv != 0) {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Server Beta: Failed to receive c12 from Server Alpha");
        close(sock_beta_alpha_con);
        return;
    }
    c12_local = mpz_class(std::string(net_buf, received_sz));

    // Receive c21
    ret_recv = recvAll(sock_beta_alpha_con, net_buf, sizeof(net_buf), &received_sz);
    if (ret_recv != 0) {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Server Beta: Failed to receive c21 from Server Alpha");
        close(sock_beta_alpha_con);
        return;
    }
    c21_local = mpz_class(std::string(net_buf, received_sz));

    // Receive c22
    ret_recv = recvAll(sock_beta_alpha_con, net_buf, sizeof(net_buf), &received_sz);
    if (ret_recv != 0) {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Server Beta: Failed to receive c22 from Server Alpha");
        close(sock_beta_alpha_con);
        return;
    }
    c22_local = mpz_class(std::string(net_buf, received_sz));

    // Receive c31
    ret_recv = recvAll(sock_beta_alpha_con, net_buf, sizeof(net_buf), &received_sz);
    if (ret_recv != 0) {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Server Beta: Failed to receive c31 from Server Alpha");
        close(sock_beta_alpha_con);
        return;
    }
    c31_local = mpz_class(std::string(net_buf, received_sz));

    // Receive c32
    ret_recv = recvAll(sock_beta_alpha_con, net_buf, sizeof(net_buf), &received_sz);
    if (ret_recv != 0) {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Server Beta: Failed to receive c32 from Server Alpha");
        close(sock_beta_alpha_con);
        return;
    }
    c32_local = mpz_class(std::string(net_buf, received_sz));

    // Receive c41
    ret_recv = recvAll(sock_beta_alpha_con, net_buf, sizeof(net_buf), &received_sz);
    if (ret_recv != 0) {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Server Beta: Failed to receive c41 from Server Alpha");
        close(sock_beta_alpha_con);
        return;
    }
    c41_local = mpz_class(std::string(net_buf, received_sz));

    // Receive c42
    ret_recv = recvAll(sock_beta_alpha_con, net_buf, sizeof(net_buf), &received_sz);
    if (ret_recv != 0) {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Server Beta: Failed to receive c42 from Server Alpha");
        close(sock_beta_alpha_con);
        return;
    }
    c42_local = mpz_class(std::string(net_buf, received_sz));

    // Receive FHE_c1
    ret_recv = recvAll(sock_beta_alpha_con, net_buf, sizeof(net_buf), &received_sz);
    if (ret_recv != 0)
    {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Failed to receive FHE_c1 from Server Beta");
        close(sock_beta_alpha_con);
        return;
    }
    Serial::DeserializeFromString(FHE_c1_local, std::string(net_buf, received_sz));

    // Decrypt and check
    mpz_class decrypted_m1 = ElGamal_decrypt({c11_local, c12_local}, sk_E);
    mpz_class decrypted_m2 = ElGamal_decrypt({c21_local, c22_local}, sk_E);
    mpz_class decrypted_m3 = ElGamal_decrypt({c31_local, c32_local}, sk_E);
    mpz_class decrypted_m4 = ElGamal_decrypt({c41_local, c42_local}, sk_E);

    if (decrypted_m1 != m1_local) {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Server Beta: Decrypted m1: " + decrypted_m1.get_str() + " does not match with expected value: " + m1_local.get_str() + " !!");
    } else {
        PrintLog(LOG_LEVEL_INFO, __FILE__, __LINE__, "Server Beta: Decrypted m1 matches with expected value");
    }

    if (decrypted_m2 != m2_local) {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Server Beta: Decrypted m2: " + decrypted_m2.get_str() + " does not match with expected value: " + m2_local.get_str() + " !!");
    } else {
        PrintLog(LOG_LEVEL_INFO, __FILE__, __LINE__, "Server Beta: Decrypted m2 matches with expected value");
    }

    if (decrypted_m3 != m3_local) {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Server Beta: Decrypted m3: " + decrypted_m3.get_str() + " does not match with expected value: " + m3_local.get_str() + " !!");
    } else {
        PrintLog(LOG_LEVEL_INFO, __FILE__, __LINE__, "Server Beta: Decrypted m3 matches with expected value");
    }

    if (decrypted_m4 != m4_local) {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Server Beta: Decrypted m4: " + decrypted_m4.get_str() + " does not match with expected value: " + m4_local.get_str() + " !!");
    } else {
        PrintLog(LOG_LEVEL_INFO, __FILE__, __LINE__, "Server Beta: Decrypted m4 matches with expected value");
    }

    //Decrypt FHE_c1_local
    Plaintext decrypted_PIRBlock;
    FHEcryptoContext->Decrypt(sk_F, FHE_c1_local, &decrypted_PIRBlock);

    //std::vector<int64_t> decodedAddResult = plaintextAddResult->GetPackedValue<int64_t>();

    decrypted_PIRBlock->SetLength(NUM_FHE_BLOCKS_PER_PIR_BLOCK);//TODO: What to set?

    for (int i = 0; i < NUM_FHE_BLOCKS_PER_PIR_BLOCK; ++i) {
        PrintLog(LOG_LEVEL_INFO, __FILE__, __LINE__, "Server Beta: Decrypted decrypted_PIRBlock[" + std::to_string(i) + "]: " + std::to_string(decrypted_PIRBlock->GetPackedValue()[i]));
    }

    return;
}

int main() {

    InitSrv_beta();

    TestSrv_beta();

    FinSrv_beta();

    return 0;
}
