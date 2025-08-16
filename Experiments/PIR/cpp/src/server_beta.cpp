#include <iostream>
#include <vector>
#include <random>
#include <cmath>    // std::sqrt, std::ceil
#include <algorithm> // std::swap
#include "pir_common.h"

static int sock_beta_alpha_srv = -1, sock_beta_alpha_con = -1;
static int sock_beta_gamma_srv = -1, sock_beta_gamma_con = -1;


// Function, which initializes p, q, g, GG(cyclic group) and r
void Init_p_q_GG_r(int p_bits = 3072, int q_bits = 256, int r_bits = 64);


int shuffle();

void Init_p_q_GG_r(int p_bits, int q_bits, int r_bits) {
    gmp_randclass rng(gmp_randinit_default);
    rng.seed(time(NULL));
    
    //Randomly choose q
    do {
        q = rng.get_z_bits(q_bits);
        mpz_nextprime(q.get_mpz_t(), q.get_mpz_t());
    } while (mpz_sizeinbase(q.get_mpz_t(), 2) != q_bits);

    //Accordingly choose q
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
}

int SendInitializedParamsToAllServers(){
    int ret = 0;

    //Try to send p, q, g, r, pk_E, TODO: pk_F
    // Optionally print for debug
    std::cout << "Server Beta: Sending parameters:" << std::endl;
    std::cout << "p: " << p.get_str() << std::endl;
    std::cout << "q: " << q.get_str() << std::endl;
    std::cout << "g: " << g.get_str() << std::endl;
    std::cout << "r: " << r.get_str() << std::endl;
    std::cout << "pk_E: " << pk_E.get_str() << std::endl;

    std::string msg = p.get_str() + "\n" + q.get_str() + "\n" + g.get_str() + "\n" + r.get_str() + "\n" + pk_E.get_str() + "\n" +  "\n";
    ret = send(sock_beta_alpha_con, msg.c_str(), msg.size(), 0);

    if (ret != msg.size()) {
        std::cerr << "Server Beta: Failed to send all parameters to Server Alpha" << std::endl;
        return -1;
    }

    #if 0//Now the gamma server does not exist
    ret = send(sock_beta_to_gamma, msg.c_str(), msg.size(), 0);
    if (ret != msg.size()) {
        std::cerr << "Server Beta: Failed to send all parameters to Server Gamma" << std::endl;
        return -1;
    }
    #endif

    return ret;
}

int OneTimeInitialization(){
    int ret = 0;

    ret = InitAcceptingSocket(BETA_LISTENING_TO_ALPHA_PORT, &sock_beta_alpha_srv, &sock_beta_alpha_con);

    if (ret != 0) {
        std::cerr << "Server Beta: Failed to initialize accepting socket for Server Alpha" << std::endl;
        return -1;
    }

    #if 0
    ret = InitAcceptingSocket(BETA_LISTENING_TO_GAMMA_PORT, &sock_beta_gamma_srv, &sock_beta_gamma_con);

    if (ret != 0) {
        std::cerr << "Server Beta: Failed to initialize accepting socket for Server Alpha" << std::endl;
        return -1;
    }
    #endif

    //Initialize p, q, g, GG(cyclic group) and r
    Init_p_q_GG_r(P_BITS, Q_BITS, R_BITS);

    //Initialize El-Gamal key-pair
    std::tie(pk_E, sk_E) = ElGamal_keyGen(p, q, g);

    //TODO Initialize FHE key-pair

    //TODO send {p, q, r, g, pk_E, pk_F} to other parties over network in serialized format
    ret = SendInitializedParamsToAllServers();


    //Other parties must store them in their own pir_common.cpp file

    // Optionally print initialized values for debug
    std::cout << "Initialized ElGamal params:" << std::endl;
    std::cout << "p: " << p.get_str() << std::endl;
    std::cout << "q: " << q.get_str() << std::endl;
    std::cout << "g: " << g.get_str() << std::endl;
    return ret;
}

int finalize_server_beta(){
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

    return ret;
}

int shuffle() {
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

int main() {
    OneTimeInitialization(); // example
    return 0;
}
