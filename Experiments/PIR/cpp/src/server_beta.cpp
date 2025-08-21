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
static void Init_parameters(int p_bits = 3072, int q_bits = 256, int r_bits = 64);// Initializes p, q, g, GG(cyclic group) and r
static int shuffle();
static void TestSrv_beta();
static int shuffle();
static int FinSrv_beta();
static int InitSrv_beta();
static int OneTimeInitialization();
static int SendInitializedParamsToAllServers();
static mpz_class selectRho();
static void TestBlindedExponentiation();
static void TestBlindedExponentiation1();
static void TestBlindedExponentiation2();

static void Init_parameters(int p_bits, int q_bits, int r_bits) {
 
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

    //Choose 2 as the generator of the multiplicative sub-group ZZ_q*
    g_q = mpz_class(2);
 
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
    (void)sendAll(sock_beta_alpha_con, g_q.get_str().c_str(), g_q.get_str().size());
    (void)sendAll(sock_beta_alpha_con, r.get_str().c_str(), r.get_str().size());
    (void)sendAll(sock_beta_alpha_con, pk_E.get_str().c_str(), pk_E.get_str().size());
    (void)sendAll(sock_beta_alpha_con, pk_E_q.get_str().c_str(), pk_E_q.get_str().size());
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

    // Initialize random number generator
    // good-ish seed source: std::random_device (combine two samples)
    std::random_device rd;
    unsigned long seed = (static_cast<unsigned long>(rd()) << 1) ^ rd();
    rng.seed(seed); // seed() seeds the gmp_randclass

    //Initialize p, q, g, GG(cyclic group) and r
    Init_parameters(P_BITS, Q_BITS, R_BITS);

    //Initialize El-Gamal key-pair
    std::tie(pk_E, sk_E) = ElGamal_keyGen();
    PrintLog(LOG_LEVEL_TRACE, __FILE__, __LINE__, "Generated El-Gamal key-pair");

    //Initialize El-Gamal_q key-pair
    std::tie(pk_E_q, sk_E_q) = ElGamal_q_keyGen();
    PrintLog(LOG_LEVEL_TRACE, __FILE__, __LINE__, "Generated El-Gamal key-pair in ZZ*_q");

    //TODO: Temporary, just for testing. Calling from here, so that server_alpha is not required to be executed now
    //TestBlindedExponentiation();
    //TestBlindedExponentiation1();
    TestBlindedExponentiation2();
    return 0;

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

static mpz_class selectRho() {

    mpz_class rho = ElGamal_randomGroupElement();

    // rho must not be divided by q
    while ((rho % q) == 0) {
        rho = ElGamal_randomGroupElement();
    }

    //PrintLog(LOG_LEVEL_INFO, __FILE__, __LINE__, "Selected rho: " + rho.get_str());

    return rho;
}

static void TestBlindedExponentiation() {
    int ret;

    PrintLog(LOG_LEVEL_INFO, __FILE__, __LINE__, "El-Gamal parameters p: " + p.get_str() + " q: " + q.get_str()+ " g: " + g.get_str()+ " pk_E: " + pk_E.get_str()+ " sk_E: " + sk_E.get_str());

    mpz_class Rho = selectRho();

    mpz_class h = rng.get_z_range(q-1) + 1;//To ensure that the element belongs to ZZ_q*
    mpz_class h_1;
    mpz_invert(h_1.get_mpz_t(), h.get_mpz_t(), q.get_mpz_t());//Compute h^-1 in ZZ_q*

    mpz_class alpha = rng.get_z_range(q-1) + 1;
    mpz_class alpha_1;
    mpz_invert(alpha_1.get_mpz_t(), alpha.get_mpz_t(), q.get_mpz_t());

    //Confirm that the inverses really work
    assert((h * h_1) % q == 1);
    assert((alpha * alpha_1) % q == 1);

    mpz_class I = rng.get_z_range(N-1) + 1;//I will be in the range of [1, N]
    PrintLog(LOG_LEVEL_INFO, __FILE__, __LINE__, "Chosen plain text messages are Rho: " + Rho.get_str() + " I: " + I.get_str() + " h: " + h.get_str() + " h_1: " + h_1.get_str()+ " alpha: " + alpha.get_str() + " alpha_1: " + alpha_1.get_str() );
    mpz_class RhoExpI;
    //mpz_powm(RhoExpI.get_mpz_t(), Rho.get_mpz_t(), I.get_mpz_t(), p.get_mpz_t()); previous
    //mpz_class RhoExpI_h = (RhoExpI * h) % p; previous
    mpz_powm(RhoExpI.get_mpz_t(), Rho.get_mpz_t(), I.get_mpz_t(), q.get_mpz_t());//Since this is an exponenet
    mpz_class RhoExpI_h = (RhoExpI * h) % q;

    mpz_class gExp_RhoExpI;
    mpz_powm(gExp_RhoExpI.get_mpz_t(), g.get_mpz_t(), RhoExpI.get_mpz_t(), p.get_mpz_t());
    mpz_class gExp_RhoExpI_h;
    mpz_powm(gExp_RhoExpI_h.get_mpz_t(), g.get_mpz_t(), RhoExpI_h.get_mpz_t(), p.get_mpz_t());

    //Reduce exponents to mod q
    mpz_class RhoExpI_mod_q = RhoExpI % q;
    mpz_class gExp_RhoExpI_mod_q;
    mpz_powm(gExp_RhoExpI_mod_q.get_mpz_t(), g.get_mpz_t(), RhoExpI_mod_q.get_mpz_t(), p.get_mpz_t());
    mpz_class RhoExpI_h__mod_q = RhoExpI_h % q;
    mpz_class gExp_RhoExpI_h__mod_q;
    mpz_powm(gExp_RhoExpI_h__mod_q.get_mpz_t(), g.get_mpz_t(), RhoExpI_h__mod_q.get_mpz_t(), p.get_mpz_t());

    std::pair<mpz_class, mpz_class> E_Rho = ElGamal_encrypt(Rho, pk_E);//Compute E(\rho)
    std::pair<mpz_class, mpz_class> E_h = ElGamal_encrypt(h, pk_E);//Compute E(h)
    std::pair<mpz_class, mpz_class> E_RhoExpI = ElGamal_exp_ct(E_Rho, I, pk_E);//Compute E(\rho^I)
    std::pair<mpz_class, mpz_class> E_RhoExpImulh = ElGamal_mult_ct(E_RhoExpI, E_h);//Compute E(\rho^I.h)

    mpz_class decrypted_RhoExpI = ElGamal_decrypt(E_RhoExpI, sk_E);
    decrypted_RhoExpI = decrypted_RhoExpI % q;//Since this should be in exponent
    mpz_class decrypted_RhoExpImulh = ElGamal_decrypt(E_RhoExpImulh, sk_E);
    decrypted_RhoExpImulh = decrypted_RhoExpImulh % q;//Since this should be in exponent

    std::pair<mpz_class, mpz_class> E_gExp_RhoExpI_h__mod_q = ElGamal_encrypt(gExp_RhoExpI_h__mod_q, pk_E);//Compute E(g^{\rho^{I}.h mod q})
    std::pair<mpz_class, mpz_class> E_gExp_RhoExpI_hh_1_mod_q = ElGamal_exp_ct(E_gExp_RhoExpI_h__mod_q, h_1, pk_E);//Compute E(g^{\rho^{I}.h.h_1 mod q})
    std::pair<mpz_class, mpz_class> E_gExp_RhoExpI_mulh_alpha_mod_q = ElGamal_exp_ct(E_gExp_RhoExpI_h__mod_q, alpha, pk_E);//Compute E(g^{\rho^{I}.h.alpha})
    std::pair<mpz_class, mpz_class> E_gExp_RhoExpI_mulh_alpha_alpha_1_mod_q = ElGamal_exp_ct(E_gExp_RhoExpI_mulh_alpha_mod_q, alpha_1, pk_E);//Compute E(g^{\rho^{I}.h.alpha.h^{-1}})
    std::pair<mpz_class, mpz_class> E_gExp_RhoExpI_mulh_mulalpha_alpha_1_h_1_mod_q = ElGamal_exp_ct(E_gExp_RhoExpI_mulh_alpha_alpha_1_mod_q, h_1, pk_E);//Compute E(g^{\rho^{I}.h.alpha.h^{-1}})

    mpz_class decrypted_h = ElGamal_decrypt(E_h, sk_E);
    mpz_class decrypted_Rho = ElGamal_decrypt(E_Rho, sk_E);
    mpz_class decrypted_gExp_RhoExpI_hh_1_mod_q = ElGamal_decrypt(E_gExp_RhoExpI_hh_1_mod_q, sk_E);
    mpz_class decrypted_gExp_RhoExpI_mulh_alpha_alpha_1_mod_q = ElGamal_decrypt(E_gExp_RhoExpI_mulh_alpha_alpha_1_mod_q, sk_E);
    mpz_class decrypted_gExp_RhoExpI_mulh_mulalpha_alpha_1_mod_q = ElGamal_decrypt(E_gExp_RhoExpI_mulh_alpha_alpha_1_mod_q, sk_E);
    mpz_class decrypted_gExp_RhoExpI_mulh_mulalpha_alpha_1_h_1_mod_q = ElGamal_decrypt(E_gExp_RhoExpI_mulh_mulalpha_alpha_1_h_1_mod_q, sk_E);

    if (Rho != decrypted_Rho) {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Decryption failed for Rho: expected " + Rho.get_str() + ", got " + decrypted_Rho.get_str());
        //return;
    } else {
        PrintLog(LOG_LEVEL_INFO, __FILE__, __LINE__, "Rho matched with decrypted_Rho ");
    }

    if (h != decrypted_h) {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Decryption failed for h: expected " + h.get_str() + ", got " + decrypted_h.get_str());
        //return;
    } else {
        PrintLog(LOG_LEVEL_INFO, __FILE__, __LINE__, "h matched with decrypted_h ");
    }

    if (RhoExpI != decrypted_RhoExpI) {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Decryption failed for RhoExpI: expected " + RhoExpI.get_str() + ", got " + decrypted_RhoExpI.get_str());
        //return;
    } else {
        PrintLog(LOG_LEVEL_INFO, __FILE__, __LINE__, "RhoExpI matched with decrypted_RhoExpI ");
    }

    if (RhoExpI_h__mod_q != decrypted_RhoExpImulh) {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Decryption failed for RhoExpI_h: expected " + RhoExpI_h__mod_q.get_str() + ", got " + decrypted_RhoExpImulh.get_str());
        //return;
    } else {
        PrintLog(LOG_LEVEL_INFO, __FILE__, __LINE__, "RhoExpI_h matched with decrypted_RhoExpImulh ");
    }

    if (gExp_RhoExpI != gExp_RhoExpI_mod_q) {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "g^RhoExpI does not match with g^RhoExpI mod q: expected " + gExp_RhoExpI.get_str() + ", got " + gExp_RhoExpI_mod_q.get_str());
        //return;
    } else {
        PrintLog(LOG_LEVEL_INFO, __FILE__, __LINE__, "g^RhoExpI matches with g^RhoExpI mod q"); 
    }

    if (gExp_RhoExpI_h != gExp_RhoExpI_h__mod_q) {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "g^RhoExpI_h does not match with g^RhoExpI_h mod q: expected " + gExp_RhoExpI_h.get_str() + ", got " + gExp_RhoExpI_h__mod_q.get_str());
        //return;
    } else {
        PrintLog(LOG_LEVEL_INFO, __FILE__, __LINE__, "g^RhoExpI_h matches with g^RhoExpI_h mod q");
    }

    if (decrypted_gExp_RhoExpI_hh_1_mod_q != gExp_RhoExpI_mod_q) {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "g^{Rho^{I}} mod q does not match with g^{Rho^{I}.h.h^{-1}} mod q: expected " + gExp_RhoExpI.get_str() + ", got " + decrypted_gExp_RhoExpI_hh_1_mod_q.get_str());
        //return;
    } else {
        PrintLog(LOG_LEVEL_INFO, __FILE__, __LINE__, "g^{Rho^{I}} mod q matches with g^{Rho^{I}.h.h^{-1}} mod q..!!");
    }

    if (decrypted_gExp_RhoExpI_mulh_alpha_alpha_1_mod_q != gExp_RhoExpI_h__mod_q) {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "g^{Rho^{I}.h} mod q does not match with g^{Rho^{I}.h.alpha.alpha^{-1}} mod q: expected " + gExp_RhoExpI_h.get_str() + ", got " + decrypted_gExp_RhoExpI_mulh_alpha_alpha_1_mod_q.get_str());
        //return;
    } else {
        PrintLog(LOG_LEVEL_INFO, __FILE__, __LINE__, "g^{Rho^{I}.h} mod q matches with g^{Rho^{I}.h.alpha.alpha^{-1}} mod q..!!");
    }
#if 0
    mpz_class reduced_RhoExpImulh = decrypted_RhoExpImulh % q;

    mpz_class Exp_rhoExpImulh, Exp_RedrhoExpImulh, h_inverse, final_result, reduced_final_result;
    mpz_powm(Exp_rhoExpImulh.get_mpz_t(), g.get_mpz_t(), decrypted_rhoExpImulh.get_mpz_t(), p.get_mpz_t());
    mpz_powm(Exp_RedrhoExpImulh.get_mpz_t(), g.get_mpz_t(), reduced_rhoExpImulh.get_mpz_t(), p.get_mpz_t());
    mpz_invert(h_inverse.get_mpz_t(), h.get_mpz_t(), p.get_mpz_t());
    mpz_class reduced_h_inverse = h_inverse % q;
    mpz_powm(final_result.get_mpz_t(), Exp_rhoExpImulh.get_mpz_t(), h_inverse.get_mpz_t(), p.get_mpz_t());
    mpz_powm(reduced_final_result.get_mpz_t(), Exp_RedrhoExpImulh.get_mpz_t(), reduced_h_inverse.get_mpz_t(), p.get_mpz_t());

    PrintLog(LOG_LEVEL_INFO, __FILE__, __LINE__, "Decrypted messages are decrypted_rhoExpI: " + decrypted_rhoExpI.get_str() + ", decrypted_rhoExpImulh: " + decrypted_rhoExpImulh.get_str() + ", Exp_rhoExpImulh: " + Exp_rhoExpImulh.get_str() + ", Exp_RedrhoExpImulh: " + Exp_RedrhoExpImulh.get_str());
    PrintLog(LOG_LEVEL_INFO, __FILE__, __LINE__, "h_inverse: " + h_inverse.get_str() + ", reduced_h_inverse: " + reduced_h_inverse.get_str() + ", final_result: " + final_result.get_str() + ", reduced_final_result: " + reduced_final_result.get_str());
#endif
}

static void TestBlindedExponentiation1() {
    int ret;

    PrintLog(LOG_LEVEL_INFO, __FILE__, __LINE__, "El-Gamal parameters p: " + p.get_str() + " q: " + q.get_str()+ " g: " + g.get_str()+ " pk_E: " + pk_E.get_str()+ " sk_E: " + sk_E.get_str());

    mpz_class Rho = selectRho();

    mpz_class h = rng.get_z_range(q-1) + 1;//To ensure that the element belongs to ZZ_q*
    mpz_class h_1;
    mpz_invert(h_1.get_mpz_t(), h.get_mpz_t(), q.get_mpz_t());//Compute h^-1 in ZZ_q*

    mpz_class alpha = rng.get_z_range(q-1) + 1;
    mpz_class alpha_1;
    mpz_invert(alpha_1.get_mpz_t(), alpha.get_mpz_t(), q.get_mpz_t());

    //Confirm that the inverses really work
    assert((h * h_1) % q == 1);
    assert((alpha * alpha_1) % q == 1);

    mpz_class I = rng.get_z_range(
        N-1) + 1;//I will be in the range of [1, N]
    PrintLog(LOG_LEVEL_INFO, __FILE__, __LINE__, "Chosen plain text messages are Rho: " + Rho.get_str() + " I: " + I.get_str() + " h: " + h.get_str() + " h_1: " + h_1.get_str()+ " alpha: " + alpha.get_str() + " alpha_1: " + alpha_1.get_str() );
    mpz_class RhoExpI;
    //mpz_powm(RhoExpI.get_mpz_t(), Rho.get_mpz_t(), I.get_mpz_t(), p.get_mpz_t()); previous
    //mpz_class RhoExpI_h = (RhoExpI * h) % p; previous
    mpz_powm(RhoExpI.get_mpz_t(), Rho.get_mpz_t(), I.get_mpz_t(), q.get_mpz_t());//Since this is an exponenet
    //mpz_class RhoExpI_h = (RhoExpI * h) % q;

    mpz_class gExp_RhoExpI;
    mpz_powm(gExp_RhoExpI.get_mpz_t(), g.get_mpz_t(), RhoExpI.get_mpz_t(), p.get_mpz_t());

    //Server beta broadcasts E_Rho
    std::pair<mpz_class, mpz_class> E_Rho = ElGamal_encrypt(Rho, pk_E);//Compute E(\rho)

    //Client computes E_RhoExpI
    std::pair<mpz_class, mpz_class> E_RhoExpI = ElGamal_exp_ct(E_Rho, I, pk_E);//Compute E(\rho^I)
    //Client computes E_RhoExpI_h and sends to server beta
    std::pair<mpz_class, mpz_class> E_h = ElGamal_encrypt(h, pk_E);//Compute E(h)
    std::pair<mpz_class, mpz_class> E_RhoExpI_h = ElGamal_mult_ct(E_RhoExpI, E_h);//Compute E(\rho^I.h)

    //Server beta decrypts and reduces to mod q
    mpz_class decrypted_RhoExpI_h = ElGamal_decrypt(E_RhoExpI_h, sk_E);
    mpz_class decrypted_RhoExpI_h__modq = decrypted_RhoExpI_h % q;//Reduce to mod q
    //Server beta raises to g
    mpz_class gExp_RhoExpI_h__modq;
    mpz_powm(gExp_RhoExpI_h__modq.get_mpz_t(), g.get_mpz_t(), decrypted_RhoExpI_h__modq.get_mpz_t(), p.get_mpz_t());
    std::pair<mpz_class, mpz_class> E_gExp_RhoExpI_h__modq = ElGamal_encrypt(gExp_RhoExpI_h__modq, pk_E);//Sends corresponding ciphertext

    //Server alpha removes h by raising it to h^{-1}
    std::pair<mpz_class, mpz_class> E_gExp_RhoExpI_hh_1__modq = ElGamal_exp_ct(E_gExp_RhoExpI_h__modq, h_1, pk_E);

    //If decrypted must match with gExp_RhoExpI
    mpz_class decrypted_gExp_RhoExpI_hh_1__modq = ElGamal_decrypt(E_gExp_RhoExpI_hh_1__modq, sk_E);

    assert(decrypted_gExp_RhoExpI_hh_1__modq == gExp_RhoExpI);

    if (decrypted_gExp_RhoExpI_hh_1__modq == gExp_RhoExpI) {
        PrintLog(LOG_LEVEL_INFO, __FILE__, __LINE__, "Server Beta: Decryption successful and matches with gExp_RhoExpI");
    } else {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Server Beta: Decryption failed or does not match with gExp_RhoExpI. Expected :" + gExp_RhoExpI.get_str() + " but got: " + decrypted_gExp_RhoExpI_hh_1__modq.get_str());
    }
}

static void TestBlindedExponentiation2() {
    int ret;

    PrintLog(LOG_LEVEL_INFO, __FILE__, __LINE__, "El-Gamal parameters p: " + p.get_str() + " q: " + q.get_str()+ " g: " + g.get_str()+ " pk_E: " + pk_E.get_str()+ " sk_E: " + sk_E.get_str());
    //Choose random message and random exponent
    mpz_class m1 = ElGamal_randomGroupElement();
    mpz_class m2 = ElGamal_randomGroupElement();
    mpz_class exp = rng.get_z_range(q);
    PrintLog(LOG_LEVEL_INFO, __FILE__, __LINE__, "Chosen plaintext messages are m1: " + m1.get_str() + " m2: " + m2.get_str() + " exp: " + exp.get_str());

    mpz_class m3 = (m1*m2)%p;
    mpz_class m4;
    mpz_powm(m4.get_mpz_t(), m1.get_mpz_t(), exp.get_mpz_t(), p.get_mpz_t());

    mpz_class m5 = (m4*m2)%p;//(m1^exp)*m2

    std::pair<mpz_class, mpz_class> E_m1 = ElGamal_encrypt(m1, pk_E);
    std::pair<mpz_class, mpz_class> E_m2 = ElGamal_encrypt(m2, pk_E);
    std::pair<mpz_class, mpz_class> E_m3 = ElGamal_mult_ct(E_m1, E_m2);
    std::pair<mpz_class, mpz_class> E_m4 = ElGamal_exp_ct(E_m1, exp, pk_E);
    std::pair<mpz_class, mpz_class> E_m5 = ElGamal_mult_ct(E_m4, E_m2);

    mpz_class decrypted_m1 = ElGamal_decrypt(E_m1, sk_E);
    mpz_class decrypted_m3 = ElGamal_decrypt(E_m3, sk_E);
    mpz_class decrypted_m4 = ElGamal_decrypt(E_m4, sk_E);
    mpz_class decrypted_m5 = ElGamal_decrypt(E_m5, sk_E);

    if (decrypted_m1 == m1) {
        PrintLog(LOG_LEVEL_INFO, __FILE__, __LINE__, "El-Gamal encryption works");
    } else {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "El-Gamal encryption is not working. Expected: " + m1.get_str() + " but got: " + decrypted_m1.get_str());
    }

    if (decrypted_m3 == m3) {
        PrintLog(LOG_LEVEL_INFO, __FILE__, __LINE__, "El-Gamal multiplication works");
    } else {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "El-Gamal multiplication is not working. Expected: " + m3.get_str() + " but got: " + decrypted_m3.get_str());
    }

    if (decrypted_m4 == m4) {
        PrintLog(LOG_LEVEL_INFO, __FILE__, __LINE__, "El-Gamal exponentiation works");
    } else {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "El-Gamal exponentiation is not working. Expected: " + m4.get_str() + " but got: " + decrypted_m4.get_str());
    }

    if (decrypted_m5 == m5) {
        PrintLog(LOG_LEVEL_INFO, __FILE__, __LINE__, "El-Gamal multiplication after exponentiation works");
    } else {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "El-Gamal multiplication after exponentiation is not working. Expected: " + m5.get_str() + " but got: " + decrypted_m5.get_str());
    }

    PrintLog(LOG_LEVEL_INFO, __FILE__, __LINE__, "El-Gamal in ZZ*_q parameters g_q: " + g_q.get_str()+ " pk_E_q: " + pk_E_q.get_str()+ " sk_E_q: " + sk_E_q.get_str());

    mpz_class Rho = rng.get_z_range(q-1)+1;//i.e., within ZZ_q*
    mpz_class h = rng.get_z_range(q-1)+1;//i.e., within ZZ_q*
    mpz_class h_1;
    mpz_invert(h_1.get_mpz_t(), h.get_mpz_t(), q.get_mpz_t());

    mpz_class alpha = rng.get_z_range(q-1)+1;//i.e., within ZZ_q*
    mpz_class alpha_1;
    mpz_invert(alpha_1.get_mpz_t(), alpha.get_mpz_t(), q.get_mpz_t());

    mpz_class I = rng.get_z_range(q-1)+1;//i.e., within ZZ_q*
    PrintLog(LOG_LEVEL_INFO, __FILE__, __LINE__, "Chosen plaintext messages are Rho: " + Rho.get_str() + " h: " + h.get_str()+ " h_1: " + h_1.get_str() + " I: " + I.get_str() + " alpha: " + alpha.get_str() + " alpha_1: " + alpha_1.get_str());

    mpz_class Rho_h = (Rho*h)%q;
    mpz_class Rho_pow_I;
    mpz_powm(Rho_pow_I.get_mpz_t(), Rho.get_mpz_t(), I.get_mpz_t(), q.get_mpz_t());
    mpz_class g_pow_Rho_pow_I;
    mpz_powm(g_pow_Rho_pow_I.get_mpz_t(), g.get_mpz_t(), Rho_pow_I.get_mpz_t(), p.get_mpz_t());

    mpz_class Rho_pow_I__h = (Rho_pow_I*h)%q;//(Rho^I)*h

    std::pair<mpz_class, mpz_class> E_Rho = ElGamal_q_encrypt(Rho, pk_E_q);
    std::pair<mpz_class, mpz_class> E_h = ElGamal_q_encrypt(h, pk_E_q);
    std::pair<mpz_class, mpz_class> E_Rho_h = ElGamal_q_mult_ct(E_Rho, E_h);
    std::pair<mpz_class, mpz_class> E_Rho_pow_I = ElGamal_q_exp_ct(E_Rho, I, pk_E_q);
    std::pair<mpz_class, mpz_class> E_Rho_pow_I__h = ElGamal_q_mult_ct(E_Rho_pow_I, E_h);

    mpz_class decrypted_Rho = ElGamal_q_decrypt(E_Rho, sk_E_q);
    mpz_class decrypted_Rho_h = ElGamal_q_decrypt(E_Rho_h, sk_E_q);
    mpz_class decrypted_Rho_pow_I = ElGamal_q_decrypt(E_Rho_pow_I, sk_E_q);
    mpz_class decrypted_Rho_pow_I__h = ElGamal_q_decrypt(E_Rho_pow_I__h, sk_E_q);
   
    //Server beta decrypts and perform g^{decrypted_Rho_pow_I__h} mod p
    mpz_class g_pow_Rho_pow_I__h;
    mpz_powm(g_pow_Rho_pow_I__h.get_mpz_t(), g.get_mpz_t(), decrypted_Rho_pow_I__h.get_mpz_t(), p.get_mpz_t());
    //Encrypts and sends that under ElGamal encryption in GG
    std::pair<mpz_class, mpz_class> E_g_pow_Rho_pow_I__h = ElGamal_encrypt(g_pow_Rho_pow_I__h, pk_E);

    //Server alpha semi-homomorphically raises that to alpha under ElGamal encryption in GG
    std::pair<mpz_class, mpz_class> E_g_pow_Rho_pow_I__h_alpha = ElGamal_exp_ct(E_g_pow_Rho_pow_I__h, alpha, pk_E);

    //Client semi-homomorphically raises that to h_1 under ElGamal encryption in GG
    std::pair<mpz_class, mpz_class> E_g_pow_Rho_pow_I__h_alpha_h_1 = ElGamal_exp_ct(E_g_pow_Rho_pow_I__h_alpha, h_1, pk_E);

    //Server semi-homomorphically raises that to alpha_1 under ElGamal encryption in GG
    std::pair<mpz_class, mpz_class> E_g_pow_Rho_pow_I__h_alpha_h_1_alpha_1 = ElGamal_exp_ct(E_g_pow_Rho_pow_I__h_alpha_h_1, alpha_1, pk_E);

    //Server beta decrypts that under ElGamal encryption in GG
    mpz_class decrypted_g_pow_Rho_pow_I__h_alpha_h_1_alpha_1 = ElGamal_decrypt(E_g_pow_Rho_pow_I__h_alpha_h_1_alpha_1, sk_E);

    if (decrypted_Rho == Rho) {
        PrintLog(LOG_LEVEL_INFO, __FILE__, __LINE__, "El-Gamal encryption works in ZZ*_q");
    } else {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "El-Gamal encryption is not working in ZZ*_q. Expected: " + Rho.get_str() + " but got: " + decrypted_Rho.get_str());
    }

    if (decrypted_Rho_h == Rho_h) {
        PrintLog(LOG_LEVEL_INFO, __FILE__, __LINE__, "El-Gamal multiplication works in ZZ*_q");
    } else {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "El-Gamal multiplication is not working in ZZ*_q. Expected: " + Rho_h.get_str() + " but got: " + decrypted_Rho_h.get_str());
    }

    if (decrypted_Rho_pow_I == Rho_pow_I) {
        PrintLog(LOG_LEVEL_INFO, __FILE__, __LINE__, "El-Gamal exponentiation works in ZZ*_q");
    } else {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "El-Gamal exponentiation is not working in ZZ*_q. Expected: " + Rho_pow_I.get_str() + " but got: " + decrypted_Rho_pow_I.get_str());
    }

    if (decrypted_Rho_pow_I__h == Rho_pow_I__h) {
        PrintLog(LOG_LEVEL_INFO, __FILE__, __LINE__, "El-Gamal multiplication after exponentiation works in ZZ*_q");
    } else {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "El-Gamal multiplication after exponentiation is not working in ZZ*_q. Expected: " + Rho_pow_I__h.get_str() + " but got: " + decrypted_Rho_pow_I__h.get_str());
    }

    if (decrypted_g_pow_Rho_pow_I__h_alpha_h_1_alpha_1 == g_pow_Rho_pow_I) {
        PrintLog(LOG_LEVEL_INFO, __FILE__, __LINE__, "Finding tag homomorphically works..!!..Ha ha..Thank you..:) :) :) :) :)");
    } else {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Finding tag homomorphically is not working :( Expected: " + g_pow_Rho_pow_I.get_str() + " but got: " + decrypted_g_pow_Rho_pow_I__h_alpha_h_1_alpha_1.get_str());
    }


    #if 0
    mpz_class Rho = selectRho();

    mpz_class h = rng.get_z_range(q-1) + 1;//To ensure that the element belongs to ZZ_q*
    mpz_class h_1;
    mpz_invert(h_1.get_mpz_t(), h.get_mpz_t(), q.get_mpz_t());//Compute h^-1 in ZZ_q*

    mpz_class alpha = rng.get_z_range(q-1) + 1;
    mpz_class alpha_1;
    mpz_invert(alpha_1.get_mpz_t(), alpha.get_mpz_t(), q.get_mpz_t());

    //Confirm that the inverses really work
    assert((h * h_1) % q == 1);
    assert((alpha * alpha_1) % q == 1);

    mpz_class I = rng.get_z_range(
        N-1) + 1;//I will be in the range of [1, N]
    PrintLog(LOG_LEVEL_INFO, __FILE__, __LINE__, "Chosen plain text messages are Rho: " + Rho.get_str() + " I: " + I.get_str() + " h: " + h.get_str() + " h_1: " + h_1.get_str()+ " alpha: " + alpha.get_str() + " alpha_1: " + alpha_1.get_str() );
    mpz_class RhoExpI;
    //mpz_powm(RhoExpI.get_mpz_t(), Rho.get_mpz_t(), I.get_mpz_t(), p.get_mpz_t()); previous
    //mpz_class RhoExpI_h = (RhoExpI * h) % p; previous
    mpz_powm(RhoExpI.get_mpz_t(), Rho.get_mpz_t(), I.get_mpz_t(), q.get_mpz_t());//Since this is an exponenet
    //mpz_class RhoExpI_h = (RhoExpI * h) % q;

    mpz_class gExp_RhoExpI;
    mpz_powm(gExp_RhoExpI.get_mpz_t(), g.get_mpz_t(), RhoExpI.get_mpz_t(), p.get_mpz_t());

    //Server beta broadcasts E_Rho
    std::pair<mpz_class, mpz_class> E_Rho = ElGamal_encrypt(Rho, pk_E);//Compute E(\rho)

    //Client computes E_RhoExpI
    std::pair<mpz_class, mpz_class> E_RhoExpI = ElGamal_exp_ct(E_Rho, I, pk_E);//Compute E(\rho^I)
    //Client computes E_RhoExpI_h and sends to server beta
    std::pair<mpz_class, mpz_class> E_h = ElGamal_encrypt(h, pk_E);//Compute E(h)
    std::pair<mpz_class, mpz_class> E_RhoExpI_h = ElGamal_mult_ct(E_RhoExpI, E_h);//Compute E(\rho^I.h)

    //Server beta decrypts and reduces to mod q
    mpz_class decrypted_RhoExpI_h = ElGamal_decrypt(E_RhoExpI_h, sk_E);
    mpz_class decrypted_RhoExpI_h__modq = decrypted_RhoExpI_h % q;//Reduce to mod q
    //Server beta raises to g
    mpz_class gExp_RhoExpI_h__modq;
    mpz_powm(gExp_RhoExpI_h__modq.get_mpz_t(), g.get_mpz_t(), decrypted_RhoExpI_h__modq.get_mpz_t(), p.get_mpz_t());
    std::pair<mpz_class, mpz_class> E_gExp_RhoExpI_h__modq = ElGamal_encrypt(gExp_RhoExpI_h__modq, pk_E);//Sends corresponding ciphertext

    //Server alpha removes h by raising it to h^{-1}
    std::pair<mpz_class, mpz_class> E_gExp_RhoExpI_hh_1__modq = ElGamal_exp_ct(E_gExp_RhoExpI_h__modq, h_1, pk_E);

    //If decrypted must match with gExp_RhoExpI
    mpz_class decrypted_gExp_RhoExpI_hh_1__modq = ElGamal_decrypt(E_gExp_RhoExpI_hh_1__modq, sk_E);

    assert(decrypted_gExp_RhoExpI_hh_1__modq == gExp_RhoExpI);

    if (decrypted_gExp_RhoExpI_hh_1__modq == gExp_RhoExpI) {
        PrintLog(LOG_LEVEL_INFO, __FILE__, __LINE__, "Server Beta: Decryption successful and matches with gExp_RhoExpI");
    } else {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Server Beta: Decryption failed or does not match with gExp_RhoExpI. Expected :" + gExp_RhoExpI.get_str() + " but got: " + decrypted_gExp_RhoExpI_hh_1__modq.get_str());
    }
    #endif
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
