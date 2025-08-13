#pragma once
#include <gmp.h>
#include <gmpxx.h>
#include <string>
#include <utility>
#include <sys/socket.h>
#include <netinet/in.h>
#include <arpa/inet.h>
#include <unistd.h>
#include <cstring>
#include <cassert>
#include <iostream>
#include <string>
#include <chrono>
#include <thread>


class ElGamal {
public:
    mpz_class p, q, g;
    void setup(int p_bits = 3072, int q_bits = 256);
    mpz_class randomGroupElement();
    std::pair<mpz_class, mpz_class> keyGen();
    std::pair<mpz_class, mpz_class> encrypt(const mpz_class& message, const mpz_class& publicKey);
    mpz_class decrypt(const std::pair<mpz_class, mpz_class>& ciphertext, const mpz_class& privateKey);
    std::pair<mpz_class, mpz_class> mult_ct(const std::pair<mpz_class, mpz_class>& ciphertext1, const std::pair<mpz_class, mpz_class>& ciphertext2);
    std::pair<mpz_class, mpz_class> exp_ct(const std::pair<mpz_class, mpz_class>& ciphertext, const mpz_class& exp, const mpz_class& publicKey);
};

// Implementation
inline void ElGamal::setup(int p_bits, int q_bits) {
    gmp_randclass rng(gmp_randinit_default);
    rng.seed(time(NULL));
    do {
        q = rng.get_z_bits(q_bits);
        mpz_nextprime(q.get_mpz_t(), q.get_mpz_t());
    } while (mpz_sizeinbase(q.get_mpz_t(), 2) != q_bits);
    mpz_class temp;
    do {
        temp = rng.get_z_bits(p_bits - q_bits);
        p = temp * q + 1;
    } while (!mpz_probab_prime_p(p.get_mpz_t(), 25));
    mpz_class h, exp;
    exp = (p - 1) / q;
    do {
        h = rng.get_z_range(p - 1) + 1;
        mpz_powm(g.get_mpz_t(), h.get_mpz_t(), exp.get_mpz_t(), p.get_mpz_t());
    } while (g == 1);
}
inline mpz_class ElGamal::randomGroupElement() {
    gmp_randclass rng(gmp_randinit_default);
    mpz_class r = rng.get_z_range(q);
    mpz_class result;
    mpz_powm(result.get_mpz_t(), g.get_mpz_t(), r.get_mpz_t(), p.get_mpz_t());
    return result;
}
inline std::pair<mpz_class, mpz_class> ElGamal::keyGen() {
    gmp_randclass rng(gmp_randinit_default);
    mpz_class x = rng.get_z_range(q);
    mpz_class y;
    mpz_powm(y.get_mpz_t(), g.get_mpz_t(), x.get_mpz_t(), p.get_mpz_t());
    return std::make_pair(x, y);
}
inline std::pair<mpz_class, mpz_class> ElGamal::encrypt(const mpz_class& message, const mpz_class& publicKey) {
    gmp_randclass rng(gmp_randinit_default);
    mpz_class k = rng.get_z_range(q);
    mpz_class c1, c2;
    mpz_powm(c1.get_mpz_t(), g.get_mpz_t(), k.get_mpz_t(), p.get_mpz_t());
    mpz_class temp;
    mpz_powm(temp.get_mpz_t(), publicKey.get_mpz_t(), k.get_mpz_t(), p.get_mpz_t());
    c2 = (message * temp) % p;
    return std::make_pair(c1, c2);
}
inline mpz_class ElGamal::decrypt(const std::pair<mpz_class, mpz_class>& ciphertext, const mpz_class& privateKey) {
    mpz_class c1 = ciphertext.first;
    mpz_class c2 = ciphertext.second;
    mpz_class temp, inv_temp;
    mpz_powm(temp.get_mpz_t(), c1.get_mpz_t(), privateKey.get_mpz_t(), p.get_mpz_t());
    mpz_invert(inv_temp.get_mpz_t(), temp.get_mpz_t(), p.get_mpz_t());
    return (c2 * inv_temp) % p;
}
inline std::pair<mpz_class, mpz_class> ElGamal::mult_ct(const std::pair<mpz_class, mpz_class>& ciphertext1, const std::pair<mpz_class, mpz_class>& ciphertext2) {
    mpz_class cm1 = (ciphertext1.first * ciphertext2.first) % p;
    mpz_class cm2 = (ciphertext1.second * ciphertext2.second) % p;
    return std::make_pair(cm1, cm2);
}
inline std::pair<mpz_class, mpz_class> ElGamal::exp_ct(const std::pair<mpz_class, mpz_class>& ciphertext, const mpz_class& exp, const mpz_class& publicKey) {
    mpz_class c1, c2;
    mpz_powm(c1.get_mpz_t(), ciphertext.first.get_mpz_t(), exp.get_mpz_t(), p.get_mpz_t());
    mpz_powm(c2.get_mpz_t(), ciphertext.second.get_mpz_t(), exp.get_mpz_t(), p.get_mpz_t());
    auto [cI1, cI2] = encrypt(mpz_class(1), publicKey);
    return mult_ct({c1, c2}, {cI1, cI2});
}

// Networking related functions
void FinishAcceptingSocket(int server_fd, int new_socket);
int InitAcceptingSocket(int port, int* p_server_fd, int* p_new_socket);
void InitConnectingSocket(const std::string& server_ip, int port, int* p_sock);//No corresponding finish function, only call close()
