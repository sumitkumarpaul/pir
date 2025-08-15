#include <gmp.h>
#include <gmpxx.h>
#include <ctime>

// Global ElGamal parameters
mpz_class p, q, r, g;

mpz_class ElGamal_randomGroupElement() {
    gmp_randclass rng(gmp_randinit_default);
    mpz_class r = rng.get_z_range(q);
    mpz_class result;
    mpz_powm(result.get_mpz_t(), g.get_mpz_t(), r.get_mpz_t(), p.get_mpz_t());
    return result;
}

std::pair<mpz_class, mpz_class> ElGamal_keyGen() {
    gmp_randclass rng(gmp_randinit_default);
    mpz_class x = rng.get_z_range(q);
    mpz_class y;
    mpz_powm(y.get_mpz_t(), g.get_mpz_t(), x.get_mpz_t(), p.get_mpz_t());
    return std::make_pair(x, y);
}

std::pair<mpz_class, mpz_class> ElGamal_encrypt(const mpz_class& message, const mpz_class& publicKey) {
    gmp_randclass rng(gmp_randinit_default);
    mpz_class k = rng.get_z_range(q);
    mpz_class c1, c2;
    mpz_powm(c1.get_mpz_t(), g.get_mpz_t(), k.get_mpz_t(), p.get_mpz_t());
    mpz_class temp;
    mpz_powm(temp.get_mpz_t(), publicKey.get_mpz_t(), k.get_mpz_t(), p.get_mpz_t());
    c2 = (message * temp) % p;
    return std::make_pair(c1, c2);
}

mpz_class ElGamal_decrypt(const std::pair<mpz_class, mpz_class>& ciphertext, const mpz_class& privateKey) {
    mpz_class c1 = ciphertext.first;
    mpz_class c2 = ciphertext.second;
    mpz_class temp, inv_temp;
    mpz_powm(temp.get_mpz_t(), c1.get_mpz_t(), privateKey.get_mpz_t(), p.get_mpz_t());
    mpz_invert(inv_temp.get_mpz_t(), temp.get_mpz_t(), p.get_mpz_t());
    return (c2 * inv_temp) % p;
}

std::pair<mpz_class, mpz_class> ElGamal_mult_ct(const std::pair<mpz_class, mpz_class>& ciphertext1, const std::pair<mpz_class, mpz_class>& ciphertext2) {
    mpz_class cm1 = (ciphertext1.first * ciphertext2.first) % p;
    mpz_class cm2 = (ciphertext1.second * ciphertext2.second) % p;
    return std::make_pair(cm1, cm2);
}

std::pair<mpz_class, mpz_class> ElGamal_exp_ct(const std::pair<mpz_class, mpz_class>& ciphertext, const mpz_class& exp, const mpz_class& publicKey) {
    mpz_class c1, c2;
    mpz_powm(c1.get_mpz_t(), ciphertext.first.get_mpz_t(), exp.get_mpz_t(), p.get_mpz_t());
    mpz_powm(c2.get_mpz_t(), ciphertext.second.get_mpz_t(), exp.get_mpz_t(), p.get_mpz_t());
    auto [cI1, cI2] = ElGamal_encrypt(mpz_class(1), publicKey);
    return ElGamal_mult_ct({c1, c2}, {cI1, cI2});
}
#include "pir_common.h"

int InitAcceptingSocket(int port, int* p_server_fd, int* p_new_socket) {
    struct sockaddr_in address;
    int opt = 1;
    int addrlen = sizeof(address);
    char buffer[65536] = {0};
    *p_new_socket = -1;
    *p_server_fd = -1;
    int ret = -1;

    *p_server_fd = socket(AF_INET, SOCK_STREAM, 0);
    if (*p_server_fd < 0) {
        std::cerr << "Server: Socket creation failed" << std::endl;
        return ret;
    }
    if (setsockopt(*p_server_fd, SOL_SOCKET, SO_REUSEADDR | SO_REUSEPORT, &opt, sizeof(opt)) < 0) {
        std::cerr << "Server: setsockopt failed" << std::endl;
        close(*p_server_fd);
        *p_server_fd = -1;
        return ret;
    }
    address.sin_family = AF_INET;
    address.sin_addr.s_addr = INADDR_ANY;
    address.sin_port = htons(port);
    if (bind(*p_server_fd, (struct sockaddr *)&address, sizeof(address)) < 0) {
        std::cerr << "Server: Bind failed" << std::endl;
        close(*p_server_fd);
        *p_server_fd = -1;
        return ret;
    }
    if (listen(*p_server_fd, 1) < 0) {
        std::cerr << "Server: Listen failed" << std::endl;
        close(*p_server_fd);
        *p_server_fd = -1;
        return ret;
    }
    *p_new_socket = accept(*p_server_fd, (struct sockaddr *)&address, (socklen_t*)&addrlen);
    if (*p_new_socket < 0) {
        std::cerr << "Server: Accept failed" << std::endl;
        close(*p_server_fd);
        *p_server_fd = -1;
        return ret;
    }

    return 0;
}

void FinishAcceptingSocket(int server_fd, int new_socket) {
    close(new_socket);
    close(server_fd);
}

void InitConnectingSocket(const std::string& server_ip, int port, int* p_sock) {
    // Connect to server
    struct sockaddr_in serv_addr;
    *p_sock = socket(AF_INET, SOCK_STREAM, 0);
    serv_addr.sin_family = AF_INET;
    serv_addr.sin_port = htons(port);
    inet_pton(AF_INET, server_ip.c_str(), &serv_addr.sin_addr);
    while (connect(*p_sock, (struct sockaddr *)&serv_addr, sizeof(serv_addr)) < 0) {
        std::this_thread::sleep_for(std::chrono::milliseconds(100)); // Wait for server
    }

    return;
}