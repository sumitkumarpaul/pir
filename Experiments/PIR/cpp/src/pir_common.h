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


#define N 16 // Number of elements in the plaintext database

// Global ElGamal parameters
extern mpz_class p, q, r, g;

// ElGamal cryptographic functions using global parameters
mpz_class ElGamal_randomGroupElement();
std::pair<mpz_class, mpz_class> ElGamal_keyGen();
std::pair<mpz_class, mpz_class> ElGamal_encrypt(const mpz_class& message, const mpz_class& publicKey);
mpz_class ElGamal_decrypt(const std::pair<mpz_class, mpz_class>& ciphertext, const mpz_class& privateKey);
std::pair<mpz_class, mpz_class> ElGamal_mult_ct(const std::pair<mpz_class, mpz_class>& ciphertext1, const std::pair<mpz_class, mpz_class>& ciphertext2);
std::pair<mpz_class, mpz_class> ElGamal_exp_ct(const std::pair<mpz_class, mpz_class>& ciphertext, const mpz_class& exp, const mpz_class& publicKey);

// Networking related functions
void FinishAcceptingSocket(int server_fd, int new_socket);
int InitAcceptingSocket(int port, int* p_server_fd, int* p_new_socket);
void InitConnectingSocket(const std::string& server_ip, int port, int* p_sock);//No corresponding finish function, only call close()
