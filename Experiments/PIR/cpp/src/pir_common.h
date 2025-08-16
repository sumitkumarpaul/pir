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


#define SERVER_ALPHA_IP "192.168.16.246" // IP address of the server_alpha
#define SERVER_BETA_IP  "127.0.0.1"//"192.168.16.245" // IP address of the server_beta
#define SERVER_GAMMA_IP "192.168.16.244" // IP address of the server_gamma

#define BETA_LISTENING_TO_ALPHA_PORT    1234 // Port to listen to alpha
#define BETA_LISTENING_TO_GAMMA_PORT    1235 // Port to listen to gamma
#define BETA_LISTENING_TO_CLIENT_PORT   1236 // Port to listen to client

#define N 16 // Number of elements in the plaintext database


#define P_BITS  5 // Size of p in bits
#define Q_BITS  3 // Size of q in bits
#define R_BITS  2 // Size of r in bits


// Global ElGamal parameters
extern mpz_class p, q, r, g;

extern mpz_class pk_E, sk_E;

// ElGamal cryptographic functions using global parameters
mpz_class ElGamal_randomGroupElement();
std::pair<mpz_class, mpz_class> ElGamal_keyGen(const mpz_class& p, const mpz_class& q, const mpz_class& g);
std::pair<mpz_class, mpz_class> ElGamal_encrypt(const mpz_class& message, const mpz_class& publicKey);
mpz_class ElGamal_decrypt(const std::pair<mpz_class, mpz_class>& ciphertext, const mpz_class& privateKey);
std::pair<mpz_class, mpz_class> ElGamal_mult_ct(const std::pair<mpz_class, mpz_class>& ciphertext1, const std::pair<mpz_class, mpz_class>& ciphertext2);
std::pair<mpz_class, mpz_class> ElGamal_exp_ct(const std::pair<mpz_class, mpz_class>& ciphertext, const mpz_class& exp, const mpz_class& publicKey);

// Networking related functions
void FinishAcceptingSocket(int server_fd, int new_socket);
int InitAcceptingSocket(int port, int* p_server_fd, int* p_new_socket);
void InitConnectingSocket(const std::string& server_ip, int port, int* p_sock);//No corresponding finish function, only call close()
