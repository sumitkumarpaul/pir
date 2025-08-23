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
#include <random>
#include <ctime>
#include <iomanip>
#include <fstream>
#include <iterator>
#include <openfhe.h>

#include "ciphertext-ser.h"
#include "cryptocontext-ser.h"
#include "key/key-ser.h"
#include "scheme/bgvrns/bgvrns-ser.h"

using namespace lbcrypto;

#define LOG_LEVEL_SPECIAL   0
#define LOG_LEVEL_ERROR     1
#define LOG_LEVEL_INFO      2
#define LOG_LEVEL_DEBUG     3
#define LOG_LEVEL_TRACE     4

#define SET_LOG_LEVEL LOG_LEVEL_TRACE

#define NET_BUF_SZ  8388608 //Size of the buffer used during transferring data over network

#define SERVER_ALPHA_IP "192.168.16.246" // IP address of the server_alpha
#define SERVER_BETA_IP  "127.0.0.1"//"192.168.16.245" // IP address of the server_beta
#define SERVER_GAMMA_IP "192.168.16.244" // IP address of the server_gamma

#define BETA_LISTENING_TO_ALPHA_PORT    1234 // Port to listen to alpha
#define BETA_LISTENING_TO_GAMMA_PORT    1235 // Port to listen to gamma
#define BETA_LISTENING_TO_CLIENT_PORT   1236 // Port to listen to client

#define N 5 // Number of elements in the plaintext database

/*****************************************************************
* Since, the plaintext modulus is 65537, hence upto 16-bit number
* can be represented in a single ciphertext. However, we may add
* two ciphertexts as well and the result must be within 16-bit.
* Hence, each individual plaintext must remain within 15-bit.
* ***************************************************************/
#define PLAINTEXT_PIR_BLOCK_SIZE 512 // Number of bits in a single PIR block
#define PLAINTEXT_FHE_BLOCK_SIZE 15 // Single encryptable plaintext block size is these many bits
#define NUM_FHE_BLOCKS_PER_PIR_BLOCK ((PLAINTEXT_PIR_BLOCK_SIZE + PLAINTEXT_FHE_BLOCK_SIZE - 1) / PLAINTEXT_FHE_BLOCK_SIZE) // To ensure the ceiling value

#define P_BITS  3072//5 // Size of p in bits
#define Q_BITS  (256+2)//3 // Size of q in bits. Adding two more bits, since we require to generate a safe prime and (q-1)/2 must be a prime
#define R_BITS  64//2 // Size of r in bits

// Randoms
extern gmp_randclass rng;


// Global ElGamal parameters
extern mpz_class p, q, r, g, g_q;

extern mpz_class pk_E, sk_E, pk_E_q, sk_E_q;

extern PublicKey<DCRTPoly> pk_F;
extern PrivateKey<DCRTPoly> sk_F;
extern CryptoContext<DCRTPoly> FHEcryptoContext;

//Function prototypes
//Logging related
void PrintLog(int log_level, const char* file, int line, const std::string& message);

// ElGamal cryptographic functions using global parameters
mpz_class ElGamal_randomGroupElement();
std::pair<mpz_class, mpz_class> ElGamal_keyGen();
std::pair<mpz_class, mpz_class> ElGamal_encrypt(const mpz_class& message, const mpz_class& publicKey);
mpz_class ElGamal_decrypt(const std::pair<mpz_class, mpz_class>& ciphertext, const mpz_class& privateKey);
std::pair<mpz_class, mpz_class> ElGamal_mult_ct(const std::pair<mpz_class, mpz_class>& ciphertext1, const std::pair<mpz_class, mpz_class>& ciphertext2);
std::pair<mpz_class, mpz_class> ElGamal_exp_ct(const std::pair<mpz_class, mpz_class>& ciphertext, const mpz_class& exp, const mpz_class& publicKey);

std::pair<mpz_class, mpz_class> ElGamal_q_keyGen();
std::pair<mpz_class, mpz_class> ElGamal_q_encrypt(const mpz_class& message, const mpz_class& publicKey);
mpz_class ElGamal_q_decrypt(const std::pair<mpz_class, mpz_class>& ciphertext, const mpz_class& privateKey);
std::pair<mpz_class, mpz_class> ElGamal_q_mult_ct(const std::pair<mpz_class, mpz_class>& ciphertext1, const std::pair<mpz_class, mpz_class>& ciphertext2);
std::pair<mpz_class, mpz_class> ElGamal_q_exp_ct(const std::pair<mpz_class, mpz_class>& ciphertext, const mpz_class& exp, const mpz_class& publicKey);

// Networking related functions
void FinishAcceptingSocket(int server_fd, int new_socket);
int InitAcceptingSocket(int port, int* p_server_fd, int* p_new_socket);
void InitConnectingSocket(const std::string& server_ip, int port, int* p_sock);//No corresponding finish function, only call close()
int recvAll(int sock, char* data, size_t max_sz, size_t* received_sz);
int sendAll(int sock, const char* data, size_t sz);

// FHE related functions
int FHE_keyGen();
Ciphertext<DCRTPoly> FHE_encSingleMsg(const Plaintext& pt);
int FHE_sel(int select_bit);