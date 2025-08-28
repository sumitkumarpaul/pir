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
//#define SERVER_GAMMA_IP "192.168.16.244" // IP address of the server_gamma
#define SERVER_GAMMA_IP "127.0.0.1" // For the time being

#define BETA_LISTENING_TO_ALPHA_PORT    1234 // Port of beta to listen to alpha
#define BETA_LISTENING_TO_GAMMA_PORT    1235 // Port of beta to listen to gamma
#define BETA_LISTENING_TO_CLIENT_PORT   1236 // Port of beta to listen to client
#define GAMMA_LISTENING_TO_ALPHA_PORT   1237 // Port of gamma to listen to alpha

extern std::string ready_for_epoch_message;

/* We will be experimenting with 100GB database. Each block is of size 512-bits. */
#define N       1677721600 // Number of elements in the plaintext database ((100*1024*1024*1024) / (512/8)) 
#define log_N   31    // ceil((log2(N)))
#define sqrt_N  40960//1024//40960 // ceil((sqrt(N)))

/* For quick-tag generation and experimentation with cuckoo hashing reduced the size of bits */
#define P_BITS  3072//128//5 // Size of p in bits
#define Q_BITS  (256+2)//(64+2)//3 // Size of q in bits. Adding two more bits, since we require to generate a safe prime and (q-1)/2 must be a prime
#define R_BITS  64//2 // Size of r in bits

/*****************************************************************
* Since, the plaintext modulus is 65537, hence upto 16-bit number
* can be represented in a single ciphertext. However, we may add
* two ciphertexts as well and the result must be within 16-bit.
* Hence, each individual plaintext must remain within 15-bit.
* ***************************************************************/
#define PLAINTEXT_PIR_BLOCK_DATA_SIZE       512 // Number of bits in a single PIR block
#define PLAINTEXT_FHE_BLOCK_SIZE            15 // Single encryptable plaintext block size is these many bits. Actually 16-bits, but 1 bit is kept reserved for carry forwarding during homomorphic selection processing
#define NUM_FHE_BLOCKS_PER_PIR_BLOCK        ((PLAINTEXT_PIR_BLOCK_DATA_SIZE + PLAINTEXT_FHE_BLOCK_SIZE - 1) / PLAINTEXT_FHE_BLOCK_SIZE) // To ensure the ceiling value
#define NUM_FHE_BLOCKS_PER_PIR_INDEX        ((log_N + PLAINTEXT_FHE_BLOCK_SIZE - 1) / PLAINTEXT_FHE_BLOCK_SIZE)
#define TOTAL_NUM_FHE_BLOCKS_PER_ELEMENT    (NUM_FHE_BLOCKS_PER_PIR_BLOCK + NUM_FHE_BLOCKS_PER_PIR_INDEX)
#define NUM_FHE_BLOCKS_PER_TAG              ((P_BITS + PLAINTEXT_FHE_BLOCK_SIZE - 1) / PLAINTEXT_FHE_BLOCK_SIZE) // To ensure the ceiling value


// Randoms
extern gmp_randclass rng;

// Global ElGamal parameters
extern mpz_class p, q, r, g, g_q;

extern mpz_class pk_E, sk_E, pk_E_q, sk_E_q;

extern PublicKey<DCRTPoly> pk_F;
extern PrivateKey<DCRTPoly> sk_F;
extern CryptoContext<DCRTPoly> FHEcryptoContext;
extern Ciphertext<DCRTPoly> vectorOnesforElement_ct;
extern Ciphertext<DCRTPoly> vectorOnesforTag_ct;
extern Ciphertext<DCRTPoly> fnd_ct;

typedef struct {
    /* Here only storing the tags. The ciphertext part is stored in the RAM, in serialized format */
    mpz_class tag;
    mpz_class tag_short;
} shelter_tags;

//Function prototypes
//Logging related
extern void PrintLog(int log_level, const char* file, int line, const std::string& message);

// Utility functions
extern mpz_class import_from_bytes(const std::string &bytes);

// ElGamal cryptographic functions using global parameters
extern mpz_class ElGamal_randomGroupElement();
extern std::pair<mpz_class, mpz_class> ElGamal_keyGen();
extern std::pair<mpz_class, mpz_class> ElGamal_encrypt(const mpz_class& message, const mpz_class& publicKey);
extern mpz_class ElGamal_decrypt(const std::pair<mpz_class, mpz_class>& ciphertext, const mpz_class& privateKey);
extern std::pair<mpz_class, mpz_class> ElGamal_mult_ct(const std::pair<mpz_class, mpz_class>& ciphertext1, const std::pair<mpz_class, mpz_class>& ciphertext2);
extern std::pair<mpz_class, mpz_class> ElGamal_exp_ct(const std::pair<mpz_class, mpz_class>& ciphertext, const mpz_class& exp, const mpz_class& publicKey);

extern std::pair<mpz_class, mpz_class> ElGamal_q_keyGen();
extern std::pair<mpz_class, mpz_class> ElGamal_q_encrypt(const mpz_class& message, const mpz_class& publicKey);
extern mpz_class ElGamal_q_decrypt(const std::pair<mpz_class, mpz_class>& ciphertext, const mpz_class& privateKey);
extern std::pair<mpz_class, mpz_class> ElGamal_q_mult_ct(const std::pair<mpz_class, mpz_class>& ciphertext1, const std::pair<mpz_class, mpz_class>& ciphertext2);
extern std::pair<mpz_class, mpz_class> ElGamal_q_exp_ct(const std::pair<mpz_class, mpz_class>& ciphertext, const mpz_class& exp, const mpz_class& publicKey);

// Networking related functions
extern void FinishAcceptingSocket(int server_fd, int new_socket);
extern int InitAcceptingSocket(int port, int* p_server_fd, int* p_new_socket);
extern void InitConnectingSocket(const std::string& server_ip, int port, int* p_sock);//No corresponding finish function, only call close()
extern int recvAll(int sock, char* data, size_t max_sz, size_t* received_sz);
extern int sendAll(int sock, const char* data, size_t sz);

// FHE related functions
extern int FHE_keyGen();
extern Ciphertext<DCRTPoly> FHE_Enc_DBElement(const mpz_class block_content, const mpz_class block_index);
extern void FHE_Dec_DBElement(const Ciphertext<DCRTPoly>& ct, mpz_class& block_content, mpz_class& block_index);
extern void FHE_Dec_Tag(const Ciphertext<DCRTPoly>& ct, mpz_class& tag);
extern Ciphertext<DCRTPoly> FHE_Enc_Tag(const mpz_class tag);
extern Ciphertext<DCRTPoly> FHE_SelectElement(const Ciphertext<DCRTPoly>& fnd_ct, const Ciphertext<DCRTPoly>& A_ct, const Ciphertext<DCRTPoly>& B_ct);
extern Ciphertext<DCRTPoly> FHE_SelectTag(const Ciphertext<DCRTPoly>& fnd_ct, const Ciphertext<DCRTPoly>& A_ct, const Ciphertext<DCRTPoly>& B_ct);
extern void FHE_EncOfOnes(Ciphertext<DCRTPoly>& OnesforElement_ct, Ciphertext<DCRTPoly>& OnesforTag_ct);
extern mpz_class serialized_ct_to_mpz_class(const std::string& filename);