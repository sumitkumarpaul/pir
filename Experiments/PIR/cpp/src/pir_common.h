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
#include <algorithm> // Required for std::shuffle
#include <string>
#include <chrono>
#include <thread>
#include <random>
#include <ctime>
#include <iomanip>
#include <fstream>
#include <filesystem> 
#include <iterator>
#include <openssl/sha.h>
#include <openfhe.h>
#include <kuku/kuku.h>
#if defined(__AVX2__) || defined(__AVX512F__)
  #include <immintrin.h>
#endif
#include <omp.h>

#include "fss/fss-common.h"
#include "fss/fss-server.h"
#include "fss/fss-client.h"

#include "ciphertext-ser.h"
#include "cryptocontext-ser.h"
#include "key/key-ser.h"
#include "scheme/bgvrns/bgvrns-ser.h"

using namespace lbcrypto;
using namespace kuku;


#define LOG_LEVEL_SPECIAL   0
#define LOG_LEVEL_ERROR     1
#define LOG_LEVEL_INFO      2
#define LOG_LEVEL_DEBUG     3
#define LOG_LEVEL_TRACE     4

#define SET_LOG_LEVEL LOG_LEVEL_TRACE

#define NET_BUF_SZ  12000000 //Size of the buffer used during transferring data over network

//#define SERVER_ALPHA_IP "192.168.16.126" // IP address of the server_alpha
//#define SERVER_BETA_IP  "192.168.16.34"//"192.168.16.245" // IP address of the server_beta
//#define SERVER_GAMMA_IP "192.168.16.132" // IP address of the server_gamma

#define SERVER_ALPHA_IP "127.0.0.1" // IP address of the server_alpha
#define SERVER_BETA_IP  "127.0.0.1"//"192.168.16.245" // IP address of the server_beta
#define SERVER_GAMMA_IP "127.0.0.1" // For the time being

#define BETA_LISTENING_TO_ALPHA_PORT    1234 // Port of beta to listen to alpha
#define BETA_LISTENING_TO_GAMMA_PORT    1235 // Port of beta to listen to gamma
#define GAMMA_LISTENING_TO_ALPHA_PORT   1236 // Port of gamma to listen to alpha

#define ALPHA_LISTENING_TO_CLIENT_PORT  1237 // Port of alpha to listen to the client
#define BETA_LISTENING_TO_CLIENT_PORT   1238 // Port of beta to listen to the client
#define GAMMA_LISTENING_TO_CLIENT_PORT  1239 // Port of gamma to listen to the client

extern std::string start_reinit_for_epoch_message;
extern std::string completed_reinit_for_epoch_message;

/* We will be experimenting with 100GB database. Each block is of size 512-bits. */
#define N       256//1677721600 // Number of elements in the plaintext database ((100*1024*1024*1024) / (512/8)) 
#define log_N   31    // ceil((log2(N)))
#define sqrt_N  16//40960//1024//40960 // ceil((sqrt(N))) TODO: Forcefully makig it 0, so that total size remains small and divisible by 16(number of cpu cores)

/* For quick-tag generation and experimentation with cuckoo hashing reduced the size of bits */
#define P_BITS  3072//128//5 // Size of p in bits
#define Q_BITS  (256+2)//(64+2)//3 // Size of q in bits. Adding two more bits, since we require to generate a safe prime and (q-1)/2 must be a prime
#define R_BITS  64//2 // Size of r in bits

/* Regarding optimization */
#define REDUCE_CT_SIZE  (1)
#warning The effect of compression is not verified yet

/* Regarding testing */
#define TEST_SHELTER_FOUND (0) /* Forcefully make the shelter search successful */

/*****************************************************************
* Since, the plaintext modulus is 65537, hence upto 16-bit number
* can be represented in a single ciphertext. However, we may add
* two ciphertexts as well and the result must be within 16-bit.
* Hence, each individual plaintext must remain within 15-bit.
* ***************************************************************/
#define PLAINTEXT_PIR_BLOCK_DATA_SIZE       512 // Number of bits in a single PIR block
#define PLAINTEXT_FHE_BLOCK_SIZE            14//Single encryptable plaintext block size is these many bits. Actually 16-bits, but after experimentation it is found that, homomorphic addition is not working for larger than 14-bit values
#define NUM_BYTES_PER_PDB_ELEMENT           ((PLAINTEXT_PIR_BLOCK_DATA_SIZE + 7) / 8) // Only data, no-index. To ensure ceiling value
#define NUM_BYTES_PER_SDB_ELEMENT           ((PLAINTEXT_PIR_BLOCK_DATA_SIZE + log_N + 7) / 8) // Data and index. To ensure ceiling value
#define TOTAL_NUM_FHE_BLOCKS_PER_ELEMENT    ((PLAINTEXT_PIR_BLOCK_DATA_SIZE + log_N + PLAINTEXT_FHE_BLOCK_SIZE - 1)/PLAINTEXT_FHE_BLOCK_SIZE)
#define NUM_FHE_BLOCKS_PER_TAG              ((P_BITS + PLAINTEXT_FHE_BLOCK_SIZE - 1) / PLAINTEXT_FHE_BLOCK_SIZE) // To ensure the ceiling value

/* Cuckoo hashing parameters */
#define CUCKOO_TABLE_SIZE       ((N+sqrt_N)*1.28) // Deliberately keeping it larger than M
#define CUCKOO_LOC_FUNC_COUNT   8 // This value already worked
#define CUCKOO_LOC_FUNC_SEED    make_item(0x1234567812345678, 0x1234567812345678) // For the time being using a fixed seed
#define CUCKOO_EMPTY_ITEM       make_item(0, 0)
#define CUCKOO_MAX_PROBE        1000000 // TODO, may be we can reduce it by noting down the maximum probe required during insertion
#define CUCKOO_STASH_SIZE       100 // It actually does not matter, since observed that it is not used


#define START_TIME_MEASUREMENT(t)    t = std::chrono::high_resolution_clock::now()
#define STOP_TIME_MEASUREMENT(t)    std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - t).count()

// Randoms
extern gmp_randclass rng;

// Global ElGamal parameters
extern mpz_class p, q, r, g, g_q, Rho;
extern std::pair<mpz_class, mpz_class> E_q_Rho;

extern mpz_class pk_E, sk_E, pk_E_q, sk_E_q;

extern PublicKey<DCRTPoly> pk_F;
extern PrivateKey<DCRTPoly> sk_F;
extern CryptoContext<DCRTPoly> FHEcryptoContext;
extern Ciphertext<DCRTPoly> vectorOnesforElement_ct;
extern Ciphertext<DCRTPoly> vectorOnesforTag_ct;
extern Ciphertext<DCRTPoly> fnd_ct;
extern Ciphertext<DCRTPoly> fnd_ct_element;
extern Ciphertext<DCRTPoly> fnd_ct_tag;

typedef struct [[gnu::packed]]{
    unsigned char element[NUM_BYTES_PER_PDB_ELEMENT];//Only data, no index
} plain_db_entry;

typedef struct [[gnu::packed]]{
    //char T[(P_BITS/8)];// P_BITS will require (P_BITS/8) bytes
    item_type cuckoo_key; // Instead of storing the full tag, only storing the 128-bit cuckoo hash key, generated by issuing SHA256 hash on the tag and then convert it to 128-bit key
    unsigned char element[NUM_BYTES_PER_SDB_ELEMENT];//Data and index
} shuffled_db_entry;

typedef struct {
    /* Here only storing the tags. The ciphertext part is stored in the RAM, in serialized format */
    mpz_class element_FHE_ct;
    mpz_class tag;
    mpz_class tag_short;
} shelter_element;

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
extern int InitListeningSocket(int port, int* p_server_fd);//No corresponding finish function, only call close()
extern int recvAll(int sock, char* data, size_t max_sz, size_t* received_sz);
extern int sendAll(int sock, const char* data, size_t sz);

// FHE related functions
extern int FHE_keyGen();
extern Ciphertext<DCRTPoly> FHE_Enc_SDBElement(const mpz_class block_content_and_index);
extern Ciphertext<DCRTPoly> FHE_bitwise_XOR(const Ciphertext<DCRTPoly>& A_ct, const Ciphertext<DCRTPoly>& B_ct);
extern void FHE_Dec_SDBElement(const Ciphertext<DCRTPoly>& ct, mpz_class& block_content_and_index);
extern void FHE_Dec_Tag(const Ciphertext<DCRTPoly>& ct, mpz_class& tag);
extern Ciphertext<DCRTPoly> FHE_Enc_Tag(const mpz_class tag);
extern Ciphertext<DCRTPoly> FHE_SelectElement(const Ciphertext<DCRTPoly>& fnd_ct, const Ciphertext<DCRTPoly>& A_ct, const Ciphertext<DCRTPoly>& B_ct);
extern Ciphertext<DCRTPoly> FHE_SelectTag(const Ciphertext<DCRTPoly>& fnd_ct, const Ciphertext<DCRTPoly>& A_ct, const Ciphertext<DCRTPoly>& B_ct);
extern void FHE_EncOfOnes(Ciphertext<DCRTPoly>& OnesforElement_ct, Ciphertext<DCRTPoly>& OnesforTag_ct);
extern void FHE_EncOfZeros(Ciphertext<DCRTPoly>& ZerosforElement_ct, Ciphertext<DCRTPoly>& ZerosforTag_ct);
extern void FHE_EncOfOnes(Ciphertext<DCRTPoly>& OnesforTag_ct);
extern void FHE_EncOfZeros(Ciphertext<DCRTPoly>& ZerosforTag_ct);
extern mpz_class import_from_file_to_mpz_class(const std::string& filename);
extern void export_to_file_from_mpz_class(const std::string& filename, const mpz_class& value);

// Database related functions
extern void read_sdb_entry(std::fstream& sdb, uint64_t id, shuffled_db_entry& out_entry);
extern void insert_sdb_entry(std::fstream& sdb, uint64_t id, const shuffled_db_entry& entry);
extern void read_pdb_entry(std::fstream& pdb, uint64_t id, plain_db_entry& out_entry);
extern void insert_pdb_entry(std::fstream& pdb, uint64_t id, const plain_db_entry& entry);

extern void convert_buf_to_item_type(const unsigned char* buf, size_t buf_size, item_type& out_item);
extern void convert_buf_to_item_type1(const unsigned char* buf, size_t buf_size, std::array<unsigned char, 16>& out_item);
extern void convert_buf_to_item_type2(const unsigned char* buf, size_t buf_size, std::array<unsigned char, 16>& out_item);

// FSS helper functions
extern size_t deserializeFssAndServerKeyEq(const char* buff, size_t buff_size, Fss& fss, ServerKeyEq& key);
extern size_t serializeFssAndServerKeyEq(const Fss& fss, const ServerKeyEq& key, char* buff, size_t buff_size);

extern bool save_mpz_vector(const std::vector<mpz_class>& vec, const std::string& path);
extern bool load_mpz_vector(std::vector<mpz_class>& vec, const std::string& path);
