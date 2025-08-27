#include <assert.h>
#include <stdio.h>
#include <openssl/sha.h>
#include "fss/fss-common.h"
#include "fss/fss-server.h"
#include "fss/fss-client.h"
#include <immintrin.h>  // Include AVX header
#include <omp.h>
#include <vector>
#include <gmp.h>
#include <gmpxx.h>
#include <random>

#include <iomanip>
#include <iostream>
#include <kuku/kuku.h>


#include <iterator>
#include <cstring>
#include <unistd.h>
#include "pir_common.h"

#define CUCKOO_HASH_TABLE_REHASH_TRY_COUNT 1

static int sock_alpha_to_beta = -1, sock_alpha_to_gamma = -1;
static char net_buf[NET_BUF_SZ] = {0};

// The value of N will determine the bitlength during the client initialization
#define NUM_TAG_BITS 3072 // 16 bits can represent up to 65536, which is more than enough for N=50000
#define B 512 // Block size in bits, can be adjusted as needed
// And number of bits determine the evalution time drastically
static mpz_class sh[sqrt_N][2]; // Database to store values, each entry is a pair, {Tag, Block-content}.

// Function declarations
static int InitSrv_alpha();
static int RecvInitParamsFromBeta();
static int FinSrv_alpha();
static int SelShuffDBSearchTag_alpha(std::pair<mpz_class, mpz_class>& E_T_I);
static void TestSrv_alpha();
static void TestPKEOperations_alpha();
static void TestSelShuffDBSearchTag_alpha();

using namespace kuku;
static int Test_CuckooHash(table_size_type table_size, uint64_t num_entry, table_size_type stash_size, uint8_t loc_func_count, uint64_t max_probe);

// Function definitions
static int InitSrv_alpha(){
    int ret = -1;
    // Initialize random number generation
    std::random_device rd;
    unsigned long seed = (static_cast<unsigned long>(rd()) << 1) ^ rd();
    rng.seed(seed); // seed() seeds the gmp_randclass    
    
    // Server_alpha only connects to other servers, it does not listen
    InitConnectingSocket(SERVER_BETA_IP, BETA_LISTENING_TO_ALPHA_PORT, &sock_alpha_to_beta);
    PrintLog(LOG_LEVEL_TRACE, __FILE__, __LINE__, "Established connection with Server Beta");

    ret = RecvInitParamsFromBeta();

    if (ret != 0) {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Server Alpha: Failed to receive initialization parameters from Server Beta");
        close(sock_alpha_to_beta);
        return -1;
    } else {
        PrintLog(LOG_LEVEL_TRACE, __FILE__, __LINE__, "Server Alpha: Successfully received initialization parameters from Server Beta");
    }

    InitConnectingSocket(SERVER_GAMMA_IP, GAMMA_LISTENING_TO_ALPHA_PORT, &sock_alpha_to_gamma);

    PrintLog(LOG_LEVEL_TRACE, __FILE__, __LINE__, "Established connection with Server Gamma");

    PrintLog(LOG_LEVEL_INFO, __FILE__, __LINE__, "Server Alpha initialization complete");

    return 0;
}

static int RecvInitParamsFromBeta() {
    size_t received_sz = 0;
    int ret_recv = 0;

    // Receive p
    ret_recv = recvAll(sock_alpha_to_beta, net_buf, sizeof(net_buf), &received_sz);
    if (ret_recv != 0)
    {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Failed to receive p from Server Beta");
        return -1;
    }
    p = mpz_class(std::string(net_buf, received_sz));

    // Receive q
    ret_recv = recvAll(sock_alpha_to_beta, net_buf, sizeof(net_buf), &received_sz);
    if (ret_recv != 0)
    {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Failed to receive q from Server Beta");
        return -1;
    }
    q = mpz_class(std::string(net_buf, received_sz));

    // Receive g
    ret_recv = recvAll(sock_alpha_to_beta, net_buf, sizeof(net_buf), &received_sz);
    if (ret_recv != 0)
    {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Failed to receive g from Server Beta");
        return -1;
    }
    g = mpz_class(std::string(net_buf, received_sz));

    // Receive g_q
    ret_recv = recvAll(sock_alpha_to_beta, net_buf, sizeof(net_buf), &received_sz);
    if (ret_recv != 0)
    {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Failed to receive g_q from Server Beta");
        return -1;
    }
    g_q = mpz_class(std::string(net_buf, received_sz));

    // Receive r
    ret_recv = recvAll(sock_alpha_to_beta, net_buf, sizeof(net_buf), &received_sz);
    if (ret_recv != 0)
    {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Failed to receive r from Server Beta");
        return -1;
    }
    r = mpz_class(std::string(net_buf, received_sz));

    // Receive pk_E
    ret_recv = recvAll(sock_alpha_to_beta, net_buf, sizeof(net_buf), &received_sz);
    if (ret_recv != 0)
    {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Failed to receive pk_E from Server Beta");
        return -1;
    }
    pk_E = mpz_class(std::string(net_buf, received_sz));

    // Receive pk_E_q
    ret_recv = recvAll(sock_alpha_to_beta, net_buf, sizeof(net_buf), &received_sz);
    if (ret_recv != 0)
    {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Failed to receive pk_E_q from Server Beta");
        return -1;
    }
    pk_E_q = mpz_class(std::string(net_buf, received_sz));

    // Receive FHEcryptoContext
    ret_recv = recvAll(sock_alpha_to_beta, net_buf, sizeof(net_buf), &received_sz);
    if (ret_recv != 0)
    {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Failed to receive FHEcryptoContext from Server Beta");
        return -1;
    }
    Serial::DeserializeFromString(FHEcryptoContext, std::string(net_buf, received_sz));

    // Receive pk_F
    ret_recv = recvAll(sock_alpha_to_beta, net_buf, sizeof(net_buf), &received_sz);
    if (ret_recv != 0)
    {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Failed to receive pk_F from Server Beta");
        return -1;
    }
    Serial::DeserializeFromString(pk_F, std::string(net_buf, received_sz));

    // Receive vectorOnesforElement_ct
    ret_recv = recvAll(sock_alpha_to_beta, net_buf, sizeof(net_buf), &received_sz);
    if (ret_recv != 0)
    {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Failed to receive vectorOnesforElement_ct from Server Beta");
        return -1;
    }
    Serial::DeserializeFromString(vectorOnesforElement_ct, std::string(net_buf, received_sz));

    // Receive vectorOnesforTag_ct
    ret_recv = recvAll(sock_alpha_to_beta, net_buf, sizeof(net_buf), &received_sz);
    if (ret_recv != 0)
    {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Failed to receive vectorOnesforTag_ct from Server Beta");
        return -1;
    }
    Serial::DeserializeFromString(vectorOnesforTag_ct, std::string(net_buf, received_sz));

    return 0;
}

static int FinSrv_alpha(){
    int ret = -1;

    // Close the sockets
    if (sock_alpha_to_beta != -1) {
        close(sock_alpha_to_beta);
        sock_alpha_to_beta = -1;
    }
    if (sock_alpha_to_gamma != -1) {
        close(sock_alpha_to_gamma);
        sock_alpha_to_gamma = -1;
    }

    PrintLog(LOG_LEVEL_SPECIAL, __FILE__, __LINE__, "Finalized Server Alpha");

    return ret;
}


static int SelShuffDBSearchTag_alpha(std::pair<mpz_class, mpz_class>& E_T_I){
    int ret = 0;
    size_t received_sz = 0;
    mpz_class h_alpha1 = ElGamal_randomGroupElement();
    mpz_class h_alpha1_1;
    mpz_invert(h_alpha1_1.get_mpz_t(), h_alpha1.get_mpz_t(), p.get_mpz_t());

    mpz_class h_alpha2 = ElGamal_randomGroupElement();
    mpz_class h_alpha2_1;
    mpz_invert(h_alpha2_1.get_mpz_t(), h_alpha2.get_mpz_t(), p.get_mpz_t());

    std::pair<mpz_class, mpz_class> E_T_I_h_alpha1 = ElGamal_mult_ct(E_T_I, ElGamal_encrypt(h_alpha1, pk_E));

    /* Send E(T_I.h_{\alpha 1}) */
    (void)sendAll(sock_alpha_to_beta, E_T_I_h_alpha1.first.get_str().c_str(), E_T_I_h_alpha1.first.get_str().size());
    (void)sendAll(sock_alpha_to_beta, E_T_I_h_alpha1.second.get_str().c_str(), E_T_I_h_alpha1.second.get_str().size());

    /* Receive E(T_I.h_{\alpha 1}h_{\beta 0}) */
    (void)recvAll(sock_alpha_to_beta, net_buf, sizeof(net_buf), &received_sz);
    if (ret != 0)
    {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Failed to receive E(T_I.h_{\\alpha 1}h_{\\beta 0}) from Server Beta");
        return ret;
    }

    /* Remove h_{\alpha 1} and deternube T_I_h_beta0 */
    mpz_class T_I_h_alpha1_h_beta0 = mpz_class(std::string(net_buf, received_sz));
    mpz_class T_I_h_beta0 = (T_I_h_alpha1_h_beta0 * h_alpha1_1) % p;

    PrintLog(LOG_LEVEL_INFO, __FILE__, __LINE__, "Determined T_I.h_{\\beta 0}: " + T_I_h_beta0.get_str());

    return ret;
}

static void TestSelShuffDBSearchTag_alpha(){
    int ret = -1;
    mpz_class test_T_I = ElGamal_randomGroupElement();

    PrintLog(LOG_LEVEL_INFO, __FILE__, __LINE__, "TestSelShuffDBSearchTag_alpha: Testing with T_I: " + test_T_I.get_str());

    std::pair<mpz_class, mpz_class> E_T_I = ElGamal_encrypt(test_T_I, pk_E);

    ret = SelShuffDBSearchTag_alpha(E_T_I);

    return;
}

static void TestPKEOperations_alpha(){
    int ret;

    //Choose random message and random exponent
    mpz_class m1 = ElGamal_randomGroupElement();
    mpz_class m2 = ElGamal_randomGroupElement();
    mpz_class exp = rng.get_z_range(q);

    mpz_class m3 = (m1*m2)%p;
    mpz_class m4;
    mpz_powm(m4.get_mpz_t(), m1.get_mpz_t(), exp.get_mpz_t(), p.get_mpz_t());

    auto [c11, c12] = ElGamal_encrypt(m1, pk_E);
    auto [c21, c22] = ElGamal_encrypt(m2, pk_E);
    auto [c31, c32] = ElGamal_mult_ct({c11, c12}, {c21, c22});
    auto [c41, c42] = ElGamal_exp_ct({c11, c12}, exp, pk_E);

    // Generate random tag of P_BITS bits
    mpz_class tag = rng.get_z_bits(P_BITS);
    Ciphertext<DCRTPoly> ct_tag = FHE_Enc_Tag(tag);

    (void)sendAll(sock_alpha_to_beta, m1.get_str().c_str(), m1.get_str().size());
    (void)sendAll(sock_alpha_to_beta, m2.get_str().c_str(), m2.get_str().size());
    (void)sendAll(sock_alpha_to_beta, m3.get_str().c_str(), m3.get_str().size());
    (void)sendAll(sock_alpha_to_beta, m4.get_str().c_str(), m4.get_str().size());
    (void)sendAll(sock_alpha_to_beta, c11.get_str().c_str(), c11.get_str().size());
    (void)sendAll(sock_alpha_to_beta, c12.get_str().c_str(), c12.get_str().size());
    (void)sendAll(sock_alpha_to_beta, c21.get_str().c_str(), c21.get_str().size());
    (void)sendAll(sock_alpha_to_beta, c22.get_str().c_str(), c22.get_str().size());
    (void)sendAll(sock_alpha_to_beta, c31.get_str().c_str(), c31.get_str().size());
    (void)sendAll(sock_alpha_to_beta, c32.get_str().c_str(), c32.get_str().size());
    (void)sendAll(sock_alpha_to_beta, c41.get_str().c_str(), c41.get_str().size());
    (void)sendAll(sock_alpha_to_beta, c42.get_str().c_str(), c42.get_str().size());
    (void)sendAll(sock_alpha_to_beta, tag.get_str().c_str(), tag.get_str().size());
    (void)sendAll(sock_alpha_to_beta, Serial::SerializeToString(ct_tag).c_str(), Serial::SerializeToString(ct_tag).size());

    return;
}

static void TestSrv_alpha()
{
    //TestPKEOperations_alpha();
    TestSelShuffDBSearchTag_alpha();
}

ostream &operator<<(ostream &stream, item_type item)
{
    stream << item[1] << " " << item[0];
    return stream;
}

void print_stash(const KukuTable &table)
{
    for (table_size_type i = 0; i < table.stash().size(); i++)
    {
        const auto &item = table.stash(i);
        //PrintLog(LOG_LEVEL_INFO, __FILE__, __LINE__, "Stash item " + to_string(i) + ": " + get_high_word(item) + get_low_word(item));
    }
}

/* Copied from: https://github.com/microsoft/Kuku/blob/main/examples/example.cpp */
static int Test_CuckooHash(table_size_type table_size, uint64_t num_entry, table_size_type stash_size, uint8_t loc_func_count, uint64_t max_probe)
{
    unsigned int rehash_cnt = 0;
    KukuTable *table = nullptr;
    //unsigned long hash_table_sz = (N + sqrt_N);
    unsigned long hash_table_sz = table_size;
    item_type loc_func_seed;
    item_type empty_item = make_item(0, 0);
    uint64_t high_rand, low_rand;
    mpz_class rand_val_high, rand_val_low;
    unsigned long I;
    bool success = false;

    /* Select a random rho for experimentation */
    mpz_class local_Rho = ElGamal_randomGroupElement();

    PrintLog(LOG_LEVEL_INFO, __FILE__, __LINE__, "Starting to build Cuckoo hash table of size: " + to_string(hash_table_sz) + " , number of entries: " + to_string(num_entry) + " with stash size: " + to_string(stash_size) + " with function count: " + to_string(loc_func_count)  + " and " + to_string(max_probe) + " number of maximum probes");

    while (rehash_cnt < CUCKOO_HASH_TABLE_REHASH_TRY_COUNT) // Limit the number of rehash attempts to avoid infinite loops
    {
        rand_val_high = rng.get_z_bits(64);
        rand_val_low = rng.get_z_bits(64);

        mpz_export(&high_rand, nullptr, 1, sizeof(high_rand), 0, 0, rand_val_high.get_mpz_t());
        mpz_export(&low_rand, nullptr, 1, sizeof(low_rand), 0, 0, rand_val_low.get_mpz_t());

        loc_func_seed = make_item(high_rand, low_rand);

        // Allocate a new Cuckoo hash table, larger than the number of entries, with the hope of less number of probe
        table = new KukuTable(hash_table_sz, stash_size, loc_func_count, loc_func_seed, max_probe, empty_item);

        /* From I = 1 to N+sqrt{N} */
        for (I = 1; I < (num_entry+1); I++)
        {
            /* Create a ||P||-bit tag as per our MAC-function */
            mpz_class Rho_pow_I;
            mpz_powm(Rho_pow_I.get_mpz_t(), local_Rho.get_mpz_t(), mpz_class(I).get_mpz_t(), q.get_mpz_t());

            mpz_class T_I;
            mpz_powm(T_I.get_mpz_t(), g.get_mpz_t(), Rho_pow_I.get_mpz_t(), p.get_mpz_t());

            // Convert tag to string
            std::string T_I_str = T_I.get_str();

            // Hash with SHA-256
            unsigned char hash[SHA256_DIGEST_LENGTH];
            SHA256(reinterpret_cast<const unsigned char *>(T_I_str.data()), T_I_str.size(), hash);

            // Get first 128 bits as two 64-bit words
            uint64_t high = 0, low = 0;
            memcpy(&high, hash, 8);
            memcpy(&low, hash + 8, 8);

            // Create item_type key
            item_type key = make_item(high, low);

            // Insert into the cuckoo hash table
            if (!table->insert(key))
            {
                PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Insertion failed for the tag, having value:" + T_I.get_str());
                PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Before failure, successfully inserted: " + to_string(I) + " out of " + to_string(hash_table_sz) + " items");
                PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "The size of the stash during failure: " + to_string(table->stash().size()));
                rehash_cnt++;
                PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Starting again, rehash count: " + to_string(rehash_cnt));

                /* Delete the already built table */
                delete table;
                table = nullptr;
                break;
            }
            if ((I > 0) && (I % 100000000 == 0)) {
                PrintLog(LOG_LEVEL_INFO, __FILE__, __LINE__, "Inserted " + to_string(I) + " items so far, and current stash size is: " + to_string(table->stash().size()));
            }
        }

        if (I == (num_entry+1)) {
            PrintLog(LOG_LEVEL_INFO, __FILE__, __LINE__, "Successfully inserted all " + to_string(hash_table_sz) + " items, to the Cuckoo hash table, after: " + to_string(rehash_cnt) + " rehash attempts, and the stash size is: " + to_string(table->stash().size()));
            success = true;
            break; // Successfully inserted all items
        }
    }

    if (success == true) {
        PrintLog(LOG_LEVEL_INFO, __FILE__, __LINE__, "Cuckoo hash table built successfully, now starting to process queries");

        /* From I = 1 to N+sqrt{N} */
        for (I = 1; I < (num_entry+1); I++)
        {
            /* Create a ||P||-bit tag as per our MAC-function */
            mpz_class Rho_pow_I;
            mpz_powm(Rho_pow_I.get_mpz_t(), local_Rho.get_mpz_t(), mpz_class(I).get_mpz_t(), q.get_mpz_t());

            mpz_class T_I;
            mpz_powm(T_I.get_mpz_t(), g.get_mpz_t(), Rho_pow_I.get_mpz_t(), p.get_mpz_t());

            // Convert tag to string
            std::string T_I_str = T_I.get_str();

            // Hash with SHA-256
            unsigned char hash[SHA256_DIGEST_LENGTH];
            SHA256(reinterpret_cast<const unsigned char *>(T_I_str.data()), T_I_str.size(), hash);

            // Get first 128 bits as two 64-bit words
            uint64_t high = 0, low = 0;
            memcpy(&high, hash, 8);
            memcpy(&low, hash + 8, 8);

            // Create item_type key
            item_type key = make_item(high, low);

            QueryResult res = table->query(key);
            
            if (!res)
            {
                PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Query failed for the tag, having value:" + T_I.get_str());
            }
            if ((I > 0) && (I % 100000000 == 0)) {
                PrintLog(LOG_LEVEL_INFO, __FILE__, __LINE__, "Completed " + to_string(I) + " queries so far, and current stash size is: " + to_string(table->stash().size()));
            }
        }
    } else {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Failed to build Cuckoo hash table");
    }

    PrintLog(LOG_LEVEL_INFO, __FILE__, __LINE__, "Completed all the query processing");

    delete table;

    return 0;
}


int main(int argc, char *argv[])
{
    InitSrv_alpha();
    #if 0
    TestSrv_alpha();

    FinSrv_alpha();
    #else

    unsigned long hash_table_sz = (N + sqrt_N);

    /* Probe upto half of the items, it will increase the build time, but not the lookup time */
    //Test_CuckooHash(hash_table_sz, (hash_table_sz +(hash_table_sz/2)), (sqrt_N/2), 32, (hash_table_sz/2));

    /* Keep decreasing the function count */
    //Test_CuckooHash(hash_table_sz, (hash_table_sz +(hash_table_sz/2)), (sqrt_N/2), 16, (hash_table_sz/2));

    //Test_CuckooHash(hash_table_sz, (hash_table_sz +(hash_table_sz/2)), (sqrt_N/2), 2, (hash_table_sz/2));


    Test_CuckooHash(atoi(argv[1]), atoi(argv[2]), atoi(argv[3]), atoi(argv[4]), atoi(argv[5]));
    #endif


    return 0;
}