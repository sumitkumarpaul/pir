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


#include <iterator>
#include <cstring>
#include <unistd.h>
#include "pir_common.h"

static int sock_alpha_to_beta = -1, sock_alpha_to_gamma = -1;
static char net_buf[NET_BUF_SZ] = {0};

// The value of N will determine the bitlength during the client initialization
#define NUM_TAG_BITS 3072 // 16 bits can represent up to 65536, which is more than enough for N=50000
#define B 512 // Block size in bits, can be adjusted as needed
// And number of bits determine the evalution time drastically
static shelter_element sh[sqrt_N]; // Database to store values, each entry is a tuple.

std::pair<mpz_class, mpz_class> E_T_I;

#define CUCKOO_HASH_TABLE_REHASH_TRY_COUNT 1
#define NUM_CPU_CORES 16

#define ONE_TIME_MATERIALS_LOCATION_ALPHA std::string("/mnt/sumit/PIR_ALPHA/ONE_TIME_MATERIALS/")
#define PER_EPOCH_MATERIALS_LOCATION_ALPHA std::string("/mnt/sumit/PIR_ALPHA/PER_EPOCH_MATERIALS/")
#define DATABASE_LOCATION_ALPHA std::string("/mnt/sumit/PIR_ALPHA/")
std::string sdb_filename = DATABASE_LOCATION_ALPHA+"ShuffledDB_alpha.bin";


// Function declarations
static int InitSrv_alpha();
static int OneTimeInit_alpha();
static int FinSrv_alpha();
static int SelShuffDBSearchTag_alpha();
static int PerEpochReInit_alpha();


static void TestSrv_alpha();
static void TestPKEOperations_alpha();
static void TestSelShuffDBSearchTag_alpha();
static int TestShelterDPFSearch_alpha();
static int TestClientProcessing_alpha();

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

    InitConnectingSocket(SERVER_GAMMA_IP, GAMMA_LISTENING_TO_ALPHA_PORT, &sock_alpha_to_gamma);

    PrintLog(LOG_LEVEL_TRACE, __FILE__, __LINE__, "Established connection with Server Gamma");

    PrintLog(LOG_LEVEL_INFO, __FILE__, __LINE__, "Server Alpha initialization complete");

    ret = 0;

    return ret;
}

static int OneTimeInit_alpha() {
    size_t received_sz = 0;
    int ret_recv = 0;

    // Receive all the parameters from server beta
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

    //Save parameters to local files
    export_to_file_from_mpz_class(ONE_TIME_MATERIALS_LOCATION_ALPHA + "p.bin", p);
    export_to_file_from_mpz_class(ONE_TIME_MATERIALS_LOCATION_ALPHA + "q.bin", q);
    export_to_file_from_mpz_class(ONE_TIME_MATERIALS_LOCATION_ALPHA + "g.bin", g);
    export_to_file_from_mpz_class(ONE_TIME_MATERIALS_LOCATION_ALPHA + "g_q.bin", g_q);
    export_to_file_from_mpz_class(ONE_TIME_MATERIALS_LOCATION_ALPHA + "r.bin", r);
    export_to_file_from_mpz_class(ONE_TIME_MATERIALS_LOCATION_ALPHA + "pk_E.bin", pk_E);
    export_to_file_from_mpz_class(ONE_TIME_MATERIALS_LOCATION_ALPHA + "pk_E_q.bin", pk_E_q);

    Serial::SerializeToFile(ONE_TIME_MATERIALS_LOCATION_ALPHA + "FHEcryptoContext.bin", FHEcryptoContext, SerType::BINARY);
    Serial::SerializeToFile(ONE_TIME_MATERIALS_LOCATION_ALPHA + "pk_F.bin", pk_F, SerType::BINARY);
    Serial::SerializeToFile(ONE_TIME_MATERIALS_LOCATION_ALPHA + "vectorOnesforElement_ct.bin", vectorOnesforElement_ct, SerType::BINARY);
    Serial::SerializeToFile(ONE_TIME_MATERIALS_LOCATION_ALPHA + "vectorOnesforTag_ct.bin", vectorOnesforTag_ct, SerType::BINARY);

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

static int PerEpochReInit_alpha(){
    int ret = 0;
    mpz_class T_I, d_share_alpha;
    size_t received_sz = 0;
    uint64_t M = N + sqrt_N;

    /* Receive ready message from server-beta */
    (void)recvAll(sock_alpha_to_beta, net_buf, sizeof(net_buf), &received_sz);
    if (ret != 0)
    {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Failed to receive START_REINIT_FOR_EPOCH message from Server Beta");
        return ret;
    }

    if (std::string(net_buf, received_sz) != start_reinit_for_epoch_message) {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Did not receive expected START_REINIT_FOR_EPOCH message from Server Beta");
        return -1;
    }

    for (uint64_t i = 1; i <= M; ++i){
        /* Receive the tag */
        ret = recvAll(sock_alpha_to_beta, net_buf, sizeof(net_buf), &received_sz);
        if (ret == 0)
        {
            mpz_import(T_I.get_mpz_t(), received_sz, 1, 1, 1, 0, net_buf);
        } else {
            PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Failed to receive tag from Server Beta");
            return ret;
        }

        /* Receive the secret-share */
        ret = recvAll(sock_alpha_to_beta, net_buf, sizeof(net_buf), &received_sz);
        if (ret == 0)
        {
            mpz_import(d_share_alpha.get_mpz_t(), received_sz, 1, 1, 1, 0, net_buf);
        } else {
            PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Failed to receive tag from Server Beta");
            return ret;
        }
    }

    /* Receive ready message from server-beta */
    (void)recvAll(sock_alpha_to_beta, net_buf, sizeof(net_buf), &received_sz);
    if (ret != 0)
    {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Failed to receive READY_FOR_EPOCH message from Server Beta");
        return ret;
    }

    if (std::string(net_buf, received_sz) != completed_reinit_for_epoch_message) {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Did not receive expected COMPLETED_REINIT_FOR_EPOCH message from Server Beta");
        return -1;
    } else {
        PrintLog(LOG_LEVEL_TRACE, __FILE__, __LINE__, "Server Alpha: Ready for processing client requests");
    }

    return ret;
}

static int SelShuffDBSearchTag_alpha(){
    int ret = 0;
    size_t received_sz = 0;

    PrintLog(LOG_LEVEL_TRACE, __FILE__, __LINE__, "Server Alpha: Starting SelShuffDBSearchTag sequence");

    /* 1.a.1.1 Select random h_{alpha1} */
    mpz_class h_alpha1 = ElGamal_randomGroupElement();
    /* 1.a.1.2 Also find its inverse */
    mpz_class h_alpha1_1;
    mpz_invert(h_alpha1_1.get_mpz_t(), h_alpha1.get_mpz_t(), p.get_mpz_t());

    /* 1.a.2.1 Select random h_{alpha2} */
    mpz_class h_alpha2 = ElGamal_randomGroupElement();
    /* 1.a.2.2 Also find its inverse */
    mpz_class h_alpha2_1;
    mpz_invert(h_alpha2_1.get_mpz_t(), h_alpha2.get_mpz_t(), p.get_mpz_t());

    /* 2.1 Compute E(T_I.h_{\alpha 1}) to server beta */
    std::pair<mpz_class, mpz_class> E_T_I_h_alpha1 = ElGamal_mult_ct(E_T_I, ElGamal_encrypt(h_alpha1, pk_E));

    /* 2.2 Send both the componets of E(T_I.h_{\alpha 1}) */
    (void)sendAll(sock_alpha_to_beta, E_T_I_h_alpha1.first.get_str().c_str(), E_T_I_h_alpha1.first.get_str().size());
    (void)sendAll(sock_alpha_to_beta, E_T_I_h_alpha1.second.get_str().c_str(), E_T_I_h_alpha1.second.get_str().size());

    /* 4.a.1 Receive T_I.h_{\alpha 1}h_{\beta 0} */
    (void)recvAll(sock_alpha_to_beta, net_buf, sizeof(net_buf), &received_sz);
    if (ret != 0)
    {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Failed to receive T_I.h_{\\alpha 1}h_{\\beta 0} from Server Beta");
        return ret;
    }

    /* 5. Remove h_{\alpha 1} and determine T_I_h_beta0 */
    mpz_class T_I_h_alpha1_h_beta0 = mpz_class(std::string(net_buf, received_sz));
    mpz_class T_I_h_beta0 = (T_I_h_alpha1_h_beta0 * h_alpha1_1) % p;

    //PrintLog(LOG_LEVEL_INFO, __FILE__, __LINE__, "Determined T_I.h_{\\beta 0}: " + T_I_h_beta0.get_str());

    /* 6.1 Determine T_I.h_{\\alpha 2}.h_{\\beta 0} */
    mpz_class T_I_h_alpha2_h_beta0 = (T_I_h_beta0 * h_alpha2) % p;

    /* 6.2 FHE Encrypt T_I.h_{\\alpha 2}.h_{\\beta 0} */
    Ciphertext<DCRTPoly> FHE_ct_T_I_h_alpha2_h_beta0 = FHE_Enc_Tag(T_I_h_alpha2_h_beta0);

    /* 9.a Send h_{\alpha 2} to server gamma */
    (void)sendAll(sock_alpha_to_gamma, h_alpha2.get_str().c_str(), h_alpha2.get_str().size());

    // 10.a Receive FHE Ciphertext of T_phi.h_{\\alpha 2}.h_{\\beta 0}
    ret = recvAll(sock_alpha_to_gamma, net_buf, sizeof(net_buf), &received_sz);
    if (ret != 0)
    {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Failed to receive FHE Ciphertext of T_phi.h_{\\alpha 2}.h_{\\beta 0} from Server Gamma");
        return -1;
    }
    Ciphertext<DCRTPoly> FHE_ct_T_phi_h_alpha2_h_beta0;
    Serial::DeserializeFromString(FHE_ct_T_phi_h_alpha2_h_beta0, std::string(net_buf, received_sz));

    /* 11.a.1 Homomorphically select T_*h_{\\alpha 2}h_{\\beta 0} */
    Ciphertext<DCRTPoly> FHE_ct_T_star_h_alpha2_h_beta0 = FHE_SelectTag(fnd_ct, FHE_ct_T_I_h_alpha2_h_beta0, FHE_ct_T_phi_h_alpha2_h_beta0);

    /* 11.a.2 Send FHE_ct_T_star_h_alpha2_h_beta0 to server beta for decryption */
    (void)sendAll(sock_alpha_to_beta, Serial::SerializeToString(FHE_ct_T_star_h_alpha2_h_beta0).c_str(), Serial::SerializeToString(FHE_ct_T_star_h_alpha2_h_beta0).size());

    /* 13.a Receive T_star.h_{\alpha 2} from the server beta */
    (void)recvAll(sock_alpha_to_beta, net_buf, sizeof(net_buf), &received_sz);
    if (ret != 0)
    {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Failed to receive T_star.h_{\\alpha 2} from Server Beta");
        return -1;
    }
    mpz_class T_star_h_alpha2 = mpz_class(std::string(net_buf, received_sz));

    /* 14.a.1 Extract T_* */
    mpz_class T_star = (T_star_h_alpha2 * h_alpha2_1) % p;

    /* 14.a.2 Send T_* to server gamma */
    (void)sendAll(sock_alpha_to_gamma, T_star.get_str().c_str(), T_star.get_str().size());

    return ret;
}

static void TestSelShuffDBSearchTag_alpha(){
    int ret = -1;
    mpz_class test_T_I = ElGamal_randomGroupElement();
    /* For testing purpose, either choose as ciphertext of 0s */
    //fnd_ct = FHE_Enc_Tag(mpz_class(0));
    /* Or, the ciphertext of 1 */
    fnd_ct = vectorOnesforTag_ct;

    PrintLog(LOG_LEVEL_INFO, __FILE__, __LINE__, "TestSelShuffDBSearchTag_alpha: Testing with T_I: " + test_T_I.get_str());

    /* Generate a dummy E_T_I for testing purpose */
    E_T_I = ElGamal_encrypt(test_T_I, pk_E);

    ret = SelShuffDBSearchTag_alpha();

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
    //TestSelShuffDBSearchTag_alpha();
    //TestShelterDPFSearch_alpha();
    TestClientProcessing_alpha();

    return;
}

ostream &operator<<(ostream &stream, item_type item)
{
    stream << item[1] << " " << item[0];
    return stream;
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

static int TestShelterDPFSearch_alpha() {
    // Set up variables
    Fss fClient, fServer;
    ServerKeyEq k0;
    ServerKeyEq k1;
    int ret = 0;
    size_t received_sz = 0;
    /* On average half of the shelter elements will be populated */
    int average_shelter_size = (sqrt_N/2);

    #define DPF_SEARCH_INDEX_K 6
    //#define SHELTER_STORING_LOCATION std::string("./")
    #define SHELTER_STORING_LOCATION std::string("/dev/shm/")

    // First, receive sk_F from the server Beta
    ret = recvAll(sock_alpha_to_beta, net_buf, sizeof(net_buf), &received_sz);
    if (ret != 0)
    {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Failed to receive sk_F from Server Beta");
        return -1;
    }
    Serial::DeserializeFromString(sk_F, std::string(net_buf, received_sz));

    PrintLog(LOG_LEVEL_TRACE, __FILE__, __LINE__, "Starting to randomly populate a shelter of average size: " + to_string(average_shelter_size));

    /* Populate the shelter, with random elements */
    for(size_t k = 0; k < average_shelter_size; k++) {
        //if (!std::filesystem::exists(SHELTER_STORING_LOCATION + "sh[" + std::to_string(k) + "].ct")) {
        if (1) {
            // Generate random block_content of PLAINTEXT_PIR_BLOCK_DATA_SIZE bits of random | k as the block index
            Ciphertext<DCRTPoly> tmp_ct = FHE_Enc_DBElement(rng.get_z_bits(PLAINTEXT_PIR_BLOCK_DATA_SIZE), mpz_class(k));
            /* Store the ciphertexts to serialized form to a file, which resides in the RAM */
            if (Serial::SerializeToFile(SHELTER_STORING_LOCATION + "sh[" + std::to_string(k) + "].ct", tmp_ct, SerType::BINARY) == true)
            {
#if 0 /* Already performed this error checking, while doing experimentation for serveral times. Hence ommited */
            Ciphertext<DCRTPoly> deserialized_tmp_ct;
            if (Serial::DeserializeFromFile(SHELTER_STORING_LOCATION + "sh[" + std::to_string(k) + "].ct", deserialized_tmp_ct, SerType::BINARY) == false) {
                PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Failed to deserialize element FHE ciphertext from file");
            } else {
                if (*(tmp_ct) != *(deserialized_tmp_ct)) {
                    PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Deserialized element FHE ciphertext does not match original");
                }
            }
#endif
            }
            else
            {
                PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Failed to serialize element FHE ciphertext to file");
            }
        }
        
        sh[k].element_FHE_ct = import_from_file_to_mpz_class(SHELTER_STORING_LOCATION + "sh[" + std::to_string(k) + "].ct");

        /* Generate the tags and keep them in the variable, which will be used for DPF search */
        sh[k].tag = ElGamal_randomGroupElement(); // Create a random tag
        sh[k].tag_short = sh[k].tag % r; // Create a random short tag
    }

    PrintLog(LOG_LEVEL_TRACE, __FILE__, __LINE__, "Starting to test DPF-search on the shelter");

    /* Suppose we want to search for k = 2864, a random location */
    mpz_class T_sh_short = sh[DPF_SEARCH_INDEX_K].tag_short;

    // Initialize client, use 64 bits in domain as example
    initializeClient(&fClient, R_BITS, 2); // If bit length is not set properly, then incorrect answer will be returned

    // Equality FSS test
    generateTreeEq(&fClient, &k0, &k1, T_sh_short, 1);//So that the point function will evaluate as 1 at location i, and zero elsewhere

    // Initialize server
    initializeServer(&fServer, &fClient);

    mpz_class ans0, ans1, fin;
    ans0 = 0;
    ans1 = 0;

    std::vector<mpz_class> thread_sums(NUM_CPU_CORES);
    for (size_t k = 0; k < average_shelter_size; k += NUM_CPU_CORES)
    {
        for (int t = 0; t < NUM_CPU_CORES; ++t)
            thread_sums[t] = 0;

#pragma omp parallel for
        for (int j = 0; j < NUM_CPU_CORES; ++j)
        {
            if ((k + j) < average_shelter_size)
            {
                //mpz_class y = (evaluateEq(&fServer, &k0, sh[k + j].tag_short)) % mpz_class(2);// Evaluate the FSS on the short tag
                //PrintLog(LOG_LEVEL_TRACE, __FILE__, __LINE__, "For party 0, DP.Eval at: " + to_string(k + j)+ " is: " + y.get_str());
                if (evaluateEq(&fServer, &k0, sh[k + j].tag_short)) {
                    //import_from_file_to_mpz_class(SHELTER_STORING_LOCATION + "sh[" + std::to_string(k + j) + "].ct");
                    mpz_xor(thread_sums[j].get_mpz_t(), thread_sums[j].get_mpz_t(), sh[k+j].element_FHE_ct.get_mpz_t());
                    //thread_sums[j] += import_from_file_to_mpz_class(SHELTER_STORING_LOCATION + "sh[" + std::to_string(k + j) + "].ct");
                    //thread_sums[j] += sh[k + j].serialized_element_ct; // Multiply the result with the block content
                }
            }
        }
        for (int t = 0; t < NUM_CPU_CORES; ++t){
            //ans0 += thread_sums[t];
            mpz_xor(ans0.get_mpz_t(), ans0.get_mpz_t(), thread_sums[t].get_mpz_t());
        }

    }

    PrintLog(LOG_LEVEL_TRACE, __FILE__, __LINE__, "Completed DPF evaluation over average number of shelter elements");

    for (size_t k = 0; k < average_shelter_size; k += NUM_CPU_CORES)
    {
        for (int t = 0; t < NUM_CPU_CORES; ++t)
            thread_sums[t] = 0;

#pragma omp parallel for
        for (int j = 0; j < NUM_CPU_CORES; ++j)
        {
            if ((k + j) < average_shelter_size)
            {
                //mpz_class y = (evaluateEq(&fServer, &k1, sh[k + j].tag_short)) % mpz_class(2);// Evaluate the FSS on the short tag
                //PrintLog(LOG_LEVEL_TRACE, __FILE__, __LINE__, "For party 1, DP.Eval at: " + to_string(k + j)+ " is: " + y.get_str());

                if (evaluateEq(&fServer, &k1, sh[k + j].tag_short)) {
                    //import_from_file_to_mpz_class(SHELTER_STORING_LOCATION + "sh[" + std::to_string(k + j) + "].ct");
                    mpz_xor(thread_sums[j].get_mpz_t(), thread_sums[j].get_mpz_t(), sh[k+j].element_FHE_ct.get_mpz_t());
                    //thread_sums[j] += import_from_file_to_mpz_class(SHELTER_STORING_LOCATION + "sh[" + std::to_string(k + j) + "].ct");
                    //thread_sums[j] += sh[k + j].serialized_element_ct; // Multiply the result with the block content
                }
            }
        }
        for (int t = 0; t < NUM_CPU_CORES; ++t){
            //ans1 += thread_sums[t];
            mpz_xor(ans1.get_mpz_t(), ans1.get_mpz_t(), thread_sums[t].get_mpz_t());
        }
    }

    //fin = ans0 - ans1;
    mpz_xor(fin.get_mpz_t(), ans1.get_mpz_t(), ans0.get_mpz_t());

    mpz_class expected_result = import_from_file_to_mpz_class(SHELTER_STORING_LOCATION + "sh[" + std::to_string(DPF_SEARCH_INDEX_K) + "].ct");
    //mpz_class expected_result = 0;

    if (fin != expected_result) {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Cannot reform the FHE-ciphertext after DPF-search and combination");
    } else {
        PrintLog(LOG_LEVEL_INFO, __FILE__, __LINE__, "Exact same ciphertext is formed");

        export_to_file_from_mpz_class("/dev/shm/fin.ct", fin);
        // Deserialize the crypto context
        mpz_class dec_block_content, dec_block_index;
        Ciphertext<DCRTPoly> fin_ct;
        if (!Serial::DeserializeFromFile("/dev/shm/fin.ct", fin_ct, SerType::BINARY)) {
            PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Cannot read serialization from " + std::string("/dev/shm/fin.ct"));
        }  

        FHE_Dec_DBElement(fin_ct, dec_block_content, dec_block_index);

        if (dec_block_index != mpz_class(DPF_SEARCH_INDEX_K)) {
            PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Decrypted block index does not match the expected index. Expected: " + std::to_string(DPF_SEARCH_INDEX_K) + ", but got: " + dec_block_index.get_str());
        } else {
            PrintLog(LOG_LEVEL_INFO, __FILE__, __LINE__, "Decrypted block index matches the expected index");
        }

    }

    return 1;
}


int main(int argc, char *argv[])
{
    int ret = -1;

    /* Perform the basic initialization */
    InitSrv_alpha();

    /* Process as per the command line arguments */
    if (argc >= 2) {
        if (std::string("one_time_init").compare(std::string(argv[1]))==0) {
            // Perform one-time initialization for server alpha
            ret = OneTimeInit_alpha();
        } else if (std::string("per_epoch_init").compare(std::string(argv[1]))==0) {
            // Perform per-epoch initialization for server alpha
            ret = PerEpochReInit_alpha();
        } else if (std::string("clear_epoch_state").compare(std::string(argv[1]))==0) {
            // Clear the existing state of current epoch, start as if this is the first request of the epoch
            // TODO: Clear the existing state of current epoch, start as if this is the first request of the epoch
        } else if (std::string("continue").compare(std::string(argv[1]))==0) {
            // Start from last saved state
            // TODO: Start from last saved state
        } else {
            PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Unknown command line argument:"+ std::string(argv[1]));
            PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Improper command line arguments. Usage: server_alpha <one_time_init|per_epoch_init|clear_epoch_state|continue>");
        }
    } else {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Improper command line arguments. Usage: server_alpha <one_time_init|per_epoch_init|clear_epoch_state|continue>");
    }

    if (ret != 0) {
        TestSrv_alpha();
    }

    FinSrv_alpha();

    return 0;
}

static int TestClientProcessing_alpha(){
    int ret = 0;
    size_t received_sz = 0;
    size_t total_network_bytes = 0;

    PrintLog(LOG_LEVEL_TRACE, __FILE__, __LINE__, "Client timing part 1 starts here");

    /* 1.a.1.1 Select random h_{alpha1} */
    mpz_class h_alpha1 = ElGamal_randomGroupElement();
    /* 1.a.1.2 Also find its inverse */
    mpz_class h_alpha1_1;
    mpz_invert(h_alpha1_1.get_mpz_t(), h_alpha1.get_mpz_t(), p.get_mpz_t());

    /* 1.a.2.1 Select random h_{alpha2} */
    mpz_class h_alpha2 = ElGamal_randomGroupElement();
    /* 1.a.2.2 Also find its inverse */
    mpz_class h_alpha2_1;
    mpz_invert(h_alpha2_1.get_mpz_t(), h_alpha2.get_mpz_t(), p.get_mpz_t());

    /* To simulate the timing for another E() */
    ElGamal_encrypt(h_alpha1, pk_E);
    /* 2.1 Compute E(T_I.h_{\alpha 1}) to server beta */
    std::pair<mpz_class, mpz_class> E_T_I_h_alpha1 = ElGamal_mult_ct(E_T_I, ElGamal_encrypt(h_alpha1, pk_E));

    /* 2.2 Send both the componets of E(T_I.h_{\alpha 1}) */
    (void)sendAll(sock_alpha_to_beta, E_T_I_h_alpha1.first.get_str().c_str(), E_T_I_h_alpha1.first.get_str().size());
    (void)sendAll(sock_alpha_to_beta, E_T_I_h_alpha1.second.get_str().c_str(), E_T_I_h_alpha1.second.get_str().size());

    total_network_bytes += E_T_I_h_alpha1.first.get_str().size() + E_T_I_h_alpha1.second.get_str().size();
    /* To approximate the size of Eq()+2*E() */
    total_network_bytes = (total_network_bytes*3);

    PrintLog(LOG_LEVEL_INFO, __FILE__, __LINE__, "Client timing part 1 ends here");

    /* 4.a.1 Receive T_I.h_{\alpha 1}h_{\beta 0} */
    (void)recvAll(sock_alpha_to_beta, net_buf, sizeof(net_buf), &received_sz);
    if (ret != 0)
    {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Failed to receive T_I.h_{\\alpha 1}h_{\\beta 0} from Server Beta");
        return ret;
    }

    /* 5. Remove h_{\alpha 1} and determine T_I_h_beta0 */
    mpz_class T_I_h_alpha1_h_beta0 = mpz_class(std::string(net_buf, received_sz));
    mpz_class T_I_h_beta0 = (T_I_h_alpha1_h_beta0 * h_alpha1_1) % p;

    //PrintLog(LOG_LEVEL_INFO, __FILE__, __LINE__, "Determined T_I.h_{\\beta 0}: " + T_I_h_beta0.get_str());

    /* 6.1 Determine T_I.h_{\\alpha 2}.h_{\\beta 0} */
    mpz_class T_I_h_alpha2_h_beta0 = (T_I_h_beta0 * h_alpha2) % p;

    /************************************************************************************************ */
    PrintLog(LOG_LEVEL_INFO, __FILE__, __LINE__, "Client timing part 2 starts here");
    
    Ciphertext<DCRTPoly> tmp_ct = FHE_Enc_DBElement(rng.get_z_bits(PLAINTEXT_PIR_BLOCK_DATA_SIZE), mpz_class(2864));
    /* Store the ciphertexts to serialized form to a file, which resides in the RAM */
    if (Serial::SerializeToFile("/dev/shm/tmp.ct", tmp_ct, SerType::BINARY) != true){
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Error writing serialization of tmp_ct");
        return 1;
    }

    // Open the serialized file and send its contents
    std::ifstream file("/dev/shm/tmp.ct", std::ios::binary | std::ios::ate);
    if (!file)
    {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Error opening /dev/shm/tmp.ct for sending");
        return 1;
    }
    std::streamsize size = file.tellg();
    file.seekg(0, std::ios::beg);

    std::vector<char> buffer(size);
    if (!file.read(buffer.data(), size))
    {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Error reading /dev/shm/tmp.ct for sending");
        return 1;
    }

    // Send the serialized file contents
    if (sendAll(sock_alpha_to_beta, buffer.data(), size) != 0)
    {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Error sending serialized ciphertext to server beta");
        return 1;
    }

    size_t FHE_ct_sz = size;
    total_network_bytes += size;

    /* Receiving a dummy block */
    (void)recvAll(sock_alpha_to_beta, net_buf, sizeof(net_buf), &received_sz);
    if (ret != 0)
    {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Failed to receive T_I.h_{\\alpha 1}h_{\\beta 0} from Server Beta");
        return ret;
    }

    total_network_bytes += received_sz;

    /* Simulating the time requirement for removal of mask by adding two random numbers */
    mpz_class temp = rng.get_z_bits(PLAINTEXT_PIR_BLOCK_DATA_SIZE+P_BITS) + rng.get_z_bits(PLAINTEXT_PIR_BLOCK_DATA_SIZE+P_BITS);

    PrintLog(LOG_LEVEL_INFO, __FILE__, __LINE__, "Client timing part 2 ends here, size of FHE-ciphertext is: " + std::to_string(FHE_ct_sz) + ", while total network transfer size is: " + std::to_string(total_network_bytes));
    /************************************************************************************* */

    return 0;
}
