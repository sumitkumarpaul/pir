#include <assert.h>
#include <stdio.h>
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
static int sock_alpha_client_srv = -1, sock_alpha_client_con = -1;
static char net_buf[NET_BUF_SZ] = {0};

// The value of N will determine the bitlength during the client initialization
#define NUM_TAG_BITS 3072 // 16 bits can represent up to 65536, which is more than enough for N=50000
#define B 512 // Block size in bits, can be adjusted as needed
// And number of bits determine the evalution time drastically
static shelter_element sh[sqrt_N]; // Database to store values, each entry is a tuple.
static uint64_t K; // Current number of entries in the shelter, or the number of processed requests
static KukuTable *HTable = nullptr;

std::pair<mpz_class, mpz_class> E_T_I;
static mpz_class a;

#define CUCKOO_HASH_TABLE_REHASH_TRY_COUNT 1
#define NUM_CPU_CORES 16

#define DPF_SEARCH_INDEX_K 1
//#define SHELTER_STORING_LOCATION std::string("./")
#define SHELTER_STORING_LOCATION std::string("/mnt/sumit/dummy_shelter/")

#define ONE_TIME_MATERIALS_LOCATION_ALPHA std::string("/mnt/sumit/PIR_ALPHA/ONE_TIME_MATERIALS/")
#define PER_EPOCH_MATERIALS_LOCATION_ALPHA std::string("/mnt/sumit/PIR_ALPHA/PER_EPOCH_MATERIALS/")
#define DATABASE_LOCATION_ALPHA std::string("/mnt/sumit/PIR_ALPHA/")
std::string L_filename = PER_EPOCH_MATERIALS_LOCATION_ALPHA+"L_alpha.bin";
std::string DK_filename = PER_EPOCH_MATERIALS_LOCATION_ALPHA+"DK_alpha.bin";//TODO: The key and data are seperated, unlike the description of the paper
std::string sdb_filename = PER_EPOCH_MATERIALS_LOCATION_ALPHA+"ShuffledDB_alpha.bin";
std::string HTable_filename = PER_EPOCH_MATERIALS_LOCATION_ALPHA+"H_alpha.bin";


// Function declarations
static int InitSrv_alpha();
static int OneTimeInit_alpha();
static int FinSrv_alpha();
static int PerEpochOperations_alpha();
static int ProcessClientRequest_alpha();
static int SelShuffDBSearchTag_alpha();
static int ShelterTagDetermination_alpha();
static int ObliviouslySearchShelter_alpha();

static void TestSrv_alpha();
static void TestPKEOperations_alpha();
static int TestHTableSerDser_alpha();
static void TestSelShuffDBSearchTag_alpha();
static int TestShelterDPFSearch_alpha();
static int TestClientProcessing_alpha();

// Function definitions
static int InitSrv_alpha(){
    int ret = -1;
    // Initialize random number generation
    std::random_device rd;
    unsigned long seed = (static_cast<unsigned long>(rd()) << 1) ^ rd();
    rng.seed(seed); // seed() seeds the gmp_randclass    
    
    // Server_alpha only connects to other servers, it does not listen to other servers
    InitConnectingSocket(SERVER_BETA_IP, BETA_LISTENING_TO_ALPHA_PORT, &sock_alpha_to_beta);

    PrintLog(LOG_LEVEL_TRACE, __FILE__, __LINE__, "Established connection with Server Beta");

    InitConnectingSocket(SERVER_GAMMA_IP, GAMMA_LISTENING_TO_ALPHA_PORT, &sock_alpha_to_gamma);

    PrintLog(LOG_LEVEL_TRACE, __FILE__, __LINE__, "Established connection with Server Gamma");

    PrintLog(LOG_LEVEL_INFO, __FILE__, __LINE__, "Server Alpha initialization complete");

    ret = 0;

exit:
    if (ret != 0){
        FinSrv_alpha();
    }

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

    PrintLog(LOG_LEVEL_TRACE, __FILE__, __LINE__, "Received all the one-time initialized parameters from Server Beta and exported all of them into file");
    
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
    if (sock_alpha_client_srv != -1) {
        close(sock_alpha_client_srv);
        sock_alpha_client_srv = -1;
    }
    if (sock_alpha_client_con != -1) {
        close(sock_alpha_client_con);
        sock_alpha_client_con = -1;
    }

    PrintLog(LOG_LEVEL_SPECIAL, __FILE__, __LINE__, "Finalized Server Alpha");

    return ret;
}

/* TODO: This function, currently does not matches with the diagram of the paper. It has to be updated or the diagram has to be modified. */
static int PerEpochOperations_alpha(){
    int ret = 0;
    size_t received_sz = 0;
    shuffled_db_entry sdb_entry;
    uint64_t M = (N + sqrt_N);
    item_type Kuku_key;    
    QueryResult res;
    std::fstream L;
    std::fstream DK;
    std::fstream sdb;
    std::ofstream exportedHFile;    

    PrintLog(LOG_LEVEL_INFO, __FILE__, __LINE__, "Server Alpha: Starting PerEpochOperations sequence");

    /* First of all retrieve all the one-time initialized materials from the saved location */
    p = import_from_file_to_mpz_class(ONE_TIME_MATERIALS_LOCATION_ALPHA + "p.bin");
    q = import_from_file_to_mpz_class(ONE_TIME_MATERIALS_LOCATION_ALPHA + "q.bin");
    g = import_from_file_to_mpz_class(ONE_TIME_MATERIALS_LOCATION_ALPHA + "g.bin");
    g_q = import_from_file_to_mpz_class(ONE_TIME_MATERIALS_LOCATION_ALPHA + "g_q.bin");
    r = import_from_file_to_mpz_class(ONE_TIME_MATERIALS_LOCATION_ALPHA + "r.bin");
    pk_E = import_from_file_to_mpz_class(ONE_TIME_MATERIALS_LOCATION_ALPHA + "pk_E.bin");
    pk_E_q = import_from_file_to_mpz_class(ONE_TIME_MATERIALS_LOCATION_ALPHA + "pk_E_q.bin");
    Serial::DeserializeFromFile(ONE_TIME_MATERIALS_LOCATION_ALPHA + "FHEcryptoContext.bin", FHEcryptoContext, SerType::BINARY);
    Serial::DeserializeFromFile(ONE_TIME_MATERIALS_LOCATION_ALPHA + "pk_F.bin", pk_F, SerType::BINARY);
    Serial::DeserializeFromFile(ONE_TIME_MATERIALS_LOCATION_ALPHA + "vectorOnesforElement_ct.bin", vectorOnesforElement_ct, SerType::BINARY);
    Serial::DeserializeFromFile(ONE_TIME_MATERIALS_LOCATION_ALPHA + "vectorOnesforTag_ct.bin", vectorOnesforTag_ct, SerType::BINARY);

    PrintLog(LOG_LEVEL_TRACE, __FILE__, __LINE__, "Server Alpha: Loaded one-time initialization materials");

    /* Wait for receiving the ready message from server-beta */
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

    /* Receive completed message from server-beta */
    (void)recvAll(sock_alpha_to_beta, net_buf, sizeof(net_buf), &received_sz);
    
    if (ret != 0)
    {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Failed to receive COMPLETED_REINIT_FOR_EPOCH message from Server Beta");
        return -1;
    }

    if (std::string(net_buf, received_sz) != completed_reinit_for_epoch_message) {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Did not receive expected COMPLETED_REINIT_FOR_EPOCH message from Server Beta");
        return -1;
    } else {
        PrintLog(LOG_LEVEL_INFO, __FILE__, __LINE__, "Server Alpha: Completed re-initialization for new epoch, now ready to process client-requests..!!");
    }

    /* TODO: Unlike the paper. The L file only contains the secret shared data part */
    /* The file should already contain all the secret shares */
    L.open(L_filename, std::ios::in | std::ios::binary);
    if (!L) {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Failed to open L file at location: " + L_filename);
        return -1;
    }

    /* TODO: Unlike the paper. The K file contains the Kuku keys, generated from the Tags */
    /* The file should already contain all the keys */
    DK.open(DK_filename, std::ios::in | std::ios::binary);
    if (!DK) {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Failed to open K file at location: " + DK_filename);
        ret = -1;
        goto exit;
    }

    /* 13.a.1 This file stores the pair: (Kuku key, the share of the element) */
    sdb.open(sdb_filename, std::ios::in | std::ios::out | std::ios::binary | std::ios::trunc);
    if (!sdb) {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Failed to open SDB file at location: " + sdb_filename);
        ret = -1;
        goto exit;
    }

    // 12.a.1 Allocate a new Cuckoo hash table
    // Keeping the number of entries, larger than the number of elements to place. The reason is, it will reduce the number of probe during placement and make the per epoch operations faster
    HTable = new KukuTable(CUCKOO_TABLE_SIZE, CUCKOO_STASH_SIZE, CUCKOO_LOC_FUNC_COUNT, CUCKOO_LOC_FUNC_SEED, CUCKOO_MAX_PROBE, CUCKOO_EMPTY_ITEM);
    if (HTable == nullptr) {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Failed to allocate memory for Cuckoo hash table");
        ret = -1;
        goto exit;
    }

    /* 12.a.2 Prepare the entire Kuku hash table, based on all the keys */
    for (uint64_t i = 0; i < M; i++){
        /* Read the next Kuku key from the DK file */
        DK.read(reinterpret_cast<char*>(Kuku_key.data()), sizeof(item_type));

        if (!HTable->insert(Kuku_key))
        {
            PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Insertion failed. Before failure, successfully inserted: " + to_string(i) + " out of " + to_string(M) + " items");
            PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "The size of the stash during failure: " + to_string(HTable->stash().size()) + " and max-probe count is: " + to_string(HTable->max_probe()));

            /* Delete the already built table */
            delete HTable;
            HTable = nullptr;
            ret = -1;
            goto exit;
        }

        if (((i+1) % 100000000) == 0){
            PrintLog(LOG_LEVEL_TRACE, __FILE__, __LINE__, "Inserted " + to_string(i+1) + " items into the cuckoo hash HTable. Current stash size: " + to_string(HTable->stash().size()) + " and total probe count is: " + to_string (HTable->total_probe_count_));
        }
    }

    PrintLog(LOG_LEVEL_INFO, __FILE__, __LINE__, "Cuckoo hash table creation complete. Current stash size: " + to_string(HTable->stash().size()) + " and total probe count is: " + to_string (HTable->total_probe_count_));

    //13.a.2 Reset the read pointers to the beginning
    DK.seekg(0, std::ios::beg);
    L.seekg(0, std::ios::beg);
    sdb.seekp(0, std::ios::beg);

    //13.a.3 Now place all the elements from the temporary list(L) to the shuffled database(SDB) according to the Kuku hash table
    for (uint64_t i = 0; i < M; i++){
        /* Read the next key */
        DK.read(reinterpret_cast<char*>(Kuku_key.data()), sizeof(item_type));
        /* Read the next secret share */
        L.read(net_buf, NUM_BYTES_PER_SDB_ELEMENT);

        /* 13.a.4: Prepare them to a tuple of shuffled database */
        memcpy(sdb_entry.cuckoo_key.data(), Kuku_key.data(), sizeof(item_type));
        memcpy(sdb_entry.element, net_buf, NUM_BYTES_PER_SDB_ELEMENT);

        /* 13.a.5: Query and find the location */
        res = HTable->query(Kuku_key);
        if (!res)
        {
            PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Query failed for the item number: " + to_string(i));
            ret = -1;
            goto exit;
        }
        else {
            /* 13.a.5 Insert at the location of the shuffled database, determined by the query result */
            insert_sdb_entry(sdb, res.location(), sdb_entry);
            PrintLog(LOG_LEVEL_TRACE, __FILE__, __LINE__, "Inserted item: " + to_string(i) + " of L to SDB at location: " + std::to_string(res.location()));
        }

        if (((i+1) % 100000000) == 0){
            PrintLog(LOG_LEVEL_TRACE, __FILE__, __LINE__, "Inserted " + to_string(i+1) + " items into the shuffled database");
        }
    }

    PrintLog(LOG_LEVEL_INFO, __FILE__, __LINE__, "Shuffled database creation complete");
    PrintLog(LOG_LEVEL_INFO, __FILE__, __LINE__, "TODO: Check whether ith item of DK and L are really placed in proper location of SDB");

    // 14.a. Nothing is required to be done for clearing the shelter content

    //Store the cuckoo table in the disk
    exportedHFile.open(HTable_filename, std::ios::binary);
    HTable->serialize(exportedHFile);
    exportedHFile.close();

exit:
    L.close();
    DK.close();
    sdb.close();

    return ret;
}

static int ShelterTagDetermination_alpha(){
    int ret = 0;
    size_t received_sz = 0;
    int ret_recv = 0;
    std::pair<mpz_class, mpz_class> E_g_pow_Rho_pow_I__mul__h_C;
    mpz_class Rho_pow_I__mul__h_C;
    mpz_class g_pow_Rho_pow_I__mul__h_C;
    mpz_class h_alpha0, h_alpha0_1;
    std::pair<mpz_class, mpz_class> E_g_pow_Rho_pow_I__mul__h_C_h_alpha0;
    std::pair<mpz_class, mpz_class> E_g_pow_Rho_pow_I__mul_h_alpha0;
    std::pair<mpz_class, mpz_class> E_g_pow_Rho_pow_I;
    std::pair<mpz_class, mpz_class> E_a;
    std::pair<mpz_class, mpz_class> E_g_pow_Rho_pow_I__mul_a;    

    /* Step 4.3.1 of the sequence diagram */
    // Receive the first component of E_g_pow_Rho_pow_I__mul__h_C
    ret_recv = recvAll(sock_alpha_to_beta, net_buf, sizeof(net_buf), &received_sz);
    if (ret_recv != 0)
    {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Failed to receive E_g_pow_Rho_pow_I__mul__h_C.first from the Server Beta");
        return -1;
    }
    PrintLog(LOG_LEVEL_TRACE, __FILE__, __LINE__, "Server Alpha starts client request processing from this point");

    E_g_pow_Rho_pow_I__mul__h_C.first = mpz_class(std::string(net_buf, received_sz));

    // Step 4.3.2 Receive the second component of E_g_pow_Rho_pow_I__mul__h_C
    ret_recv = recvAll(sock_alpha_to_beta, net_buf, sizeof(net_buf), &received_sz);
    if (ret_recv != 0)
    {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Failed to receive E_g_pow_Rho_pow_I__mul__h_C.second from the Server Beta");
        return -1;
    }

    E_g_pow_Rho_pow_I__mul__h_C.second = mpz_class(std::string(net_buf, received_sz));

    // Step 5.1 Choose h_{\alpha 0} and its inverse
    h_alpha0 = rng.get_z_range(q - 1) + 1; // i.e., within ZZ_q*
    mpz_invert(h_alpha0_1.get_mpz_t(), h_alpha0.get_mpz_t(), q.get_mpz_t());

    // Step 5.2 Semi-homomorphically raises that to the received ciphertext (under ElGamal encryption in GG)
    E_g_pow_Rho_pow_I__mul__h_C_h_alpha0 = ElGamal_exp_ct(E_g_pow_Rho_pow_I__mul__h_C, h_alpha0, pk_E);

    // Step 6.1 Send the first part to the client
    (void)sendAll(sock_alpha_client_con, E_g_pow_Rho_pow_I__mul__h_C_h_alpha0.first.get_str().c_str(), E_g_pow_Rho_pow_I__mul__h_C_h_alpha0.first.get_str().size());

    // Step 6.2 Send the second part to the client
    (void)sendAll(sock_alpha_client_con, E_g_pow_Rho_pow_I__mul__h_C_h_alpha0.second.get_str().c_str(), E_g_pow_Rho_pow_I__mul__h_C_h_alpha0.second.get_str().size());

    // Step 7.3.1 Receive the first component of E_g_pow_Rho_pow_I__mul_h_alpha0
    ret_recv = recvAll(sock_alpha_client_con, net_buf, sizeof(net_buf), &received_sz);
    if (ret_recv != 0)
    {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Failed to receive E_g_pow_Rho_pow_I__mul_h_alpha0.first from the Client");
        return -1;
    }

    E_g_pow_Rho_pow_I__mul_h_alpha0.first = mpz_class(std::string(net_buf, received_sz));

    // Step 7.3.2 Receive the second component of E_g_pow_Rho_pow_I__mul_h_alpha0
    ret_recv = recvAll(sock_alpha_client_con, net_buf, sizeof(net_buf), &received_sz);
    if (ret_recv != 0)
    {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Failed to receive E_g_pow_Rho_pow_I__mul_h_alpha0.second from the Client");
        return -1;
    }

    E_g_pow_Rho_pow_I__mul_h_alpha0.second = mpz_class(std::string(net_buf, received_sz));

    // Step 8. Semi-homomorphically raises that to h_alpha0_1 under ElGamal encryption in GG to remove h_alpha0
    E_g_pow_Rho_pow_I = ElGamal_exp_ct(E_g_pow_Rho_pow_I__mul_h_alpha0, h_alpha0_1, pk_E);

    // Step 9.1. Compute the ciphertext of a
    E_a = ElGamal_encrypt(a, pk_E);

    // Step 9.2. Multiply homomorphically
    E_g_pow_Rho_pow_I__mul_a = ElGamal_mult_ct(E_g_pow_Rho_pow_I, E_a);

    // Step 9.3.1 Send the first part of the ciphertext to the Server Gamma
    (void)sendAll(sock_alpha_to_gamma, E_g_pow_Rho_pow_I__mul_a.first.get_str().c_str(), E_g_pow_Rho_pow_I__mul_a.first.get_str().size());

    // Step 9.4.1 Send the second part of the ciphertext to the Server Gamma
    (void)sendAll(sock_alpha_to_gamma, E_g_pow_Rho_pow_I__mul_a.second.get_str().c_str(), E_g_pow_Rho_pow_I__mul_a.second.get_str().size());

    /* This is only for experimentation purpose */
    //PrintLog(LOG_LEVEL_TRACE, __FILE__, __LINE__, "Chosen a is: " + a.get_str());

    return ret;
}

static int ObliviouslySearchShelter_alpha() {
    // Set up variables
    Fss fServer;
    ServerKeyEq K_alpha;
    int ret = 0;
    size_t received_sz = 0;
    int ret_recv;
    size_t dserializedFssSize;
    Ciphertext<DCRTPoly> fnd_ct_element, fnd_alpha_ct_element, fnd_gamma_ct_element;
    Ciphertext<DCRTPoly> fnd_ct_tag, fnd_alpha_ct_tag, fnd_gamma_ct_tag;
    mpz_class d_ct_alpha = 0;
    mpz_class d_ct_gamma, d_ct;
    Ciphertext<DCRTPoly> d_ct_FHE;
    std::vector<bool> thread_fnd(NUM_CPU_CORES, false);
    bool fnd_alpha = false;
    std::vector<mpz_class> thread_sums(NUM_CPU_CORES);

    // First, receive FSS parameters from the server Beta
    ret = recvAll(sock_alpha_to_beta, net_buf, sizeof(net_buf), &received_sz);
    if (ret != 0)
    {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Failed to receive FSS parameters from Server Beta");
        return -1;
    }

    dserializedFssSize = deserializeFssAndServerKeyEq(net_buf, received_sz, fServer, K_alpha);
    PrintLog(LOG_LEVEL_TRACE, __FILE__, __LINE__, "Total size of the dserialized data: " + std::to_string(dserializedFssSize));

    /* TODO Dummy print to verify the gmp_class is actually got transferred */
    PrintLog(LOG_LEVEL_TRACE, __FILE__, __LINE__, "Received values from server Beta is: ");
    PrintLog(LOG_LEVEL_TRACE, __FILE__, __LINE__, "fServer.numBits: " + std::to_string(fServer.numBits));
    PrintLog(LOG_LEVEL_TRACE, __FILE__, __LINE__, "fServer.prime: " + (fServer.prime).get_str());
    PrintLog(LOG_LEVEL_TRACE, __FILE__, __LINE__, "fServer.numParties: " + std::to_string(fServer.numParties));
    PrintLog(LOG_LEVEL_TRACE, __FILE__, __LINE__, "fServer.numKeys: " + std::to_string(fServer.numParties));
    PrintLog(LOG_LEVEL_TRACE, __FILE__, __LINE__, "K_alpha.w: " + (K_alpha.w).get_str());

    /* For the verification purpose, set a particular location with special tag printed from server beta, to make the DPF search successful */
#if 1
    mpz_class special_tag;

    printf("Enter the value of the set search tag (base 10): ");
    mpz_inp_str(special_tag.get_mpz_t(), stdin, 10);

    sh[0].tag_short = special_tag;
#endif


    PrintLog(LOG_LEVEL_TRACE, __FILE__, __LINE__, "Starting to test DPF-search on the shelter");

    for (size_t k = 0; k < K; k += NUM_CPU_CORES)
    {
        for (int t = 0; t < NUM_CPU_CORES; ++t)
            thread_sums[t] = 0;

#pragma omp parallel for
        for (int j = 0; j < NUM_CPU_CORES; ++j)
        {
            if ((k + j) < K)
            {
                if (evaluateEq(&fServer, &K_alpha, sh[k + j].tag_short)) {
                    mpz_xor(thread_sums[j].get_mpz_t(), thread_sums[j].get_mpz_t(), sh[k+j].element_FHE_ct.get_mpz_t());

                    /* Same as XORing */
                    thread_fnd[j] = !thread_fnd[j];
                    printf("1 ");
                }
                else{
                    printf("0 ");
                }
            }
        }
        for (int t = 0; t < NUM_CPU_CORES; ++t)
        {
            mpz_xor(d_ct_alpha.get_mpz_t(), d_ct_alpha.get_mpz_t(), thread_sums[t].get_mpz_t());
        }
    }
    for (int t = 0; t < NUM_CPU_CORES; ++t)
    {
        fnd_alpha ^= thread_fnd[t];
    }

    PrintLog(LOG_LEVEL_TRACE, __FILE__, __LINE__, "Completed DPF evaluation. Value of fnd_alpha: " + std::to_string(fnd_alpha));

    // 7.2.2 Receive fnd_gamma_ct_element
    ret_recv = recvAll(sock_alpha_to_gamma, net_buf, sizeof(net_buf), &received_sz);
    if (ret_recv != 0)
    {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Failed to receive fnd_gamma_ct_element from Server Gamma");
        return -1;
    }
    Serial::DeserializeFromString(fnd_gamma_ct_element, std::string(net_buf, received_sz));

    // 7.3.2 Receive fnd_gamma_ct_tag
    ret_recv = recvAll(sock_alpha_to_gamma, net_buf, sizeof(net_buf), &received_sz);
    if (ret_recv != 0)
    {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Failed to receive fnd_gamma_ct_tag from Server Gamma");
        return -1;
    }
    Serial::DeserializeFromString(fnd_gamma_ct_tag, std::string(net_buf, received_sz));

    // 7.4.2 Receive d_ct_gamma
    ret_recv = recvAll(sock_alpha_to_gamma, net_buf, sizeof(net_buf), &received_sz);
    if (ret_recv != 0)
    {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Failed to receive d_ct_gamma from Server Gamma");
        return -1;
    }
    d_ct_gamma = mpz_class(std::string(net_buf, received_sz));

    PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "HERE");

    /* Step 8.1 Determine fnd_alpha_ct */
    if (fnd_alpha == true){
        FHE_EncOfOnes(fnd_alpha_ct_element, fnd_alpha_ct_tag);
    }else{
        FHE_EncOfZeros(fnd_alpha_ct_element, fnd_alpha_ct_tag);
    }
    PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "HERE");

    /* Step 8.2 Homomorphically compute fnd_ct = fnd_alpha_ct XOR fnd_alpha_ct = (fnd_alpha_ct + fnd_alpha_ct) - 2*(fnd_alpha_ct*fnd_alpha_ct) */
    /* Only for debugging purpose receive private key from server beta */
    ret = recvAll(sock_alpha_to_beta, net_buf, sizeof(net_buf), &received_sz);
    if (ret != 0)
    {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Failed to receive sk_F from Server Beta");
        return -1;
    }
    Serial::DeserializeFromString(sk_F, std::string(net_buf, received_sz));

    PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "HERE");

    fnd_ct_element = FHE_bitwise_XOR(fnd_alpha_ct_element, fnd_gamma_ct_element);
    fnd_ct_tag = FHE_bitwise_XOR(fnd_alpha_ct_tag, fnd_gamma_ct_tag);

    /* Step 8.2.1 Compute d_ct in mpz_class */
    mpz_xor(d_ct.get_mpz_t(), d_ct_gamma.get_mpz_t(), d_ct_alpha.get_mpz_t());
    PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "HERE");

    /* Setp 8.3 Convert from mpz_class to FHE ciphertext */
    export_to_file_from_mpz_class("/dev/shm/d.ct", d_ct);

    if (!Serial::DeserializeFromFile("/dev/shm/d.ct", d_ct_FHE, SerType::BINARY)) {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Cannot convert to the d_ct to Ciphertext<DCRTPoly>");
    }
    
    /* This one only for experimentation purpose */
    mpz_class dec_block_content, dec_block_index, dec_content_and_index;
    mpz_class dec_fnd_ct_element, dec_fnd_ct_tag;

    PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "HERE");
    
    FHE_Dec_SDBElement(d_ct_FHE, dec_content_and_index);
    PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "HERE");

    dec_block_content = (dec_content_and_index >> log_N);
    dec_block_index = (dec_content_and_index & ((1U << log_N) - 1U));

    PrintLog(LOG_LEVEL_TRACE, __FILE__, __LINE__, "The value of decrypted data: " + dec_block_content.get_str()+ " index: " + dec_block_index.get_str());

    FHE_Dec_SDBElement(fnd_ct_element, dec_fnd_ct_element);
    PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "HERE");

    PrintLog(LOG_LEVEL_TRACE, __FILE__, __LINE__, "The value of decrypted fnd_ct_element: " + dec_fnd_ct_element.get_str());

    FHE_Dec_Tag(fnd_ct_tag, dec_fnd_ct_tag);
    PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "HERE");

    PrintLog(LOG_LEVEL_TRACE, __FILE__, __LINE__, "The value of decrypted dec_fnd_ct_tag: " + dec_fnd_ct_tag.get_str());    

    return 0;
}

//TODO static int ResumeRequestProcessing_alpha(){
static int ProcessClientRequest_alpha(){
    int ret = -1;
    //struct sockaddr_in address;
    //int opt = 1;
    //int addrlen = sizeof(address);
    //int accepted_socket = -1;
    shuffled_db_entry sdb_entry;
    uint64_t M = (N + sqrt_N);
    item_type Kuku_key;    
    QueryResult res;
    std::fstream L;
    std::fstream DK;
    std::fstream sdb;
    std::ifstream importedHFile;

    PrintLog(LOG_LEVEL_INFO, __FILE__, __LINE__, "Server Alpha: Starting Request processing sequence");

    /* First of all retrieve all the one-time initialized materials from the saved location */
    p = import_from_file_to_mpz_class(ONE_TIME_MATERIALS_LOCATION_ALPHA + "p.bin");
    q = import_from_file_to_mpz_class(ONE_TIME_MATERIALS_LOCATION_ALPHA + "q.bin");
    g = import_from_file_to_mpz_class(ONE_TIME_MATERIALS_LOCATION_ALPHA + "g.bin");
    g_q = import_from_file_to_mpz_class(ONE_TIME_MATERIALS_LOCATION_ALPHA + "g_q.bin");
    r = import_from_file_to_mpz_class(ONE_TIME_MATERIALS_LOCATION_ALPHA + "r.bin");
    pk_E = import_from_file_to_mpz_class(ONE_TIME_MATERIALS_LOCATION_ALPHA + "pk_E.bin");
    pk_E_q = import_from_file_to_mpz_class(ONE_TIME_MATERIALS_LOCATION_ALPHA + "pk_E_q.bin");
    Serial::DeserializeFromFile(ONE_TIME_MATERIALS_LOCATION_ALPHA + "FHEcryptoContext.bin", FHEcryptoContext, SerType::BINARY);
    Serial::DeserializeFromFile(ONE_TIME_MATERIALS_LOCATION_ALPHA + "pk_F.bin", pk_F, SerType::BINARY);
    Serial::DeserializeFromFile(ONE_TIME_MATERIALS_LOCATION_ALPHA + "vectorOnesforElement_ct.bin", vectorOnesforElement_ct, SerType::BINARY);
    Serial::DeserializeFromFile(ONE_TIME_MATERIALS_LOCATION_ALPHA + "vectorOnesforTag_ct.bin", vectorOnesforTag_ct, SerType::BINARY);

    PrintLog(LOG_LEVEL_TRACE, __FILE__, __LINE__, "Server Alpha: Loaded one-time initialization materials");

    //Always initialize them
    K = 0;
    a = 1;

    importedHFile.open(HTable_filename, std::ios::binary);
    if (!importedHFile) {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Failed to open H file at location: " + HTable_filename);
        ret = -1;
        goto exit;
    }

    HTable = KukuTable::deserialize(importedHFile).release();
    importedHFile.close();
    if (HTable == nullptr) {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Cannot import the hash table");
        ret = -1;
        goto exit;
    }
    PrintLog(LOG_LEVEL_TRACE, __FILE__, __LINE__, "Server Alpha: Loaded hash table into the RAM");

    /* For debugging loading the shelter into the RAM, in real situation the shelter will always remain in the RAM */
    for(size_t k = 0; k < sqrt_N; k++) {
        sh[k].element_FHE_ct = import_from_file_to_mpz_class(SHELTER_STORING_LOCATION + "sh[" + std::to_string(k) + "].ct");

        /* Generate the tags and keep them in the variable, which will be used for DPF search */
        sh[k].tag = import_from_file_to_mpz_class(SHELTER_STORING_LOCATION + "sh[" + std::to_string(k) + "].T_hat");
        sh[k].tag_short = import_from_file_to_mpz_class(SHELTER_STORING_LOCATION + "sh[" + std::to_string(k) + "].t_hat"); // Create a random short tag
    }

    while (K < sqrt_N){
        PrintLog(LOG_LEVEL_INFO, __FILE__, __LINE__, "Waiting for processing the PIR request number: "+ std::to_string(K+1) +" from the client..!!");

        ret = InitAcceptingSocket(ALPHA_LISTENING_TO_CLIENT_PORT, &sock_alpha_client_srv, &sock_alpha_client_con);

        if (ret != 0)
        {
            PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Cannot establish communication with the client!!");
            ret = -1;
            goto exit;
        }

        ret = ShelterTagDetermination_alpha();
        if (ret != 0)
        {
            PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Problem while determining the shelter tag..!!");
            ret = -1;
            goto exit;
        }

        /* For the first request, the shelter is not required to be searched */
        if (K > 0){
            ret = ObliviouslySearchShelter_alpha();
            if (ret != 0)
            {
                PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Problem during the shelter search operation..!!");
                ret = -1;
                goto exit;
            }
        }

        /* Other sequences */
        a = rng.get_z_range(p); /* TODO, this will be part of shelter update */

        /* Close the connection with existing client */
        close(sock_alpha_client_srv);
        close(sock_alpha_client_con);  
        K++;
    }

    PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Current epoch is completed. Please re-perform the per-epoch initialization");


exit:
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
    TestShelterDPFSearch_alpha();
    //TestClientProcessing_alpha();
    //TestHTableSerDser_alpha();

    return;
}

ostream &operator<<(ostream &stream, item_type item)
{
    stream << item[1] << " " << item[0];
    return stream;
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

    /* First of all retrieve all the one-time initialized materials from the saved location */
    p = import_from_file_to_mpz_class(ONE_TIME_MATERIALS_LOCATION_ALPHA + "p.bin");
    q = import_from_file_to_mpz_class(ONE_TIME_MATERIALS_LOCATION_ALPHA + "q.bin");
    g = import_from_file_to_mpz_class(ONE_TIME_MATERIALS_LOCATION_ALPHA + "g.bin");
    g_q = import_from_file_to_mpz_class(ONE_TIME_MATERIALS_LOCATION_ALPHA + "g_q.bin");
    r = import_from_file_to_mpz_class(ONE_TIME_MATERIALS_LOCATION_ALPHA + "r.bin");
    pk_E = import_from_file_to_mpz_class(ONE_TIME_MATERIALS_LOCATION_ALPHA + "pk_E.bin");
    pk_E_q = import_from_file_to_mpz_class(ONE_TIME_MATERIALS_LOCATION_ALPHA + "pk_E_q.bin");
    Serial::DeserializeFromFile(ONE_TIME_MATERIALS_LOCATION_ALPHA + "FHEcryptoContext.bin", FHEcryptoContext, SerType::BINARY);
    Serial::DeserializeFromFile(ONE_TIME_MATERIALS_LOCATION_ALPHA + "pk_F.bin", pk_F, SerType::BINARY);
    Serial::DeserializeFromFile(ONE_TIME_MATERIALS_LOCATION_ALPHA + "vectorOnesforElement_ct.bin", vectorOnesforElement_ct, SerType::BINARY);
    Serial::DeserializeFromFile(ONE_TIME_MATERIALS_LOCATION_ALPHA + "vectorOnesforTag_ct.bin", vectorOnesforTag_ct, SerType::BINARY);

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
            Ciphertext<DCRTPoly> tmp_ct = FHE_Enc_SDBElement((rng.get_z_bits(PLAINTEXT_PIR_BLOCK_DATA_SIZE) << log_N) | mpz_class(k));
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
        mpz_class dec_block_content, dec_block_index, dec_content_and_index;
        Ciphertext<DCRTPoly> fin_ct;
        if (!Serial::DeserializeFromFile("/dev/shm/fin.ct", fin_ct, SerType::BINARY)) {
            PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Cannot read serialization from " + std::string("/dev/shm/fin.ct"));
        }
        PrintLog(LOG_LEVEL_TRACE, __FILE__, __LINE__, "Here");
        
        FHE_Dec_SDBElement(fin_ct, dec_content_and_index);
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Here");
        
        dec_block_content = (dec_content_and_index >> log_N);
        dec_block_index = (dec_content_and_index & ((1U << log_N) - 1U)); 

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
        } else if (std::string("per_epoch_operations").compare(std::string(argv[1]))==0) {
            // Perform per-epoch initialization for server alpha
            ret = PerEpochOperations_alpha();
        } else if (std::string("clear_epoch_state").compare(std::string(argv[1]))==0) {
            // Clear the existing state of current epoch, start as if this is the first request of the epoch
            // Delete shelter content and set K = 0
        } else if (std::string("process_request").compare(std::string(argv[1]))==0) {
            // Start from last saved state
            ret = ProcessClientRequest_alpha();
        } else if (std::string("test").compare(std::string(argv[1]))==0) {
            TestSrv_alpha();
        } else {
            PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Unknown command line argument:"+ std::string(argv[1]));
            PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Improper command line arguments. Usage: server_alpha <one_time_init|per_epoch_operations|clear_epoch_state|process_request>");
        }
    } else {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Improper command line arguments. Usage: server_alpha <one_time_init|per_epoch_operations|clear_epoch_state|process_request>");
    }

    if (ret == 0) {
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
    
    Ciphertext<DCRTPoly> tmp_ct = FHE_Enc_SDBElement((rng.get_z_bits(PLAINTEXT_PIR_BLOCK_DATA_SIZE) << log_N) | mpz_class(2864));
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

static int TestHTableSerDser_alpha(){
    uint64_t M = (N + sqrt_N);
    QueryResult res;
    std::fstream DK;
    item_type Kuku_key;

    DK.open(DK_filename, std::ios::in | std::ios::binary);
    if (!DK) {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Failed to open K file at location: " + DK_filename);
        goto exit;
    }    

    /* Print all the mappings and verify it with running of the program with last executing with per_epoch_operations */
    for (uint64_t i = 0; i < M; i++){
        DK.read(reinterpret_cast<char*>(Kuku_key.data()), sizeof(item_type));
        res = HTable->query(Kuku_key);
        if (!res)
        {
            PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Query failed for the item number: " + to_string(i));
        }
        else {
            PrintLog(LOG_LEVEL_TRACE, __FILE__, __LINE__, "Item: " + to_string(i) + " of L is mapped to SDB at location: " + std::to_string(res.location()));
        }
    }

exit:
    DK.close();

    return 0;
}