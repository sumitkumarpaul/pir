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


static int sock_gamma_to_beta = -1, sock_gamma_to_alpha = -1, sock_gamma_to_alpha_con = -1;
static int sock_gamma_client_srv = -1, sock_gamma_client_con = -1;
static char net_buf[NET_BUF_SZ] = {0};

// The value of N will determine the bitlength during the client initialization
#define NUM_TAG_BITS 3072 // 16 bits can represent up to 65536, which is more than enough for N=50000
#define B 512 // Block size in bits, can be adjusted as needed
// And number of bits determine the evalution time drastically
static mpz_class sh[sqrt_N][2]; // Database to store values, each entry is a pair, {Tag, Block-content}.
mpz_class T_star;
static uint64_t K; // Current number of entries in the shelter, or the number of processed requests
static KukuTable *HTable = nullptr;

#define ONE_TIME_MATERIALS_LOCATION_GAMMA std::string("/mnt/sumit/PIR_GAMMA/ONE_TIME_MATERIALS/")
#define PER_EPOCH_MATERIALS_LOCATION_GAMMA std::string("/mnt/sumit/PIR_GAMMA/PER_EPOCH_MATERIALS/")
#define DATABASE_LOCATION_GAMMA std::string("/mnt/sumit/PIR_GAMMA/")
std::string L_filename = PER_EPOCH_MATERIALS_LOCATION_GAMMA+"L_gamma.bin";
std::string DK_filename = PER_EPOCH_MATERIALS_LOCATION_GAMMA+"DK_gamma.bin";//TODO: The key and data are seperated, unlike the description of the paper
std::string sdb_filename = PER_EPOCH_MATERIALS_LOCATION_GAMMA+"ShuffledDB_gamma.bin";
std::string HTable_filename = PER_EPOCH_MATERIALS_LOCATION_GAMMA+"H_gamma.bin";

// Function declarations
static int InitSrv_gamma();
static int OneTimeInit_gamma();
static int SelShuffDBSearchTag_gamma();
static int PerEpochOperations_gamma();
static int ProcessClientRequest_gamma();

static int FinSrv_gamma();

static void TestSrv_gamma();
static void TestPKEOperations_gamma();
static int TestHTableSerDser_gamma();
static void TestSelShuffDBSearchTag_gamma();


using namespace kuku;
static int Test_CuckooHash(table_size_type table_size, table_size_type stash_size, uint8_t loc_func_count, uint64_t max_probe);

// Function definitions
static int InitSrv_gamma(){
    int ret = -1;
    // Initialize random number generation
    std::random_device rd;
    unsigned long seed = (static_cast<unsigned long>(rd()) << 1) ^ rd();
    rng.seed(seed); // seed() seeds the gmp_randclass    
    
    InitConnectingSocket(SERVER_BETA_IP, BETA_LISTENING_TO_GAMMA_PORT, &sock_gamma_to_beta);
    PrintLog(LOG_LEVEL_TRACE, __FILE__, __LINE__, "Established connection with Server Beta");

    //Initialize sockets for communication with server alpha
    ret = InitAcceptingSocket(GAMMA_LISTENING_TO_ALPHA_PORT, &sock_gamma_to_alpha, &sock_gamma_to_alpha_con);

    if (ret != 0) {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Cannot open Accepting socket for Server Gamma!!");
        ret = -1;
        goto exit;
    }
    
    PrintLog(LOG_LEVEL_INFO, __FILE__, __LINE__, "Server Gamma initialization complete");

exit:
    if (ret != 0){
        FinSrv_gamma();
    }

    return ret;
}

static int OneTimeInit_gamma() {
    size_t received_sz = 0;
    int ret_recv = 0;

    // Receive all the parameters from server beta
    // Receive p
    ret_recv = recvAll(sock_gamma_to_beta, net_buf, sizeof(net_buf), &received_sz);
    if (ret_recv != 0)
    {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Failed to receive p from Server Beta");
        return -1;
    }
    p = mpz_class(std::string(net_buf, received_sz));

    // Receive q
    ret_recv = recvAll(sock_gamma_to_beta, net_buf, sizeof(net_buf), &received_sz);
    if (ret_recv != 0)
    {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Failed to receive q from Server Beta");
        return -1;
    }
    q = mpz_class(std::string(net_buf, received_sz));

    // Receive g
    ret_recv = recvAll(sock_gamma_to_beta, net_buf, sizeof(net_buf), &received_sz);
    if (ret_recv != 0)
    {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Failed to receive g from Server Beta");
        return -1;
    }
    g = mpz_class(std::string(net_buf, received_sz));

    // Receive g_q
    ret_recv = recvAll(sock_gamma_to_beta, net_buf, sizeof(net_buf), &received_sz);
    if (ret_recv != 0)
    {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Failed to receive g_q from Server Beta");
        return -1;
    }
    g_q = mpz_class(std::string(net_buf, received_sz));

    // Receive r
    ret_recv = recvAll(sock_gamma_to_beta, net_buf, sizeof(net_buf), &received_sz);
    if (ret_recv != 0)
    {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Failed to receive r from Server Beta");
        return -1;
    }
    r = mpz_class(std::string(net_buf, received_sz));

    // Receive pk_E
    ret_recv = recvAll(sock_gamma_to_beta, net_buf, sizeof(net_buf), &received_sz);
    if (ret_recv != 0)
    {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Failed to receive pk_E from Server Beta");
        return -1;
    }
    pk_E = mpz_class(std::string(net_buf, received_sz));

    // Receive pk_E_q
    ret_recv = recvAll(sock_gamma_to_beta, net_buf, sizeof(net_buf), &received_sz);
    if (ret_recv != 0)
    {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Failed to receive pk_E_q from Server Beta");
        return -1;
    }
    pk_E_q = mpz_class(std::string(net_buf, received_sz));

    // Receive FHEcryptoContext
    ret_recv = recvAll(sock_gamma_to_beta, net_buf, sizeof(net_buf), &received_sz);
    if (ret_recv != 0)
    {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Failed to receive FHEcryptoContext from Server Beta");
        return -1;
    }
    Serial::DeserializeFromString(FHEcryptoContext, std::string(net_buf, received_sz));

    // Receive pk_F
    ret_recv = recvAll(sock_gamma_to_beta, net_buf, sizeof(net_buf), &received_sz);
    if (ret_recv != 0)
    {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Failed to receive pk_F from Server Beta");
        return -1;
    }
    Serial::DeserializeFromString(pk_F, std::string(net_buf, received_sz));

    // Receive vectorOnesforElement_ct
    ret_recv = recvAll(sock_gamma_to_beta, net_buf, sizeof(net_buf), &received_sz);
    if (ret_recv != 0)
    {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Failed to receive vectorOnesforElement_ct from Server Beta");
        return -1;
    }
    Serial::DeserializeFromString(vectorOnesforElement_ct, std::string(net_buf, received_sz));

    // Receive vectorOnesforTag_ct
    ret_recv = recvAll(sock_gamma_to_beta, net_buf, sizeof(net_buf), &received_sz);
    if (ret_recv != 0)
    {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Failed to receive vectorOnesforTag_ct from Server Beta");
        return -1;
    }
    Serial::DeserializeFromString(vectorOnesforTag_ct, std::string(net_buf, received_sz));

    //Save parameters to local files
    export_to_file_from_mpz_class(ONE_TIME_MATERIALS_LOCATION_GAMMA + "p.bin", p);
    export_to_file_from_mpz_class(ONE_TIME_MATERIALS_LOCATION_GAMMA + "q.bin", q);
    export_to_file_from_mpz_class(ONE_TIME_MATERIALS_LOCATION_GAMMA + "g.bin", g);
    export_to_file_from_mpz_class(ONE_TIME_MATERIALS_LOCATION_GAMMA + "g_q.bin", g_q);
    export_to_file_from_mpz_class(ONE_TIME_MATERIALS_LOCATION_GAMMA + "r.bin", r);
    export_to_file_from_mpz_class(ONE_TIME_MATERIALS_LOCATION_GAMMA + "pk_E.bin", pk_E);
    export_to_file_from_mpz_class(ONE_TIME_MATERIALS_LOCATION_GAMMA + "pk_E_q.bin", pk_E_q);

    Serial::SerializeToFile(ONE_TIME_MATERIALS_LOCATION_GAMMA + "FHEcryptoContext.bin", FHEcryptoContext, SerType::BINARY);
    Serial::SerializeToFile(ONE_TIME_MATERIALS_LOCATION_GAMMA + "pk_F.bin", pk_F, SerType::BINARY);
    Serial::SerializeToFile(ONE_TIME_MATERIALS_LOCATION_GAMMA + "vectorOnesforElement_ct.bin", vectorOnesforElement_ct, SerType::BINARY);
    Serial::SerializeToFile(ONE_TIME_MATERIALS_LOCATION_GAMMA + "vectorOnesforTag_ct.bin", vectorOnesforTag_ct, SerType::BINARY);

    PrintLog(LOG_LEVEL_TRACE, __FILE__, __LINE__, "Received all the one-time initialized parameters from Server Beta and exported all of them into file");

    return 0;
}

/* TODO: This function, currently does not matches with the diagram of the paper. It has to be updated or the diagram has to be modified. */
static int PerEpochOperations_gamma(){
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

    PrintLog(LOG_LEVEL_INFO, __FILE__, __LINE__, "Server Gamma: Starting PerEpochOperations sequence");

    /* First of all retrieve all the one-time initialized materials from the saved location */
    p = import_from_file_to_mpz_class(ONE_TIME_MATERIALS_LOCATION_GAMMA + "p.bin");
    q = import_from_file_to_mpz_class(ONE_TIME_MATERIALS_LOCATION_GAMMA + "q.bin");
    g = import_from_file_to_mpz_class(ONE_TIME_MATERIALS_LOCATION_GAMMA + "g.bin");
    g_q = import_from_file_to_mpz_class(ONE_TIME_MATERIALS_LOCATION_GAMMA + "g_q.bin");
    r = import_from_file_to_mpz_class(ONE_TIME_MATERIALS_LOCATION_GAMMA + "r.bin");
    pk_E = import_from_file_to_mpz_class(ONE_TIME_MATERIALS_LOCATION_GAMMA + "pk_E.bin");
    pk_E_q = import_from_file_to_mpz_class(ONE_TIME_MATERIALS_LOCATION_GAMMA + "pk_E_q.bin");
    Serial::DeserializeFromFile(ONE_TIME_MATERIALS_LOCATION_GAMMA + "FHEcryptoContext.bin", FHEcryptoContext, SerType::BINARY);
    Serial::DeserializeFromFile(ONE_TIME_MATERIALS_LOCATION_GAMMA + "pk_F.bin", pk_F, SerType::BINARY);
    Serial::DeserializeFromFile(ONE_TIME_MATERIALS_LOCATION_GAMMA + "vectorOnesforElement_ct.bin", vectorOnesforElement_ct, SerType::BINARY);
    Serial::DeserializeFromFile(ONE_TIME_MATERIALS_LOCATION_GAMMA + "vectorOnesforTag_ct.bin", vectorOnesforTag_ct, SerType::BINARY);

    PrintLog(LOG_LEVEL_TRACE, __FILE__, __LINE__, "Server Gamma: Loaded one-time initialization materials");

    /* Wait for receiving the ready message from server-beta */
    (void)recvAll(sock_gamma_to_beta, net_buf, sizeof(net_buf), &received_sz);
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
    (void)recvAll(sock_gamma_to_beta, net_buf, sizeof(net_buf), &received_sz);
    
    if (ret != 0)
    {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Failed to receive COMPLETED_REINIT_FOR_EPOCH message from Server Beta");
        return -1;
    }

    if (std::string(net_buf, received_sz) != completed_reinit_for_epoch_message) {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Did not receive expected COMPLETED_REINIT_FOR_EPOCH message from Server Beta");
        return -1;
    } else {
        PrintLog(LOG_LEVEL_INFO, __FILE__, __LINE__, "Server Gamma: Completed re-initialization for new epoch, now ready to process client-requests..!!");
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

    /* 13.c.1 This file stores the pair: (Kuku key, the share of the element) */
    sdb.open(sdb_filename, std::ios::in | std::ios::out | std::ios::binary | std::ios::trunc);
    if (!sdb) {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Failed to open SDB file at location: " + sdb_filename);
        ret = -1;
        goto exit;
    }

    // 12.c.1 Allocate a new Cuckoo hash table
    // Keeping the number of entries, larger than the number of elements to place. The reason is, it will reduce the number of probe during placement and make the per epoch operations faster
    HTable = new KukuTable(CUCKOO_TABLE_SIZE, CUCKOO_STASH_SIZE, CUCKOO_LOC_FUNC_COUNT, CUCKOO_LOC_FUNC_SEED, CUCKOO_MAX_PROBE, CUCKOO_EMPTY_ITEM);
    if (HTable == nullptr) {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Failed to allocate memory for Cuckoo hash table");
        ret = -1;
        goto exit;
    }

    /* 12.c.2 Prepare the entire Kuku hash HTable, based on all the keys */
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
            PrintLog(LOG_LEVEL_TRACE, __FILE__, __LINE__, "Inserted " + to_string(i+1) + " items into the cuckoo hash table. Current stash size: " + to_string(HTable->stash().size()) + " and total probe count is: " + to_string (HTable->total_probe_count_));
        }
    }

    PrintLog(LOG_LEVEL_INFO, __FILE__, __LINE__, "Cuckoo hash HTable creation complete. Current stash size: " + to_string(HTable->stash().size()) + " and total probe count is: " + to_string (HTable->total_probe_count_));

    //13.c.2 Reset the read pointers to the beginning
    DK.seekg(0, std::ios::beg);
    L.seekg(0, std::ios::beg);
    sdb.seekp(0, std::ios::beg);

    //13.c.3 Now place all the elements from the temporary list(L) to the shuffled database(SDB) according to the Kuku hash table
    for (uint64_t i = 0; i < M; i++){
        /* Read the next key */
        DK.read(reinterpret_cast<char*>(Kuku_key.data()), sizeof(item_type));
        /* Read the next secret share */
        L.read(net_buf, NUM_BYTES_PER_SDB_ELEMENT);

        /* 13.c.4: Prepare them to a tuple of shuffled database */
        memcpy(sdb_entry.cuckoo_key.data(), Kuku_key.data(), sizeof(item_type));
        memcpy(sdb_entry.element, net_buf, NUM_BYTES_PER_SDB_ELEMENT);    

        /* 13.c.5: Query and find the location */
        res = HTable->query(Kuku_key);
        if (!res)
        {
            PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Query failed for the item number: " + to_string(i));
            ret = -1;
            goto exit;
        }
        else {
            /* 13.c.5 Insert at the location of the shuffled database, determined by the query result */
            insert_sdb_entry(sdb, res.location(), sdb_entry);
            PrintLog(LOG_LEVEL_TRACE, __FILE__, __LINE__, "Inserted item: " + to_string(i) + " of L to SDB at location: " + std::to_string(res.location()));
        }

        if (((i+1) % 100000000) == 0){
            PrintLog(LOG_LEVEL_TRACE, __FILE__, __LINE__, "Inserted " + to_string(i+1) + " items into the shuffled database");
        }
    }

    PrintLog(LOG_LEVEL_INFO, __FILE__, __LINE__, "Shuffled database creation complete");
    PrintLog(LOG_LEVEL_INFO, __FILE__, __LINE__, "TODO: Check whether ith item of DK and L are really placed in proper location of SDB");

    // 14.c. Clear the shelter count as well by setting it to zero
    K = 0;

    //TODO: Store the cuckoo table in the disk
    exportedHFile.open(HTable_filename, std::ios::binary);
    HTable->serialize(exportedHFile);
    exportedHFile.close();   

exit:
    L.close();
    DK.close();
    sdb.close();

    return ret;
}


//TODO static int ResumeRequestProcessing_gamma(){
static int ProcessClientRequest_gamma(){
    int ret = -1;
    struct sockaddr_in address;
    int opt = 1;
    int addrlen = sizeof(address);
    int accepted_socket = -1;
    size_t received_sz = 0;
    shuffled_db_entry sdb_entry;
    uint64_t M = (N + sqrt_N);
    item_type Kuku_key;    
    QueryResult res;
    std::fstream L;
    std::fstream DK;
    std::fstream sdb;
    std::ifstream importedHFile;    

    PrintLog(LOG_LEVEL_INFO, __FILE__, __LINE__, "Server Gamma: Starting Request processing sequence");

    /* First of all retrieve all the one-time initialized materials from the saved location */
    p = import_from_file_to_mpz_class(ONE_TIME_MATERIALS_LOCATION_GAMMA + "p.bin");
    q = import_from_file_to_mpz_class(ONE_TIME_MATERIALS_LOCATION_GAMMA + "q.bin");
    g = import_from_file_to_mpz_class(ONE_TIME_MATERIALS_LOCATION_GAMMA + "g.bin");
    g_q = import_from_file_to_mpz_class(ONE_TIME_MATERIALS_LOCATION_GAMMA + "g_q.bin");
    r = import_from_file_to_mpz_class(ONE_TIME_MATERIALS_LOCATION_GAMMA + "r.bin");
    pk_E = import_from_file_to_mpz_class(ONE_TIME_MATERIALS_LOCATION_GAMMA + "pk_E.bin");
    pk_E_q = import_from_file_to_mpz_class(ONE_TIME_MATERIALS_LOCATION_GAMMA + "pk_E_q.bin");
    Serial::DeserializeFromFile(ONE_TIME_MATERIALS_LOCATION_GAMMA + "FHEcryptoContext.bin", FHEcryptoContext, SerType::BINARY);
    Serial::DeserializeFromFile(ONE_TIME_MATERIALS_LOCATION_GAMMA + "pk_F.bin", pk_F, SerType::BINARY);
    Serial::DeserializeFromFile(ONE_TIME_MATERIALS_LOCATION_GAMMA + "vectorOnesforElement_ct.bin", vectorOnesforElement_ct, SerType::BINARY);
    Serial::DeserializeFromFile(ONE_TIME_MATERIALS_LOCATION_GAMMA + "vectorOnesforTag_ct.bin", vectorOnesforTag_ct, SerType::BINARY);

    PrintLog(LOG_LEVEL_TRACE, __FILE__, __LINE__, "Server Gamma: Loaded one-time initialization materials");

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
    PrintLog(LOG_LEVEL_TRACE, __FILE__, __LINE__, "Server Gamma: Loaded hash table into the RAM");

    /* TODO: Retrieve K from the serialized data */

    /* TODO: Load the shelter into the RAM */

    while (K < sqrt_N){
        PrintLog(LOG_LEVEL_INFO, __FILE__, __LINE__, "Waiting for a connection from the client..!!");

        ret = InitAcceptingSocket(GAMMA_LISTENING_TO_CLIENT_PORT, &sock_gamma_client_srv, &sock_gamma_client_con);

        if (ret != 0)
        {
            PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Cannot establish communication with the client!!");
            ret = -1;
            goto exit;
        }

        /* Close the connection with existing client */
        close(sock_gamma_client_srv);
        close(sock_gamma_client_con);
        K++;
        /* TODO: Store updated value of K in the disk */
    }

    PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Current epoch is completed. Please re-perform the per-epoch initialization");

exit:
    return ret;
}

static int FinSrv_gamma(){
    int ret = -1;

    // Close the sockets
    if (sock_gamma_to_beta != -1) {
        close(sock_gamma_to_beta);
        sock_gamma_to_beta = -1;
    }
    if (sock_gamma_to_alpha != -1) {
        close(sock_gamma_to_alpha);
        sock_gamma_to_alpha = -1;
    }
    if (sock_gamma_to_alpha_con != -1) {
        close(sock_gamma_to_alpha_con);
        sock_gamma_to_alpha_con = -1;
    }
    if (sock_gamma_client_srv != -1) {
        close(sock_gamma_client_srv);
        sock_gamma_client_srv = -1;
    }
    if (sock_gamma_client_con != -1) {
        close(sock_gamma_client_con);
        sock_gamma_client_con = -1;
    }

    PrintLog(LOG_LEVEL_SPECIAL, __FILE__, __LINE__, "Finalized Server Gamma");

    return ret;
}

static void TestPKEOperations_gamma(){
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

    (void)sendAll(sock_gamma_to_beta, m1.get_str().c_str(), m1.get_str().size());
    (void)sendAll(sock_gamma_to_beta, m2.get_str().c_str(), m2.get_str().size());
    (void)sendAll(sock_gamma_to_beta, m3.get_str().c_str(), m3.get_str().size());
    (void)sendAll(sock_gamma_to_beta, m4.get_str().c_str(), m4.get_str().size());
    (void)sendAll(sock_gamma_to_beta, c11.get_str().c_str(), c11.get_str().size());
    (void)sendAll(sock_gamma_to_beta, c12.get_str().c_str(), c12.get_str().size());
    (void)sendAll(sock_gamma_to_beta, c21.get_str().c_str(), c21.get_str().size());
    (void)sendAll(sock_gamma_to_beta, c22.get_str().c_str(), c22.get_str().size());
    (void)sendAll(sock_gamma_to_beta, c31.get_str().c_str(), c31.get_str().size());
    (void)sendAll(sock_gamma_to_beta, c32.get_str().c_str(), c32.get_str().size());
    (void)sendAll(sock_gamma_to_beta, c41.get_str().c_str(), c41.get_str().size());
    (void)sendAll(sock_gamma_to_beta, c42.get_str().c_str(), c42.get_str().size());
    (void)sendAll(sock_gamma_to_beta, tag.get_str().c_str(), tag.get_str().size());
    (void)sendAll(sock_gamma_to_beta, Serial::SerializeToString(ct_tag).c_str(), Serial::SerializeToString(ct_tag).size());

    return;
}

static void TestSrv_gamma(){
    //TestPKEOperations_gamma();
    //TestSelShuffDBSearchTag_gamma();
    //Nothing is required to be executed during measuring the performance of client(simulated)
    //TestHTableSerDser_gamma();
}


static void TestSelShuffDBSearchTag_gamma(){
    int ret = -1;

    ret = SelShuffDBSearchTag_gamma();

    PrintLog(LOG_LEVEL_TRACE, __FILE__, __LINE__, "Server Gamma: Ended SelShuffDBSearchTag sequence");
    PrintLog(LOG_LEVEL_INFO, __FILE__, __LINE__, "Retrieved tag T_*: " + T_star.get_str());   

    return;
}

ostream &operator<<(ostream &stream, item_type item)
{
    stream << item[1] << " " << item[0];
    return stream;
}

void print_stash(const KukuTable &HTable)
{
    for (table_size_type i = 0; i < HTable.stash().size(); i++)
    {
        const auto &item = HTable.stash(i);
        //PrintLog(LOG_LEVEL_INFO, __FILE__, __LINE__, "Stash item " + to_string(i) + ": " + get_high_word(item) + get_low_word(item));
    }
}

/* Copied from: https://github.com/microsoft/Kuku/blob/main/examples/example.cpp */
static int Test_CuckooHash(table_size_type table_size, table_size_type stash_size, uint8_t loc_func_count, uint64_t max_probe)
{
    unsigned int rehash_cnt = 0;
    KukuTable *HTable = nullptr;
    //unsigned long hash_table_sz = (N + sqrt_N);
    unsigned long hash_table_sz = table_size;
    item_type loc_func_seed;
    item_type empty_item = make_item(0, 0);
    uint64_t high_rand, low_rand;
    mpz_class rand_val_high, rand_val_low;
    unsigned long i;

    #define CUCKOO_HASH_TABLE_REHASH_TRY_COUNT 1

    PrintLog(LOG_LEVEL_INFO, __FILE__, __LINE__, "Starting to build Cuckoo hash table of size: " + to_string(hash_table_sz) + " with stash size: " + to_string(stash_size) + " with function count: " + to_string(loc_func_count)  + " and " + to_string(max_probe) + " number of maximum probes");

    while (rehash_cnt < CUCKOO_HASH_TABLE_REHASH_TRY_COUNT) // Limit the number of rehash attempts to avoid infinite loops
    {
        rand_val_high = rng.get_z_bits(64);
        rand_val_low = rng.get_z_bits(64);

        mpz_export(&high_rand, nullptr, 1, sizeof(high_rand), 0, 0, rand_val_high.get_mpz_t());
        mpz_export(&low_rand, nullptr, 1, sizeof(low_rand), 0, 0, rand_val_low.get_mpz_t());

        loc_func_seed = make_item(high_rand, low_rand);

        // Allocate a new Cuckoo hash table
        // KukuTable table(table_size, stash_size, loc_func_count, loc_func_seed, max_probe, empty_item);
        HTable = new KukuTable(hash_table_sz, stash_size, loc_func_count, loc_func_seed, max_probe, empty_item);

        for (i = 0; i < hash_table_sz; i++)
        {
            /* Create a random ||P||-bit tag for the moment */
            mpz_class tag = rng.get_z_bits(P_BITS);
            // Convert tag to string
            std::string tag_str = tag.get_str();

            // Hash with SHA-256
            unsigned char hash[SHA256_DIGEST_LENGTH];
            SHA256(reinterpret_cast<const unsigned char *>(tag_str.data()), tag_str.size(), hash);

            // Get first 128 bits as two 64-bit words
            uint64_t high = 0, low = 0;
            memcpy(&high, hash, 8);
            memcpy(&low, hash + 8, 8);

            // Create item_type key
            item_type key = make_item(high, low);

            // Insert into the cuckoo hash table
            if (!HTable->insert(key))
            {
                PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Insertion failed for the tag, having value:" + tag.get_str());
                PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Before failure, successfully inserted: " + to_string(i) + " out of " + to_string(hash_table_sz) + " items");
                PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "The size of the stash during failure: " + to_string(HTable->stash().size()));
                rehash_cnt++;
                PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Starting again, rehash count: " + to_string(rehash_cnt));

                /* Delete the already built HTable */
                delete HTable;
                HTable = nullptr;
                break;
            }
            if ((i > 0) && (i % 100000000 == 0)) {
                PrintLog(LOG_LEVEL_INFO, __FILE__, __LINE__, "Inserted " + to_string(i) + " items so far");
            }
        }

        if (i == hash_table_sz) {
            PrintLog(LOG_LEVEL_INFO, __FILE__, __LINE__, "Successfully inserted all " + to_string(hash_table_sz) + " items, to the Cuckoo hash table, after: " + to_string(rehash_cnt) + " rehash attempts, and the stash size is: " + to_string(HTable->stash().size()));
            break; // Successfully inserted all items
        }
    }

    return 0;
}


int main(int argc, char *argv[])
{
    int ret = -1;

    /* Perform the basic initialization */
    InitSrv_gamma();

    /* Process as per the command line arguments */
    if (argc >= 2) {
        if (std::string("one_time_init").compare(std::string(argv[1]))==0) {
            // Perform one-time initialization for server gamma
            ret = OneTimeInit_gamma();
        } else if (std::string("per_epoch_operations").compare(std::string(argv[1]))==0) {
            // Perform per-epoch initialization for server gamma
            ret = PerEpochOperations_gamma();
        } else if (std::string("clear_epoch_state").compare(std::string(argv[1]))==0) {
            // Clear the existing state of current epoch, start as if this is the first request of the epoch
            // Delete shelter content and set K = 0
        } else if (std::string("process_request").compare(std::string(argv[1]))==0) {
            // Start from last saved state
            ret = ProcessClientRequest_gamma();
        } else {
            PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Unknown command line argument:"+ std::string(argv[1]));
            PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Improper command line arguments. Usage: server_gamma <one_time_init|per_epoch_operations|clear_epoch_state|continue>");
        }
    } else {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Improper command line arguments. Usage: server_gamma <one_time_init|per_epoch_operations|clear_epoch_state|continue>");
    }

    if (ret == 0) {
        TestSrv_gamma();
    }

    FinSrv_gamma();

    return 0;
}

static int TestHTableSerDser_gamma(){
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