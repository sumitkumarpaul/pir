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
#define NUM_CPU_CORES 16
//#define SHELTER_STORING_LOCATION std::string("/mnt/sumit/dummy_shelter/")
#define SHELTER_STORING_LOCATION std::string("/dev/shm/")

// The value of N will determine the bitlength during the client initialization
#define NUM_TAG_BITS 3072 // 16 bits can represent up to 65536, which is more than enough for N=50000
#define B 512 // Block size in bits, can be adjusted as needed
// And number of bits determine the evalution time drastically
static shelter_element sh[sqrt_N];
mpz_class T_star;
static uint64_t K; // Current number of entries in the shelter, or the number of processed requests
static KukuTable *HTable = nullptr;

static mpz_class c;
static std::fstream sdb;

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
static int ShelterTagDetermination_gamma();
static int ObliviouslySearchShelter_gamma();
static int FetchCombineSelect_gamma();
static int ShelterUpdate_gamma();
static int ObliDecReturn_gamma();

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
        //memcpy(sdb_entry.cuckoo_key.data(), Kuku_key.data(), sizeof(item_type)); cuckoo key is no more present in the structure
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

    // 14.c. Nothing is required to be done for clearing the shelter content

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

static int ShelterTagDetermination_gamma(){
    int ret = 0;
    size_t received_sz = 0;
    int ret_recv = 0;    
    std::pair<mpz_class, mpz_class> E_g_pow_Rho_pow_I__mul_a;
    std::pair<mpz_class, mpz_class> E_g_pow_Rho_pow_I__mul_a_mul_c;
    std::pair<mpz_class, mpz_class> E_c;

    // Step 9.3.2 Receive the first component of E_g_pow_Rho_pow_I__mul_a
    PrintLog(LOG_LEVEL_TRACE, __FILE__, __LINE__, "Waiting to receive data from server Alpha on socket: " + std::to_string(sock_gamma_to_alpha));
    ret_recv = recvAll(sock_gamma_to_alpha_con, net_buf, sizeof(net_buf), &received_sz);
    if (ret_recv != 0)
    {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Failed to receive E_g_pow_Rho_pow_I__mul_a.first from the Server Alpha");
        return -1;
    }
    PrintLog(LOG_LEVEL_TRACE, __FILE__, __LINE__, "Server Gamma starts client request processing from this point");

    E_g_pow_Rho_pow_I__mul_a.first = mpz_class(std::string(net_buf, received_sz));

    // Step 9.4.2 Receive the second component of E_g_pow_Rho_pow_I__mul_a
    ret_recv = recvAll(sock_gamma_to_alpha_con, net_buf, sizeof(net_buf), &received_sz);
    if (ret_recv != 0)
    {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Failed to receive E_g_pow_Rho_pow_I__mul_a.second from the Server Alpha");
        return -1;
    }
    E_g_pow_Rho_pow_I__mul_a.second = mpz_class(std::string(net_buf, received_sz));

    // Step 10.1 Compute the ciphertext of c
    E_c = ElGamal_encrypt(c, pk_E);

    // Step 10.2. Multiply homomorphically
    E_g_pow_Rho_pow_I__mul_a_mul_c = ElGamal_mult_ct(E_g_pow_Rho_pow_I__mul_a, E_c);

    // Step 10.3.1 Send the first part to the server Beta
    (void)sendAll(sock_gamma_to_beta, E_g_pow_Rho_pow_I__mul_a_mul_c.first.get_str().c_str(), E_g_pow_Rho_pow_I__mul_a_mul_c.first.get_str().size());

    // Step 10.4.1 Send the second part to the server Beta
    (void)sendAll(sock_gamma_to_beta, E_g_pow_Rho_pow_I__mul_a_mul_c.second.get_str().c_str(), E_g_pow_Rho_pow_I__mul_a_mul_c.second.get_str().size());

    /* This is only for experimentation purpose */
    //PrintLog(LOG_LEVEL_TRACE, __FILE__, __LINE__, "Chosen c is: " + c.get_str());

    return ret;
}


static int ObliviouslySearchShelter_gamma() {
    // Set up variables
    Fss fServer;
    ServerKeyEq K_gamma;
    int ret = 0;
    size_t received_sz = 0;
    size_t dserializedFssSize;
    Ciphertext<DCRTPoly> fnd_gamma_ct_element, fnd_gamma_ct_tag;
    mpz_class d_ct_gamma = 0;
    std::vector<bool> thread_fnd(NUM_CPU_CORES, false);
    bool fnd_gamma = false;
    std::vector<mpz_class> thread_sums(NUM_CPU_CORES);

    // First, receive sk_F from the server Beta
    ret = recvAll(sock_gamma_to_beta, net_buf, sizeof(net_buf), &received_sz);
    if (ret != 0)
    {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Failed to receive sk_F from Server Beta");
        return -1;
    }

    dserializedFssSize = deserializeFssAndServerKeyEq(net_buf, received_sz, fServer, K_gamma);

#if TEST_SHELTER_FOUND
    mpz_class special_tag, special_tag_location;

    PrintLog(LOG_LEVEL_SPECIAL, __FILE__, __LINE__, "Enter the value of the set search tag (base 10): ");
    mpz_inp_str(special_tag.get_mpz_t(), stdin, 10);

    PrintLog(LOG_LEVEL_SPECIAL, __FILE__, __LINE__, "Enter the index within the shelter, where this special tag must be placed (set the same value as specified in server_alpha): ");
    mpz_inp_str(special_tag_location.get_mpz_t(), stdin, 10);

    sh[special_tag_location.get_ui()].tag_short = special_tag;

    if (special_tag_location.get_ui() >= K) {
        PrintLog(LOG_LEVEL_SPECIAL, __FILE__, __LINE__, "Since the entered position is greater than the current size of the shelter, there will not be any shelter hit.");
    }
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
                if (evaluateEq(&fServer, &K_gamma, sh[k + j].tag_short)) {
                    mpz_xor(thread_sums[j].get_mpz_t(), thread_sums[j].get_mpz_t(), sh[k+j].element_FHE_ct.get_mpz_t());

                    /* Same as XORing */
                    thread_fnd[j] = !thread_fnd[j];
                }
                else{
                }
            }
        }
        for (int t = 0; t < NUM_CPU_CORES; ++t)
        {
            mpz_xor(d_ct_gamma.get_mpz_t(), d_ct_gamma.get_mpz_t(), thread_sums[t].get_mpz_t());
        }
    }
    for (int t = 0; t < NUM_CPU_CORES; ++t)
    {
        fnd_gamma ^= thread_fnd[t];
    }

    PrintLog(LOG_LEVEL_TRACE, __FILE__, __LINE__, "Completed DPF evaluation. Value of fnd_gamma: " + std::to_string(fnd_gamma));

    // Step 7.1 Compute
    if (fnd_gamma == true){
        FHE_EncOfOnes(fnd_gamma_ct_element, fnd_gamma_ct_tag);
    }else{
        FHE_EncOfZeros(fnd_gamma_ct_element, fnd_gamma_ct_tag);
    }

    // Step 7.2.1 Send ciphertext of fnd_gamma_ct_element
    (void)sendAll(sock_gamma_to_alpha_con, Serial::SerializeToString(fnd_gamma_ct_element).c_str(), Serial::SerializeToString(fnd_gamma_ct_element).size());
    
    // Step 7.3.1 Send ciphertext of fnd_gamma_ct_tag
    (void)sendAll(sock_gamma_to_alpha_con, Serial::SerializeToString(fnd_gamma_ct_tag).c_str(), Serial::SerializeToString(fnd_gamma_ct_tag).size());

    // Step 7.4.1 Send the computed share of the search result
    (void)sendAll(sock_gamma_to_alpha_con, d_ct_gamma.get_str().c_str(), d_ct_gamma.get_str().size());

    return 0;
}

static int FetchCombineSelect_gamma(){
    int ret = 0;
    size_t received_sz = 0;
    int ret_recv = 0;
    QueryResult Qres;
    item_type Kuku_key;
    shuffled_db_entry SR_D_gamma;
    uint64_t L_i;
    mpz_t tmp;
    mpz_class SR_D_gamma_mpz;
    mpz_init(tmp);
    Ciphertext<DCRTPoly> SR_D_alpha_ct, SR_D_gamma_ct, SR_D_ct;

    /* 1.a.1 Convert T_* to cuckoo hash key */
    mpz_export(net_buf, NULL, 1, 1, 1, 0, T_star.get_mpz_t());
    // Create a temporary array to satisfy the function signature
    std::array<unsigned char, 16> temp;

    convert_buf_to_item_type2((const unsigned char *)net_buf, (P_BITS / 8), temp);
    // Copy the bytes into the Kuku_key variable
    std::memcpy(&Kuku_key, temp.data(), sizeof(Kuku_key));

    /* 1.a.2 Lookup the location according to the Kuku table */
    Qres = HTable->query(Kuku_key);
    if (!Qres)
    {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Query failed for the request number: " + to_string(K));
        ret = -1;
        goto exit;
    }

    /* 1.a.3 Read the entry from that location of the shuffled database */
    L_i = Qres.location();
    read_sdb_entry(sdb, L_i, SR_D_gamma);
    PrintLog(LOG_LEVEL_INFO, __FILE__, __LINE__, "SDB touch location: " + std::to_string(L_i));

    /* 2.4 Receive the FHE ciphertext SR_D_alpha_ct  */
    ret = recvAll(sock_gamma_to_alpha_con, net_buf, sizeof(net_buf), &received_sz);
    if (ret != 0)
    {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Failed to receive FHE Ciphertext SR_D_gamma_ct from Server Alpha");
        goto exit;
    }
    Serial::DeserializeFromString(SR_D_alpha_ct, std::string(net_buf, received_sz));

    /* 3.1 First convert from shuffled_db_entry to mpz_class */
    mpz_import(tmp, NUM_BYTES_PER_SDB_ELEMENT, 1, 1, 1, 0, SR_D_gamma.element);
    SR_D_gamma_mpz = mpz_class(tmp);
    
    /* 3.2 Compute SR_D_gamma_ct */
    SR_D_gamma_ct = FHE_Enc_SDBElement(SR_D_gamma_mpz);

    PrintLog(LOG_LEVEL_TRACE, __FILE__, __LINE__, "SR_D_gamma_mpz is: " + SR_D_gamma_mpz.get_str(16));
    
    /* 4. Homomorphically combine them */
    SR_D_ct = (SR_D_alpha_ct + SR_D_gamma_ct);

    /* 5.1.2 Receive fnd_ct_element  */
    ret = recvAll(sock_gamma_to_alpha_con, net_buf, sizeof(net_buf), &received_sz);
    if (ret != 0)
    {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Failed to receive FHE Ciphertext fnd_ct_element from Server Alpha");
        goto exit;
    }
    Serial::DeserializeFromString(fnd_ct_element, std::string(net_buf, received_sz));

    /* 5.2.2 Receive SR_sh_ct  */
    ret = recvAll(sock_gamma_to_alpha_con, net_buf, sizeof(net_buf), &received_sz);
    if (ret != 0)
    {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Failed to receive FHE Ciphertext SR_sh_ct from Server Alpha");
        goto exit;
    }
    Serial::DeserializeFromString(SR_sh_ct, std::string(net_buf, received_sz));

    /* 6. Select the ciphertext of the requested element */
    requested_element_ct = FHE_SelectElement(fnd_ct_element, SR_D_ct, SR_sh_ct);


    /* Updated flow to cope up with cihpetext refresh related modification.
       Moved the step 7.1 of sending requested_element_ct to server_Alpha
       from here to the first step of the ShelterUpdate_gamma() function. */

exit:

    return ret;
}

static int ObliDecReturn_gamma(){
    int ret = -1;
    size_t received_sz = 0;
    int ret_recv = 0;
    Ciphertext<DCRTPoly> m_C_ct, masked_requested_element_ct;

    /* Step 2.3: Receive the ciphertext of the mask */
    ret_recv = recvAll(sock_gamma_client_con, net_buf, sizeof(net_buf), &received_sz);
    if (ret_recv != 0)
    {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Failed to receive m_C_ct from Client");
        goto exit;
    }
    Serial::DeserializeFromString(m_C_ct, std::string(net_buf, received_sz));

    /* Step 3: Homomorphically apply the mask. In paper it is mentioned +, but here we are using - */
    /* Consequently, the client is using +, instead of - */
    masked_requested_element_ct = requested_element_ct - m_C_ct;

    /* Step 4.1: Send the masked result to server beta */
    (void)sendAll(sock_gamma_to_beta, Serial::SerializeToString(masked_requested_element_ct).c_str(), Serial::SerializeToString(masked_requested_element_ct).size());

    /* Additional step */
    ret_recv = recvAll(sock_gamma_client_con, net_buf, sizeof(net_buf), &received_sz);
    if (ret_recv != 0)
    {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Failed to receive refreshed requested_element_ct from Client");
        goto exit;
    }
    Serial::DeserializeFromString(requested_element_ct, std::string(net_buf, received_sz));    

    ret = 0;

exit:

    return ret;
}


static int ShelterUpdate_gamma(){
    int ret = 0;
    size_t received_sz = 0;
    int ret_recv = 0;
    mpz_class c_dashed, h_gamma0, c_dashed_h_gamma0, T_star_a_dashed_b_dashed_c_dashed_h_gamma0, h_gamma0_1, T_star_hat, t_star_hat, h_alpha3, Del_c, c_1, Del_a_Del_b_Del_c_h_alpha3, h_alpha3_1, Del_a_Del_b_Del_c;
    std::pair<mpz_class, mpz_class> E_T_star_a_dashed, E_c_dashed_h_gamma0, E_T_star_a_dashed_c_dashed_h_gamma0, E_Del_a_h_alpha3, E_Del_a_Del_c_h_alpha3, E_Del_c;

    /* !!!!! [Updated flow to refresh ciphertext] This was actually step 7.1 of FetchCombineSelect_alpha in the diagram */
    (void)sendAll(sock_gamma_to_alpha_con, Serial::SerializeToString(requested_element_ct).c_str(), Serial::SerializeToString(requested_element_ct).size());

start:
    /* 1.c.1 Randomly select c_dashed */
    c_dashed = ElGamal_randomGroupElement();

    /* 1.c.2 Randomly select h_gamma0 */
    h_gamma0 = ElGamal_randomGroupElement();

    /* 2.2.2 Receive first component of E(T_*.a_dashed) */
    ret = recvAll(sock_gamma_to_alpha_con, net_buf, sizeof(net_buf), &received_sz);
    if (ret != 0) {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Failed to receive first component of E(T_*.a_dashed) from Server Alpha");
        close(sock_gamma_to_alpha_con);
        goto exit;
    }
    E_T_star_a_dashed.first = mpz_class(std::string(net_buf, received_sz));

    /* 2.3.2 Receive second component of E(T_*.a_dashed) */
    ret = recvAll(sock_gamma_to_alpha_con, net_buf, sizeof(net_buf), &received_sz);
    if (ret != 0) {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Failed to receive second component of E(T_*.a_dashed) from Server Alpha");
        close(sock_gamma_to_alpha_con);
        goto exit;
    }
    E_T_star_a_dashed.second = mpz_class(std::string(net_buf, received_sz));

    /* 3.1 Compute E(c'.h_gamma0) */
    c_dashed_h_gamma0 = ((c_dashed * h_gamma0) % p);
    E_c_dashed_h_gamma0 = ElGamal_encrypt(c_dashed_h_gamma0, pk_E);

    /* 3.2 Compute E(T_*.a'.c'.h_gamma0) */
    E_T_star_a_dashed_c_dashed_h_gamma0 = ElGamal_mult_ct(E_T_star_a_dashed, E_c_dashed_h_gamma0);

    /* 3.3.1 Send the first componet of E(T_*.a'.c'.h_gamma0) to Server beta */
    (void)sendAll(sock_gamma_to_beta, E_T_star_a_dashed_c_dashed_h_gamma0.first.get_str().c_str(), E_T_star_a_dashed_c_dashed_h_gamma0.first.get_str().size());

    /* 3.4.1 Send the second componet of E(T_*.a'.c'.h_gamma0) to Server beta */
    (void)sendAll(sock_gamma_to_beta, E_T_star_a_dashed_c_dashed_h_gamma0.second.get_str().c_str(), E_T_star_a_dashed_c_dashed_h_gamma0.second.get_str().size());

    /* 4.3.2 Receive T_*.a'.b'.c'.h_gamma0 from the server beta */
    ret = recvAll(sock_gamma_to_beta, net_buf, sizeof(net_buf), &received_sz);
    if (ret != 0) {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Failed to receive T_*.a'.b'.c'.h_gamma0 from Server Beta");
        close(sock_gamma_to_beta);
        goto exit;
    }

    T_star_a_dashed_b_dashed_c_dashed_h_gamma0 = mpz_class(std::string(net_buf, received_sz));

    /* 5.1 Compute h_gamma0^{-1}  */
    mpz_invert(h_gamma0_1.get_mpz_t(), h_gamma0.get_mpz_t(), p.get_mpz_t());
    
    /* 5.2 Compute T_star^  */
    T_star_hat = (T_star_a_dashed_b_dashed_c_dashed_h_gamma0 * h_gamma0_1) % p;

    /* 5.3 Compute t_star^  */
    t_star_hat = (T_star_hat % r);
    
    /* 6.1.1 Send T_star_hat to Server Alpha */
    (void)sendAll(sock_gamma_to_alpha_con, T_star_hat.get_str().c_str(), T_star_hat.get_str().size());

    /* 6.2.1 Send t_star_hat to Server Alpha */
    (void)sendAll(sock_gamma_to_alpha_con, t_star_hat.get_str().c_str(), t_star_hat.get_str().size());

    /* 7.c.1 Convert the ciphertext of the requested element in mpz format and then append that to the shelter */
    if (Serial::SerializeToFile("/dev/shm/tmp.ct", requested_element_ct, SerType::BINARY) == true){
        sh[K].element_FHE_ct = import_from_file_to_mpz_class("/dev/shm/tmp.ct");        
    }else{
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Failed to convert and then save the ciphertext of the requested element in mpz format");
        goto exit;
    }
    /* 7.c.2 Append the shelter tags */
    sh[K].tag = T_star_hat;
    sh[K].tag_short = t_star_hat;

    
    /* 8.3.2 Receive h_alpha3 */
    ret = recvAll(sock_gamma_to_alpha_con, net_buf, sizeof(net_buf), &received_sz);
    if (ret != 0) {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Failed to receive h_alpha3 from Server Alpha");
        close(sock_gamma_to_alpha_con);
        goto exit;
    }
    h_alpha3 = mpz_class(std::string(net_buf, received_sz));

    /* 8.4.1.2 Receive first component of E(Del_a_h_alpha3) */
    ret = recvAll(sock_gamma_to_alpha_con, net_buf, sizeof(net_buf), &received_sz);
    if (ret != 0) {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Failed to receive first component of E(Del_a_h_alpha3) from Server Alpha");
        close(sock_gamma_to_alpha_con);
        goto exit;
    }
    E_Del_a_h_alpha3.first = mpz_class(std::string(net_buf, received_sz));

    /* 8.4.2.2 Receive second component of E(Del_a_h_alpha3) */
    ret = recvAll(sock_gamma_to_alpha_con, net_buf, sizeof(net_buf), &received_sz);
    if (ret != 0) {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Failed to receive second component of E(Del_a_h_alpha3) from Server Alpha");
        close(sock_gamma_to_alpha_con);
        goto exit;
    }
    E_Del_a_h_alpha3.second = mpz_class(std::string(net_buf, received_sz));

    /* 9.1 Compute c^{-1} */
    mpz_invert(c_1.get_mpz_t(), c.get_mpz_t(), p.get_mpz_t());

    /* 9.2 Compute c'.c^{-1} */
    Del_c = (c_dashed * c_1) % p;

    /* 9.3 Compute c'.c^{-1} */
    E_Del_c = ElGamal_encrypt(Del_c, pk_E);

    /* 9.4 Compute E(Del_a.Del_c.h_alpha3) */
    E_Del_a_Del_c_h_alpha3 = ElGamal_mult_ct(E_Del_a_h_alpha3, E_Del_c);

    /* 9.4.1. Send the first componet of E(Del_a.Del_c.h_alpha3) to Server Beta */
    (void)sendAll(sock_gamma_to_beta, E_Del_a_Del_c_h_alpha3.first.get_str().c_str(), E_Del_a_Del_c_h_alpha3.first.get_str().size());

    /* 9.5.1 Send the second componet of E(Del_a.Del_c.h_alpha3) to Server Beta */
    (void)sendAll(sock_gamma_to_beta, E_Del_a_Del_c_h_alpha3.second.get_str().c_str(), E_Del_a_Del_c_h_alpha3.second.get_str().size());    

    /* 11.2.2 Receive Del_a_Del_b_Del_c_h_alpha3 to Server Beta */
    ret = recvAll(sock_gamma_to_beta, net_buf, sizeof(net_buf), &received_sz);
    if (ret != 0) {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Failed to receive Del_a_Del_b_Del_c_h_alpha3 from Server Beta");
        close(sock_gamma_to_beta);
        goto exit;
    }

    Del_a_Del_b_Del_c_h_alpha3 = mpz_class(std::string(net_buf, received_sz));

    /* 12.c.1 Compute h_alpha3^{-1} */
    mpz_invert(h_alpha3_1.get_mpz_t(), h_alpha3.get_mpz_t(), p.get_mpz_t());

    /* 12.c.2 Extract Del_a_Del_b_Del_c */
    Del_a_Del_b_Del_c = ((Del_a_Del_b_Del_c_h_alpha3 * h_alpha3_1) % p);

    for (unsigned int i = 0; i < K; i++) {
        sh[i].tag = ((sh[i].tag * Del_a_Del_b_Del_c) % p);
        sh[i].tag_short = sh[i].tag % r;
        PrintLog(LOG_LEVEL_TRACE, __FILE__, __LINE__, "sh["+ std::to_string(i) + "].tag_short: " + sh[i].tag_short.get_str());
        
        if (t_star_hat == sh[i].tag_short) {
            PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Tag collision.. Restarting the shelter update process..!!");
            (void)sendAll(sock_gamma_to_beta, reinit_shelter_update_message.c_str(), reinit_shelter_update_message.size());
            goto start;
        }
    }

    (void)sendAll(sock_gamma_to_beta, completed_request_processing_message.c_str(), completed_request_processing_message.size());

    ret = 0;
    c = c_dashed;

exit:

    return ret;
}

//TODO static int ResumeRequestProcessing_gamma(){
static int ProcessClientRequest_gamma(){
    int ret = -1;
    uint64_t M = (N + sqrt_N);
    item_type Kuku_key;    
    QueryResult res;
    std::fstream L;
    std::fstream DK;
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
    
    //Always initialize them
    K = 0;
    c = 1;/* Initialize with 1, so that, during the first iteration c = c'*(c^-1) = c'*1 = c' */

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

    /* Open the shuffled and secret shared database */
    sdb.open(sdb_filename, std::ios::in | std::ios::binary);
    if (!sdb) {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Failed to open SDB file at location: " + sdb_filename);
        ret = -1;
        goto exit;
    }

    while (K < sqrt_N){
        PrintLog(LOG_LEVEL_INFO, __FILE__, __LINE__, "Waiting for processing the PIR request number: "+ std::to_string(K+1) +" from the client..!!");

        ret = InitAcceptingSocket(GAMMA_LISTENING_TO_CLIENT_PORT, &sock_gamma_client_srv, &sock_gamma_client_con);

        if (ret != 0)
        {
            PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Cannot establish communication with the client!!");
            ret = -1;
            goto exit;
        }

        /* Even if for the request number 1, the determined shelter tag is not required to be used.
           Still, we cannot move this function within if (K > 0) block, since the client interaction is involved in this function. */        
        ret = ShelterTagDetermination_gamma();
        if (ret != 0)
        {
            PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Problem while determining the shelter tag..!!");
            ret = -1;
            goto exit;
        }

        /* For the first request, the shelter is not required to be searched */
        if (K > 0){
            ret = ObliviouslySearchShelter_gamma();
            if (ret != 0)
            {
                PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Problem during the shelter search operation..!!");
                ret = -1;
                goto exit;
            }
        }

        ret = SelShuffDBSearchTag_gamma();

        if (ret != 0){
            PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Problem during the selecting shuffled database tag..!!");
            ret = -1;
            goto exit;
        }

        ret = FetchCombineSelect_gamma();
        if (ret != 0){
            PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Problem during the Fetch_Combine_and_Select stage..!!");
            ret = -1;
            goto exit;
        }

        /* !!!!! [Updated flow to refresh ciphertext] First returning the data and then updating the shelter */
        ret = ObliDecReturn_gamma();
        if (ret != 0){
            PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Problem during the Shelter Update stage..!!");
            ret = -1;
            goto exit;
        }

        
        ret = ShelterUpdate_gamma();
        if (ret != 0){
            PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Problem during the Shelter Update stage..!!");
            ret = -1;
            goto exit;
        }
        
        /* Close the connection with existing client */
        close(sock_gamma_client_srv);
        close(sock_gamma_client_con);
        K++;
    }

    PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Current epoch is completed. Please re-perform the per-epoch initialization");

exit:
    sdb.close();
    /* If there are any dangling connection with client, close that */
    close(sock_gamma_client_srv);
    close(sock_gamma_client_con);

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

    //ret = SelShuffDBSearchTag_gamma();

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

static int Perf_avg_online_server_time_gamma() {
    // Set up variables
    Fss fClient, fServer;
    ServerKeyEq k0;
    ServerKeyEq k1;
    int ret = 0;
    size_t received_sz = 0;
    /* On average half of the shelter elements will be populated */
    int average_shelter_size = (sqrt_N/2);
    std::string DPF_search_test_shelter_location = std::string("/dev/shm/");
    /* Generate a dummy delta value to update the shelter tags */
    mpz_class Del_abc = rng.get_z_bits(P_BITS);
    /* Suppose we want to search for a random tag */
    mpz_class tmp = rng.get_z_range(average_shelter_size);
    uint64_t dpf_random_test_index = tmp.get_ui();
    mpz_class T_sh_short = sh[dpf_random_test_index].tag_short;


    /* First of all retrieve all the one-time initialized materials from the saved location */
    p = import_from_file_to_mpz_class(ONE_TIME_MATERIALS_LOCATION_GAMMA + "p.bin");
    q = import_from_file_to_mpz_class(ONE_TIME_MATERIALS_LOCATION_GAMMA + "q.bin");
    g = import_from_file_to_mpz_class(ONE_TIME_MATERIALS_LOCATION_GAMMA + "g.bin");
    r = import_from_file_to_mpz_class(ONE_TIME_MATERIALS_LOCATION_GAMMA + "r.bin");
    Serial::DeserializeFromFile(ONE_TIME_MATERIALS_LOCATION_GAMMA + "FHEcryptoContext.bin", FHEcryptoContext, SerType::BINARY);
    Serial::DeserializeFromFile(ONE_TIME_MATERIALS_LOCATION_GAMMA + "pk_F.bin", pk_F, SerType::BINARY);

    PrintLog(LOG_LEVEL_TRACE, __FILE__, __LINE__, "Starting to randomly populate a shelter of average size: " + to_string(average_shelter_size));

    /* Populate the shelter, with random elements */
    for(size_t k = 0; k < average_shelter_size; k++) {
        // Generate random block_content of PLAINTEXT_PIR_BLOCK_DATA_SIZE bits of random | k as the block index
        Ciphertext<DCRTPoly> tmp_ct = FHE_Enc_SDBElement((rng.get_z_bits(PLAINTEXT_PIR_BLOCK_DATA_SIZE) << log_N) | mpz_class(k));
        /* Store the ciphertexts to serialized form to a file, which resides in the RAM */
        if (Serial::SerializeToFile(DPF_search_test_shelter_location + "sh[" + std::to_string(k) + "].ct", tmp_ct, SerType::BINARY) != true)
        {
            PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Failed to serialize element FHE ciphertext to file");
        }

        sh[k].element_FHE_ct = import_from_file_to_mpz_class(DPF_search_test_shelter_location + "sh[" + std::to_string(k) + "].ct");

        /* Generate the tags and keep them in the variable, which will be used for DPF search */
        sh[k].tag = ElGamal_randomGroupElement(); // Create a random tag
        sh[k].tag_short = sh[k].tag % r; // Create a random short tag
    }

    PrintLog(LOG_LEVEL_TRACE, __FILE__, __LINE__, "Starting to test DPF-search on the shelter");

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
                if (evaluateEq(&fServer, &k0, sh[k + j].tag_short)) {
                    mpz_xor(thread_sums[j].get_mpz_t(), thread_sums[j].get_mpz_t(), sh[k+j].element_FHE_ct.get_mpz_t());
                }
                /* Simulate the time required for shelter tag update operation */
                sh[k + j].tag = (sh[k + j].tag * Del_abc) % p;
                sh[k + j].tag_short = sh[k + j].tag % r;                
            }
        }
        for (int t = 0; t < NUM_CPU_CORES; ++t){
            mpz_xor(ans0.get_mpz_t(), ans0.get_mpz_t(), thread_sums[t].get_mpz_t());
        }

    }

    PrintLog(LOG_LEVEL_TRACE, __FILE__, __LINE__, "Completed DPF evaluation over average number of shelter elements");


    return 1;
}

static int SelShuffDBSearchTag_gamma(){
    int ret = -1;
    size_t received_sz = 0;    

    PrintLog(LOG_LEVEL_TRACE, __FILE__, __LINE__, "Server Gamma: Starting SelShuffDBSearchTag sequence");

    // 8.c.1 Receive T_phi.h_{\\beta 0} from the server beta
    ret = recvAll(sock_gamma_to_beta, net_buf, sizeof(net_buf), &received_sz);

    if (ret != 0) {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Server Gamma: Failed to receive T_phi.h_{\\beta 0} from Server Beta");
        close(sock_gamma_to_beta);
        return ret;
    }

    mpz_class T_phi_h_beta_0 = mpz_class(std::string(net_buf, received_sz));

    // 9.c Receive h_{\\alpha 2} from the server alpha
    ret = recvAll(sock_gamma_to_alpha_con, net_buf, sizeof(net_buf), &received_sz);

    if (ret != 0) {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Server Gamma: Failed to receive h_{\\alpha 2} from Server Alpha");
        close(sock_gamma_to_alpha_con);
        return ret;
    }

    mpz_class h_alpha2 = mpz_class(std::string(net_buf, received_sz));

    /* 10.c.1 Compute T_phi.h_{\alpha 2}.h_{\beta 0} */
    mpz_class T_phi_h_alpha2_h_beta0 = (T_phi_h_beta_0*h_alpha2) % p;

    /* 10.c.2 FHE Encrypt T_phi.h_{\\alpha 2}.h_{\\beta 0} */
    Ciphertext<DCRTPoly> FHE_ct_T_phi_h_alpha2_h_beta0 = FHE_Enc_Tag(T_phi_h_alpha2_h_beta0);

    /* 10.c.3 Send the ciphertext to server alpha */
    (void)sendAll(sock_gamma_to_alpha_con, Serial::SerializeToString(FHE_ct_T_phi_h_alpha2_h_beta0).c_str(), Serial::SerializeToString(FHE_ct_T_phi_h_alpha2_h_beta0).size());

    // 14.c.1 Receive T_* from the server alpha
    ret = recvAll(sock_gamma_to_alpha_con, net_buf, sizeof(net_buf), &received_sz);

    if (ret != 0) {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Server Gamma: Failed to receive T_* from Server Alpha");
        close(sock_gamma_to_alpha_con);
        return ret;
    }

    T_star = mpz_class(std::string(net_buf, received_sz));
    PrintLog(LOG_LEVEL_INFO, __FILE__, __LINE__, "Selected T_* is: " + T_star.get_str());

    return ret;
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
        } else if (std::string("test").compare(std::string(argv[1]))==0) {
            TestSrv_gamma();
        } else if (std::string("perf").compare(std::string(argv[1]))==0) {
            if (argc < 3){
                PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Performance measurement option requires at least three command line parameters. Usage: server_gamma perf [srv_avg_online_time]");
            }else{
                if (std::string("srv_avg_online_time").compare(std::string(argv[2]))==0){
                    (void)Perf_avg_online_server_time_gamma();
                }
            }
        } else {
            PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Unknown command line argument:"+ std::string(argv[1]));
            PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Improper command line arguments. Usage: server_gamma <one_time_init|per_epoch_operations|clear_epoch_state|process_request>");
        }
    } else {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Improper command line arguments. Usage: server_gamma <one_time_init|per_epoch_operations|clear_epoch_state|process_request>");
    }

    if (ret == 0) {
        //TestSrv_gamma();
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