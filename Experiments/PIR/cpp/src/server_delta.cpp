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
#include <filesystem> // Required for std::filesystem

#include <iomanip>
#include <iostream>


#include <iterator>
#include <cstring>
#include <unistd.h>
#include "pir_common.h"

static int sock_delta_to_beta = -1, sock_delta_to_alpha = -1;
static int sock_alpha_client_srv = -1, sock_alpha_client_con = -1;
static char net_buf[NET_BUF_SZ] = {0};

#define B 512 // Block size in bits, can be adjusted as needed
// And number of bits determine the evalution time drastically
static uint64_t K; // Current number of entries in the shelter, or the number of processed requests

#define NUM_CPU_CORES 16

#define DPF_SEARCH_INDEX_K 1
#define SHELTER_STORING_LOCATION std::string("/dev/shm/")

#define ONE_TIME_MATERIALS_LOCATION_DELTA std::string("/mnt/sumit/PIR_DELTA/ONE_TIME_MATERIALS/")
#define PER_EPOCH_MATERIALS_LOCATION_DELTA std::string("/mnt/sumit/PIR_DELTA/PER_EPOCH_MATERIALS/")
#define MASK_LOCATION_DELTA std::string("/mnt/sumit/PIR_DELTA/")
std::string mdb_filename = PER_EPOCH_MATERIALS_LOCATION_DELTA+"MaskDB_alpha.bin";


// Function declarations
static int InitSrv_delta();
static int OneTimeInit_delta();
static int FinSrv_delta();
static int PerEpochOperations_delta();
static int ProcessClientRequest_delta();
static int ObliviouslySearchShelter_delta();
static void TestSrv_delta();

static int Perf_avg_online_server_time_delta();

// Function definitions
static int InitSrv_delta(){
    int ret = -1;
    // Initialize random number generation
    std::random_device rd;
    unsigned long seed = (static_cast<unsigned long>(rd()) << 1) ^ rd();
    rng.seed(seed); // seed() seeds the gmp_randclass    
    
    // Server_alpha only connects to other servers, it does not listen to other servers
    InitConnectingSocket(SERVER_BETA_IP, BETA_LISTENING_TO_DELTA_PORT, &sock_delta_to_beta);

    PrintLog(LOG_LEVEL_TRACE, __FILE__, __LINE__, "Established connection with Server Beta");

    InitConnectingSocket(SERVER_ALPHA_IP, ALPHA_LISTENING_TO_DELTA_PORT, &sock_delta_to_alpha);

    PrintLog(LOG_LEVEL_TRACE, __FILE__, __LINE__, "Established connection with Server Alpha");

    PrintLog(LOG_LEVEL_INFO, __FILE__, __LINE__, "Server Delta initialization complete");

    ret = 0;

exit:
    if (ret != 0){
        FinSrv_delta();
    }

    return ret;
}

static int OneTimeInit_delta() {
    size_t received_sz = 0;
    int ret_recv = 0;

    // Receive all the parameters from server beta
    // Receive FHEcryptoContext
    ret_recv = recvAll(sock_delta_to_beta, net_buf, sizeof(net_buf), &received_sz);
    if (ret_recv != 0)
    {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Failed to receive FHEcryptoContext from Server Beta");
        return -1;
    }
    Serial::DeserializeFromString(FHEcryptoContext, std::string(net_buf, received_sz));

    // Receive pk_F
    ret_recv = recvAll(sock_delta_to_beta, net_buf, sizeof(net_buf), &received_sz);
    if (ret_recv != 0)
    {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Failed to receive pk_F from Server Beta");
        return -1;
    }
    Serial::DeserializeFromString(pk_F, std::string(net_buf, received_sz));

    //Save parameters to local files
    if (!Serial::SerializeToFile(ONE_TIME_MATERIALS_LOCATION_DELTA + "FHEcryptoContext.bin", FHEcryptoContext, SerType::BINARY)){
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Failed to serialize the received FHEcryptoContext from Server Beta");
        return -1;
    }
    
    if (!Serial::SerializeToFile(ONE_TIME_MATERIALS_LOCATION_DELTA + "pk_F.bin", pk_F, SerType::BINARY)){
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Failed to serialize the received FHEcryptoContext from Server Beta");
        return -1;
    }

    PrintLog(LOG_LEVEL_TRACE, __FILE__, __LINE__, "Received all the one-time initialized parameters from Server Beta and exported all of them into file");
    
    return 0;
}

static int FinSrv_delta(){
    int ret = -1;

    // Close the sockets
    if (sock_delta_to_beta != -1) {
        close(sock_delta_to_beta);
        sock_delta_to_beta = -1;
    }
    if (sock_delta_to_alpha != -1) {
        close(sock_delta_to_alpha);
        sock_delta_to_alpha = -1;
    }

    PrintLog(LOG_LEVEL_SPECIAL, __FILE__, __LINE__, "Finalized Server Delta");

    return ret;
}

static int PerEpochOperations_delta(){
    int ret = 0;
    size_t received_sz = 0;

    PrintLog(LOG_LEVEL_INFO, __FILE__, __LINE__, "Server Delta: Starting PerEpochOperations sequence");

    /* Wait for receiving the ready message from server-beta */
    ret = recvAll(sock_delta_to_beta, net_buf, sizeof(net_buf), &received_sz);
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
    ret = recvAll(sock_delta_to_beta, net_buf, sizeof(net_buf), &received_sz);

    /* 4.d.2  Skipping the reception of mask. We are manually transferring them in chunks. */
    
    if (ret != 0)
    {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Failed to receive COMPLETED_REINIT_FOR_EPOCH message from Server Beta");
        return -1;
    }

    if (std::string(net_buf, received_sz) != completed_reinit_for_epoch_message) {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Did not receive expected COMPLETED_REINIT_FOR_EPOCH message from Server Beta");
        return -1;
    } else {
        PrintLog(LOG_LEVEL_INFO, __FILE__, __LINE__, "Server Delta: Completed re-initialization for new epoch, now ready to process client-requests..!!");
    }

    return ret;
}

static int ObliviouslySearchShelter_delta() {
#if 0/* TODO: Implement later */
    // Set up variables
    Fss fServer;
    ServerKeyEq K_alpha;
    int ret = 0;
    size_t received_sz = 0;
    int ret_recv;
    size_t dserializedFssSize;
    Ciphertext<DCRTPoly> fnd_alpha_ct_element, fnd_gamma_ct_element;
    Ciphertext<DCRTPoly> fnd_alpha_ct_tag, fnd_gamma_ct_tag;
    Ciphertext<DCRTPoly> random_ct;
    mpz_class d_ct_alpha = 0, random_pt, tmp_pt;
    mpz_class d_ct_gamma;
    std::vector<bool> thread_fnd(NUM_CPU_CORES, false);
    bool fnd_alpha = false;
    std::vector<mpz_class> thread_sums(NUM_CPU_CORES);
    // First, receive FSS parameters from the server Beta
    ret = recvAll(sock_delta_to_beta, net_buf, sizeof(net_buf), &received_sz);
    if (ret != 0)
    {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Failed to receive FSS parameters from Server Beta");
        return -1;
    }

    dserializedFssSize = deserializeFssAndServerKeyEq(net_buf, received_sz, fServer, K_alpha);

    /* For the verification purpose, set a particular location with special tag printed from server beta, to make the DPF search successful */
#if TEST_SHELTER_FOUND
    mpz_class special_tag, special_tag_location;

    PrintLog(LOG_LEVEL_SPECIAL, __FILE__, __LINE__, "Enter the value of the set search tag (base 10): ");
    mpz_inp_str(special_tag.get_mpz_t(), stdin, 10);

    PrintLog(LOG_LEVEL_SPECIAL, __FILE__, __LINE__, "Enter the index within the shelter, where this special tag must be placed (set the same value in server_gamma as well): ");
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
                if (evaluateEq(&fServer, &K_alpha, sh[k + j].tag_short)) {
                    mpz_xor(thread_sums[j].get_mpz_t(), thread_sums[j].get_mpz_t(), sh[k+j].element.get_mpz_t());

                    /* Same as XORing */
                    thread_fnd[j] = !thread_fnd[j];
                }
                else{
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
    ret_recv = recvAll(sock_delta_to_alpha, net_buf, sizeof(net_buf), &received_sz);
    if (ret_recv != 0)
    {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Failed to receive fnd_gamma_ct_element from Server Gamma");
        return -1;
    }
    Serial::DeserializeFromString(fnd_gamma_ct_element, std::string(net_buf, received_sz));

    // 7.3.2 Receive fnd_gamma_ct_tag
    ret_recv = recvAll(sock_delta_to_alpha, net_buf, sizeof(net_buf), &received_sz);
    if (ret_recv != 0)
    {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Failed to receive fnd_gamma_ct_tag from Server Gamma");
        return -1;
    }
    Serial::DeserializeFromString(fnd_gamma_ct_tag, std::string(net_buf, received_sz));

    // 7.4.2 Receive d_ct_gamma
    ret_recv = recvAll(sock_delta_to_alpha, net_buf, sizeof(net_buf), &received_sz);
    if (ret_recv != 0)
    {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Failed to receive d_ct_gamma from Server Gamma");
        return -1;
    }
    d_ct_gamma = mpz_class(std::string(net_buf, received_sz));

    /* Step 8.1 Determine fnd_alpha_ct */
    if (fnd_alpha == true){
        FHE_EncOfOnes(fnd_alpha_ct_element, fnd_alpha_ct_tag);
    }else{
        FHE_EncOfZeros(fnd_alpha_ct_element, fnd_alpha_ct_tag);
    }

    /* Step 8.2 Homomorphically compute fnd_ct = fnd_alpha_ct XOR fnd_alpha_ct = (fnd_alpha_ct + fnd_alpha_ct) - 2*(fnd_alpha_ct*fnd_alpha_ct) */
    fnd_ct_element = FHE_bitwise_XOR(fnd_alpha_ct_element, fnd_gamma_ct_element);
    fnd_ct_tag = FHE_bitwise_XOR(fnd_alpha_ct_tag, fnd_gamma_ct_tag);

    /* Step 8.2.1 Compute d_ct in mpz_class */
    mpz_xor(SR_sh_ct_mpz.get_mpz_t(), d_ct_gamma.get_mpz_t(), d_ct_alpha.get_mpz_t());

    /* Hack */
    #if 1
    if (SR_sh_ct_mpz == 0){
        if (Serial::SerializeToFile("/dev/shm/dummy_element.ct", vectorZeroesforElement_ct, SerType::BINARY) == true){
            SR_sh_ct_mpz = import_from_file_to_mpz_class("/dev/shm/dummy_element.ct");        
        }else{
            PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Failed to serialize the vectorZeroesforElement_ct ciphertext");
        }
    }
    #endif

    /* Setp 8.3 Convert from mpz_class to FHE ciphertext */
    export_to_file_from_mpz_class("/dev/shm/d.ct", SR_sh_ct_mpz);

    if (!Serial::DeserializeFromFile("/dev/shm/d.ct", SR_sh_ct, SerType::BINARY)) {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Cannot convert to the d_ct to Ciphertext<DCRTPoly>");
    }
    
    /****************** Refresh fnd_ct_element ******************/
    /* Create a random SDBElement having PLAINTEXT_PIR_BLOCK_DATA_SIZE-bit data and log_N bit index */
    random_pt = rng.get_z_bits((PLAINTEXT_PIR_BLOCK_DATA_SIZE +  log_N));
    random_ct = FHE_Enc_SDBElement(random_pt);

    /* Mask the actual ciphertext using the random */
    fnd_ct_element = (fnd_ct_element + random_ct);

    /* Send to server beta for a refresh operation */
    (void)sendAll(sock_delta_to_beta, Serial::SerializeToString(fnd_ct_element).c_str(), Serial::SerializeToString(fnd_ct_element).size());

    /* Receive the refreshed ciphertext */
    ret_recv = recvAll(sock_delta_to_beta, net_buf, sizeof(net_buf), &received_sz);
    if (ret_recv != 0)
    {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Failed to receive refreshed fnd_ct_element from Server Beta");
        return -1;
    }
    Serial::DeserializeFromString(fnd_ct_element, std::string(net_buf, received_sz));

    /* Remove the random to get back the usable ciphertext */
    fnd_ct_element = (fnd_ct_element - random_ct); 

    /****************** Refresh fnd_ct_tag ******************/
    /* Generate random tag and encrypt that */
    random_pt = rng.get_z_bits(P_BITS);
    random_ct = FHE_Enc_Tag(random_pt);
    
    /* Homomorphically mask the ciphertext  */
    fnd_ct_tag = (fnd_ct_tag + random_ct);
    (void)sendAll(sock_delta_to_beta, Serial::SerializeToString(fnd_ct_tag).c_str(), Serial::SerializeToString(fnd_ct_tag).size());
    
    /* Receive refreshed fnd_ct_tag */
    ret_recv = recvAll(sock_delta_to_beta, net_buf, sizeof(net_buf), &received_sz);
    if (ret_recv != 0)
    {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Failed to receive refreshed fnd_ct_tag from Server Beta");
        return -1;
    }
    Serial::DeserializeFromString(fnd_ct_tag, std::string(net_buf, received_sz));

    /* Remove the random to get back the usable ciphertext */
    fnd_ct_tag = (fnd_ct_tag - random_ct); 
#endif
    return 0;
}

static int ProcessClientRequest_delta(){
    int ret = -1;

    //struct sockaddr_in address;
    //int opt = 1;
    //int addrlen = sizeof(address);
    //int accepted_socket = -1;
    uint64_t M = (N + sqrt_N);
    item_type Kuku_key;    
    QueryResult res;
    std::fstream L;
    std::fstream DK;
    std::ifstream importedHFile;

    PrintLog(LOG_LEVEL_INFO, __FILE__, __LINE__, "Server Delta: Starting Request processing sequence");

    /* First of all retrieve all the one-time initialized materials from the saved location */
    Serial::DeserializeFromFile(ONE_TIME_MATERIALS_LOCATION_DELTA + "FHEcryptoContext.bin", FHEcryptoContext, SerType::BINARY);
    Serial::DeserializeFromFile(ONE_TIME_MATERIALS_LOCATION_DELTA + "pk_F.bin", pk_F, SerType::BINARY);

    PrintLog(LOG_LEVEL_TRACE, __FILE__, __LINE__, "Server Delha: Loaded one-time initialization materials");

    //Always initialize them
    K = 0;
#if 0/* TODO: Implement later */
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

    /* Open the shuffled and secret shared database */
    sdb.open(mdb_filename, std::ios::in | std::ios::binary);
    if (!sdb) {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Failed to open SDB file at location: " + mdb_filename);
        ret = -1;
        goto exit;
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

        /* Even if for the request number 1, the determined shelter tag is not required to be used.
           Still, we cannot move this function within if (K > 0) block, since the client interaction is involved in this function. */
        ret = ShelterTagDetermination_alpha();
        if (ret != 0)
        {
            PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Problem while determining the shelter tag..!!");
            ret = -1;
            goto exit;
        }

        /* For the first request, the shelter is not required to be searched */
        if (K > 0){
            ret = ObliviouslySearchShelter_delta();
            if (ret != 0)
            {
                PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Problem during the shelter search operation..!!");
                ret = -1;
                goto exit;
            }
        } else {
            /* For the first request set the fnd_ct to encryption of Zeros */
            FHE_EncOfZeros(fnd_ct_element, fnd_ct_tag);
            SR_sh_ct = vectorZeroesforElement_ct;
        }

        ret = SelShuffDBSearchTag_alpha();

        if (ret != 0){
            PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Problem during the selecting shuffled database tag..!!");
            ret = -1;
            goto exit;
        }

        ret = FetchCombineSelect_alpha();
        if (ret != 0){
            PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Problem during the Fetch_Combine_and_Select stage..!!");
            ret = -1;
            goto exit;
        }

        ret = ShelterUpdate_alpha();
        if (ret != 0){
            PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Problem during the Shelter Update stage..!!");
            ret = -1;
            goto exit;
        }

        /* Close the connection with existing client */
        close(sock_alpha_client_srv);
        close(sock_alpha_client_con);  
        K++;
    }

    PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Current epoch is completed. Please re-perform the per-epoch initialization");


exit:
    sdb.close();
    /* If there are any dangling connection with client, close that */
    close(sock_alpha_client_srv);
    close(sock_alpha_client_con);    
#endif

    return ret;
}


static void TestSrv_delta()
{
    //TestPKEOperations_alpha();
    //TestSelShuffDBSearchTag_alpha();
    //TestShelterDPFSearch_alpha();
    //TestClientProcessing_alpha();
    //TestHTableSerDser_alpha();

    return;
}

ostream &operator<<(ostream &stream, item_type item)
{
    stream << item[1] << " " << item[0];
    return stream;
}

/*******************************************************************************************************
    On average, the shelter will be half-full. So its performance can be measured by observing the
    required time for (\sqrt{N}/2)st request processing. Or this function can be called, which 
    simulates the same situation by pre-populating the shelter with (\sqrt{N}/2)-random elements
    and then measure the performance.
********************************************************************************************************/
static int Perf_avg_online_server_time_delta() {
    #if 0/* TODO: Implement later */
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
    std::chrono::high_resolution_clock::time_point t0;
    double processingTime_us = 0.0;    


    /* First of all retrieve all the one-time initialized materials from the saved location */
    p = import_from_file_to_mpz_class(ONE_TIME_MATERIALS_LOCATION_ALPHA + "p.bin");
    q = import_from_file_to_mpz_class(ONE_TIME_MATERIALS_LOCATION_ALPHA + "q.bin");
    g = import_from_file_to_mpz_class(ONE_TIME_MATERIALS_LOCATION_ALPHA + "g.bin");
    r = import_from_file_to_mpz_class(ONE_TIME_MATERIALS_LOCATION_ALPHA + "r.bin");
    Serial::DeserializeFromFile(ONE_TIME_MATERIALS_LOCATION_ALPHA + "FHEcryptoContext.bin", FHEcryptoContext, SerType::BINARY);
    Serial::DeserializeFromFile(ONE_TIME_MATERIALS_LOCATION_ALPHA + "pk_F.bin", pk_F, SerType::BINARY);


    PrintLog(LOG_LEVEL_TRACE, __FILE__, __LINE__, "Starting to randomly populate a shelter");

    /* Populate the shelter, with random elements */
    for(size_t k = 0; k < average_shelter_size; k++) {
        // Generate random block_content of PLAINTEXT_PIR_BLOCK_DATA_SIZE bits of random | k as the block index
        Ciphertext<DCRTPoly> tmp_ct = FHE_Enc_SDBElement((rng.get_z_bits(PLAINTEXT_PIR_BLOCK_DATA_SIZE) << log_N) | mpz_class(k));
        /* Store the ciphertexts to serialized form to a file, which resides in the RAM */
        if (Serial::SerializeToFile(DPF_search_test_shelter_location + "sh[" + std::to_string(k) + "].ct", tmp_ct, SerType::BINARY) != true)
        {
            PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Failed to serialize element FHE ciphertext to file");
        }
        /* Just notedown the ciphertext size */
        if (k == 0){
            PrintLog(LOG_LEVEL_TRACE, __FILE__, __LINE__, "Size of the FHE-ciphertext: "+ std::to_string(std::filesystem::file_size(DPF_search_test_shelter_location + "sh[" + std::to_string(k) + "].ct")));
        }


        sh[k].element_FHE_ct = import_from_file_to_mpz_class(DPF_search_test_shelter_location + "sh[" + std::to_string(k) + "].ct");

        /* Generate the tags and keep them in the variable, which will be used for DPF search */
        sh[k].tag = ElGamal_randomGroupElement(); // Create a random tag
        sh[k].tag_short = sh[k].tag % r; // Create a random short tag
    }

    PrintLog(LOG_LEVEL_TRACE, __FILE__, __LINE__, "Starting to test DPF-search on the shelter of size: "  + to_string(average_shelter_size));

    t0 = std::chrono::high_resolution_clock::now();

    /*******************************************************************************************************
        In this experiment, only considering the required time for two expensive operations.
        Executing DPF search on large number of shelter elements and updating all the existing shelter tags.
        In comparison, other operations takes very less amount of time and does not depend on N.
        Approximate the time required for those operations from other experiments.
    ********************************************************************************************************/

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

    {
        auto t1 = std::chrono::high_resolution_clock::now();
        auto elapsed_us = std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0).count();
        processingTime_us = static_cast<double>(elapsed_us); // microseconds
    }

    /* Consider the additional time required for other non-intensive tasks. Those can be found from other experiments */

    std::cout << "Processing time is: " << processingTime_us << "us" << std::endl;
    PrintLog(LOG_LEVEL_SPECIAL, __FILE__, __LINE__, "Online time consumption by the server in the average scenario is: " + std::to_string(processingTime_us) + "us");
#endif

    return 1;
}

int main(int argc, char *argv[])
{
    int ret = -1;

    /* Perform the basic initialization */
    InitSrv_delta();

    /* Process as per the command line arguments */
    if (argc >= 2) {
        if (std::string("one_time_init").compare(std::string(argv[1]))==0) {
            // Perform one-time initialization for server alpha
            ret = OneTimeInit_delta();
        } else if (std::string("per_epoch_operations").compare(std::string(argv[1]))==0) {
            // Perform per-epoch initialization for server alpha
            ret = PerEpochOperations_delta();
        } else if (std::string("clear_epoch_state").compare(std::string(argv[1]))==0) {
            // Clear the existing state of current epoch, start as if this is the first request of the epoch
            // Delete shelter content and set K = 0
        } else if (std::string("process_request").compare(std::string(argv[1]))==0) {
            // Start from last saved state
            ret = ProcessClientRequest_delta();
        } else if (std::string("test").compare(std::string(argv[1]))==0) {
            TestSrv_delta();
        } else if (std::string("perf").compare(std::string(argv[1]))==0) {
            if (argc < 3){
                PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Performance measurement option requires at least three command line parameters. Usage: server_alpha perf [srv_avg_online_time]");
            }else{
                if (std::string("srv_avg_online_time").compare(std::string(argv[2]))==0){
                    (void)Perf_avg_online_server_time_delta();
                }
            }
        } else {
            PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Unknown command line argument:"+ std::string(argv[1]));
            PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Improper command line arguments. Usage: server_alpha <one_time_init|per_epoch_operations|clear_epoch_state|process_request|perf>");
        }
    } else {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Improper command line arguments. Usage: server_alpha <one_time_init|per_epoch_operations|clear_epoch_state|process_request|perf [srv_avg_online_time]>");
    }

    if (ret == 0) {
        //TestSrv_delta();
    }

    FinSrv_delta();

    return 0;
}