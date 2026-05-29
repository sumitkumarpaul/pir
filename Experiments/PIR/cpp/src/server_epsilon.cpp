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

static int sock_epsilon_to_beta = -1, sock_epsilon_to_alpha = -1, sock_epsilon_to_gamma = -1;
static char net_buf[NET_BUF_SZ] = {0};
static mpz_class M[sqrt_N]; /* This arrary extracts entire mask database into a RAM array. Note this is in mpz_class format. */
static char y_gamma_bits_buf[(sqrt_N+7)/8];
static std::fstream mdb;

static uint64_t K; // Current number of entries in the shelter, or the number of processed requests

#define NUM_CPU_CORES 16
#define ONE_TIME_MATERIALS_LOCATION_EPSILON std::string("/mnt/sumit/PIR_EPSILON/ONE_TIME_MATERIALS/")
#define PER_EPOCH_MATERIALS_LOCATION_EPSILON std::string("/mnt/sumit/PIR_EPSILON/PER_EPOCH_MATERIALS/")
#define MASK_LOCATION_EPSILON std::string("/mnt/sumit/PIR_EPSILON/")
std::string mdb_filename = PER_EPOCH_MATERIALS_LOCATION_EPSILON+"MaskDB.bin";
#define TMP_FILE std::string("/dev/shm/tmp_epsilon")


// Function declarations
static int InitSrv_epsilon();
static int OneTimeInit_epsilon();
static int FinSrv_epsilon();
static int PerEpochOperations_epsilon();
static int ProcessClientRequest_epsilon();
static int ObliviouslySearchShelter_epsilon();
static void TestSrv_epsilon();

static int Perf_avg_online_server_time_epsilon();

// Function definitions
static int InitSrv_epsilon(){
    int ret = -1;
    // Initialize random number generation
    std::random_device rd;
    unsigned long seed = (static_cast<unsigned long>(rd()) << 1) ^ rd();
    rng.seed(seed); // seed() seeds the gmp_randclass    
    
    // Server_alpha only connects to other servers, it does not listen to other servers
    InitConnectingSocket(SERVER_BETA_IP, BETA_LISTENING_TO_EPSILON_PORT, &sock_epsilon_to_beta);

    PrintLog(LOG_LEVEL_TRACE, __FILE__, __LINE__, "Established connection with Server Beta");

    InitConnectingSocket(SERVER_ALPHA_IP, ALPHA_LISTENING_TO_EPSILON_PORT, &sock_epsilon_to_alpha);

    PrintLog(LOG_LEVEL_TRACE, __FILE__, __LINE__, "Established connection with Server Alpha");

    InitConnectingSocket(SERVER_GAMMA_IP, GAMMA_LISTENING_TO_EPSILON_PORT, &sock_epsilon_to_gamma);

    PrintLog(LOG_LEVEL_TRACE, __FILE__, __LINE__, "Established connection with Server Gamma");

    PrintLog(LOG_LEVEL_INFO, __FILE__, __LINE__, "Server Epsilon initialization complete");

    ret = 0;

exit:
    if (ret != 0){
        FinSrv_epsilon();
    }

    return ret;
}

static int OneTimeInit_epsilon() {
    size_t received_sz = 0;
    int ret_recv = 0;

    // Receive all the parameters from server beta
    // Receive FHEcryptoContext
    ret_recv = recvFile(sock_epsilon_to_beta, net_buf, sizeof(net_buf), TMP_FILE);
    if (ret_recv != 0)
    {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Failed to receive FHEcryptoContext from Server Beta");
        return -1;
    }
    Serial::DeserializeFromFile(TMP_FILE, FHEcryptoContext, SerType::BINARY);

    // Receive pk_F
    ret_recv = recvFile(sock_epsilon_to_beta, net_buf, sizeof(net_buf), TMP_FILE);
    if (ret_recv != 0)
    {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Failed to receive pk_F from Server Beta");
        return -1;
    }
    Serial::DeserializeFromFile(TMP_FILE, pk_F, SerType::BINARY);

    //Save parameters to local files
    if (!Serial::SerializeToFile(ONE_TIME_MATERIALS_LOCATION_EPSILON + "FHEcryptoContext.bin", FHEcryptoContext, SerType::BINARY)){
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Failed to serialize the received FHEcryptoContext from Server Beta");
        return -1;
    }
    
    if (!Serial::SerializeToFile(ONE_TIME_MATERIALS_LOCATION_EPSILON + "pk_F.bin", pk_F, SerType::BINARY)){
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Failed to serialize the received FHEcryptoContext from Server Beta");
        return -1;
    }

    PrintLog(LOG_LEVEL_TRACE, __FILE__, __LINE__, "Received all the one-time initialized parameters from Server Beta and exported all of them into file");
    
    return 0;
}

static int FinSrv_epsilon(){
    int ret = -1;

    // Close the sockets
    if (sock_epsilon_to_alpha != -1) {
        close(sock_epsilon_to_alpha);
        sock_epsilon_to_alpha = -1;
    }
    if (sock_epsilon_to_beta != -1) {
        close(sock_epsilon_to_beta);
        sock_epsilon_to_beta = -1;
    }
    if (sock_epsilon_to_gamma != -1) {
        close(sock_epsilon_to_gamma);
        sock_epsilon_to_gamma = -1;
    }

    PrintLog(LOG_LEVEL_INFO, __FILE__, __LINE__, "Finalized Server Epsilon");

    return ret;
}

static int PerEpochOperations_epsilon(){
    int ret = 0;
    size_t received_sz = 0;

    PrintLog(LOG_LEVEL_INFO, __FILE__, __LINE__, "Server Epsilon: Starting PerEpochOperations sequence");

    /* Wait for receiving the ready message from server-beta */
    ret = recvAll(sock_epsilon_to_beta, net_buf, sizeof(net_buf), &received_sz);
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
    ret = recvAll(sock_epsilon_to_beta, net_buf, sizeof(net_buf), &received_sz);

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
        PrintLog(LOG_LEVEL_INFO, __FILE__, __LINE__, "Server Epsilon: Completed re-initialization for new epoch, now ready to process client-requests..!!");
    }

    return ret;
}

static int ObliviouslySearchShelter_epsilon() {
    int ret = 0;
    size_t received_sz = 0;
    int ret_recv;
    Ciphertext<DCRTPoly> m_epsilon_ct;

    // 4.b Initialize with zeros
    mpz_class m_epsilon = 0;
    std::vector<bool> fnd_epsilon_thread(NUM_CPU_CORES, false);
    bool fnd_epsilon = false;
    std::vector<mpz_class> m_epsilon_thread(NUM_CPU_CORES);

    // 6.a.3 Receive the entire bit array from the server gamma
    ret = recvAll(sock_epsilon_to_gamma, y_gamma_bits_buf, sizeof(y_gamma_bits_buf), &received_sz);
    if (ret != 0)
    {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Failed to receive bit array from Server Gamma");
        return -1;
    }

    PrintLog(LOG_LEVEL_TRACE, __FILE__, __LINE__, "Received the array of bits from Server Gamma, the received size is: " + std::to_string(received_sz));

    if (SET_LOG_LEVEL >= LOG_LEVEL_TRACE)
    {
        printf("Received bits are: ");
        for (size_t k = 0; k < received_sz; k++)
        {
            printf("%02x", y_gamma_bits_buf[k]);
        }
        printf("\n");
    }    


    for (size_t k = 0; k < K; k += NUM_CPU_CORES)
    {
        for (int t = 0; t < NUM_CPU_CORES; ++t){
            m_epsilon_thread[t] = 0;
            fnd_epsilon_thread[t] = false;
        }

        #pragma omp parallel for
        for (int j = 0; j < NUM_CPU_CORES; ++j)
        {
            if ((k + j) < K)
            {
                if (y_gamma_bits_buf[(k + j) / 8] & (1 << ((k + j) % 8))) {
                    mpz_xor(m_epsilon_thread[j].get_mpz_t(), m_epsilon_thread[j].get_mpz_t(), M[k+j].get_mpz_t());

                    #pragma omp critical
                    {
                        fnd_epsilon_thread[j] = fnd_epsilon_thread[j]^true;
                    }
                }
            }
        }
        for (int t = 0; t < NUM_CPU_CORES; ++t)
        {
            mpz_xor(m_epsilon.get_mpz_t(), m_epsilon.get_mpz_t(), m_epsilon_thread[t].get_mpz_t());
            fnd_epsilon ^= fnd_epsilon_thread[t];
        }
    }

    PrintLog(LOG_LEVEL_TRACE, __FILE__, __LINE__, "Completed processing the mask database. Value of fnd_epsilon: " + std::to_string(fnd_epsilon));
    PrintLog(LOG_LEVEL_TRACE, __FILE__, __LINE__, "Value of the S_epsilon's share of the mask is: " + m_epsilon.get_str(16));

    // 12.1 Compute the FHE ciphertext of m_epsilon
    m_epsilon_ct = FHE_bitwise_Enc_SDBElement(m_epsilon);

    /* 12.2 Send the ciphertext to server alpha */
    Serial::SerializeToFile(TMP_FILE, m_epsilon_ct, SerType::BINARY);
    (void)sendFile(sock_epsilon_to_alpha, net_buf, sizeof(net_buf), TMP_FILE);    

    return 0;
}

static int ProcessClientRequest_epsilon(){
    int ret = -1;
    shuffled_db_entry tmp;

    PrintLog(LOG_LEVEL_INFO, __FILE__, __LINE__, "Server Epsilon: Starting Request processing sequence");

    /* First of all retrieve all the one-time initialized materials from the saved location */
    Serial::DeserializeFromFile(ONE_TIME_MATERIALS_LOCATION_EPSILON + "FHEcryptoContext.bin", FHEcryptoContext, SerType::BINARY);
    Serial::DeserializeFromFile(ONE_TIME_MATERIALS_LOCATION_EPSILON + "pk_F.bin", pk_F, SerType::BINARY);
    /* Load the mask database into the RAM location for faster access */
    mdb.open(mdb_filename, std::ios::in | std::ios::binary);

    for (uint64_t iter = 0; iter < sqrt_N; iter++) {
        read_mdb_entry(mdb, iter, tmp);
        mpz_import(M[iter].get_mpz_t(), NUM_BYTES_PER_SDB_ELEMENT, 1, 1, 1, 0, tmp.element);
    }

    PrintLog(LOG_LEVEL_TRACE, __FILE__, __LINE__, "Server Epsilon: Loaded one-time initialization materials");

    //Always initialize them
    K = 0;
    while (K < sqrt_N){
        /* For the first request, the shelter is not required to be searched */
        if (K > 0){
            ret = ObliviouslySearchShelter_epsilon();
            if (ret != 0)
            {
                PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Problem during the shelter search operation..!!");
                ret = -1;
                goto exit;
            }
        }
        K++;
    }

    PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Current epoch is completed. Please re-perform the per-epoch initialization");

exit:
    mdb.close();
    /* Close any dangling connections */
    close(sock_epsilon_to_alpha);
    close(sock_epsilon_to_beta);
    close(sock_epsilon_to_gamma);

    return ret;
}


static void TestSrv_epsilon()
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
static int Perf_avg_online_server_time_epsilon() {
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
    /* Generate a dummy epsilon value to update the shelter tags */
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
        Ciphertext<DCRTPoly> tmp_ct = FHE_bitwise_Enc_SDBElement((rng.get_z_bits(PLAINTEXT_PIR_BLOCK_DATA_SIZE) << log_N) | mpz_class(k));
        /* Store the ciphertexts to serialized form to a file, which resides in the RAM */
        if (Serial::SerializeToFile(DPF_search_test_shelter_location + "sh[" + std::to_string(k) + "].ct", tmp_ct, SerType::BINARY) != true)
        {
            PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Failed to serialize element FHE ciphertext to file");
        }
        /* Just notedown the ciphertext size */
        if (k == 0){
            PrintLog(LOG_LEVEL_TRACE, __FILE__, __LINE__, "Size of the FHE-ciphertext: "+ std::to_string(std::filesystem::file_size(DPF_search_test_shelter_location + "sh[" + std::to_string(k) + "].ct")));
        }


        sh[k].element = import_from_file_to_mpz_class(DPF_search_test_shelter_location + "sh[" + std::to_string(k) + "].ct");

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
                    mpz_xor(thread_sums[j].get_mpz_t(), thread_sums[j].get_mpz_t(), sh[k+j].element.get_mpz_t());
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
    PrintLog(LOG_LEVEL_INFO, __FILE__, __LINE__, "Online time consumption by the server in the average scenario is: " + std::to_string(processingTime_us) + "us");
#endif

    return 1;
}

int main(int argc, char *argv[])
{
    int ret = -1;

    /* Perform the basic initialization */
    InitSrv_epsilon();

    /* Process as per the command line arguments */
    if (argc >= 2) {
        if (std::string("one_time_init").compare(std::string(argv[1]))==0) {
            // Perform one-time initialization for server alpha
            ret = OneTimeInit_epsilon();
        } else if (std::string("per_epoch_operations").compare(std::string(argv[1]))==0) {
            // Perform per-epoch initialization for server alpha
            ret = PerEpochOperations_epsilon();
        } else if (std::string("clear_epoch_state").compare(std::string(argv[1]))==0) {
            // Clear the existing state of current epoch, start as if this is the first request of the epoch
            // Delete shelter content and set K = 0
        } else if (std::string("process_request").compare(std::string(argv[1]))==0) {
            // Start from last saved state
            ret = ProcessClientRequest_epsilon();
        } else if (std::string("test").compare(std::string(argv[1]))==0) {
            TestSrv_epsilon();
        } else if (std::string("perf").compare(std::string(argv[1]))==0) {
            if (argc < 3){
                PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Performance measurement option requires at least three command line parameters. Usage: server_alpha perf [srv_avg_online_time]");
            }else{
                if (std::string("srv_avg_online_time").compare(std::string(argv[2]))==0){
                    (void)Perf_avg_online_server_time_epsilon();
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
        //TestSrv_epsilon();
    }

    FinSrv_epsilon();

    return 0;
}