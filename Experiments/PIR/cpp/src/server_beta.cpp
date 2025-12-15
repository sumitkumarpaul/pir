#include <iostream>
#include <vector>
#include <random>
#include <cmath>    // std::sqrt, std::ceil
#include <algorithm> // std::swap
#include "pir_common.h"

//Global variables
static char net_buf[NET_BUF_SZ] = {0};
static int sock_beta_alpha_srv = -1, sock_beta_alpha_con = -1;
static int sock_beta_gamma_srv = -1, sock_beta_gamma_con = -1;
static int sock_beta_client_srv = -1, sock_beta_client_con = -1;
static std::vector<mpz_class> SetPhi;
static std::fstream pdb;
static std::fstream D_K;
static std::fstream D_alpha;
static std::fstream D_gamma;
static uint64_t K;/* The number of requests processed till now */
static mpz_class b;
static mpz_class widehat_t_I;

static shelter_element sh[sqrt_N]; /* TODO For debugging only */
#define NUM_CPU_CORES 16 /* TODO For debugging only */
//#define SHELTER_STORING_LOCATION std::string("/mnt/sumit/dummy_shelter/") /* TODO For debugging only */
#define SHELTER_STORING_LOCATION std::string("/dev/shm/") /* TODO For debugging only */

#define ONE_TIME_MATERIALS_LOCATION_BETA std::string("/mnt/sumit/PIR_BETA/ONE_TIME_MATERIALS/")
#define PER_EPOCH_MATERIALS_LOCATION_BETA std::string("/mnt/sumit/PIR_BETA/PER_EPOCH_MATERIALS/")
#define DATABASE_LOCATION_BETA std::string("/mnt/sumit/PIR_BETA/")
std::string pdb_filename = DATABASE_LOCATION_BETA+"PlaintextDB.bin";
std::string DK_filename = PER_EPOCH_MATERIALS_LOCATION_BETA+"DK.bin";
std::string D_alpha_filename = PER_EPOCH_MATERIALS_LOCATION_BETA+"D_alpha.bin";
std::string D_gamma_filename = PER_EPOCH_MATERIALS_LOCATION_BETA+"D_gamma.bin";
std::string SetPhi_filename = PER_EPOCH_MATERIALS_LOCATION_BETA+"SetPhi.bin";

#define NUM_ITEMS_IN_TMP_BUF 16//((sqrt_N) * (NUM_CPU_CORES)) //100000//16//4112//6553760//100000 // After these many item creation, everything is written to the disk. This must be multiple of number of threads used in parallel for loop and also must be a divisor of (N+sqrt_N)
//static unsigned char TMP_KEY_BUF[((sizeof(item_type))* NUM_ITEMS_IN_TMP_BUF)] = {0};
static std::array<unsigned char, 16> TMP_KEY_BUF[NUM_ITEMS_IN_TMP_BUF] = {0};
static unsigned char TMP_D_ALPHA_BUF[(NUM_BYTES_PER_SDB_ELEMENT * NUM_ITEMS_IN_TMP_BUF)] = {0};
static unsigned char TMP_D_GAMMA_BUF[(NUM_BYTES_PER_SDB_ELEMENT * NUM_ITEMS_IN_TMP_BUF)] = {0};

#if TEST_SHUFF_DB_FETCH
static uint64_t TMP_IDX_LOC_MAP[(N+sqrt_N)] = {0};// TODO: Delete this array
#endif

// Function declarations
static void Init_parameters(int p_bits = 3072, int q_bits = 256, int r_bits = 64);// Initializes p, q, g, GG(cyclic group) and r
static int InitSrv_beta();
static int FinSrv_beta();
static int OneTimeInit_beta();
static int SendInitializedParamsToAllServers();
static int SelShuffDBSearchTag_beta();
static int PerEpochOperations_beta();
static int CreateRandomDatabase();
static int ProcessClientRequest_beta();
static int ShelterTagDetermination_beta();
static int ObliviouslySearchShelter_beta();
static int ShelterUpdate_beta();

static void TestSrv_beta();

static void TestPKEOperations_beta();
static void TestBlindedExponentiation();
static void TestBlindedExponentiation1();
static void TestBlindedExponentiation2();
static void Test_FHE_DBElement();
static void TestSelShuffDBSearchTag_beta();
static int TestShelterDPFSearch_beta();
static int TestClientProcessing_beta();
#if TEST_SHUFF_DB_FETCH
static void TestShuffDBFetch_beta();
#endif


static void Perf_avg_online_server_time_beta();

static void Init_parameters(int p_bits, int q_bits, int r_bits) {
    mpz_class sg_prime;

    //Choose a safe prime q
    do{
        // First randomly choose a (q_bits - 2)-bits long Sophie Germain primes
        do
        {
            sg_prime = rng.get_z_bits(q_bits - 2);
            mpz_nextprime(sg_prime.get_mpz_t(), sg_prime.get_mpz_t());
        } while (mpz_sizeinbase(sg_prime.get_mpz_t(), 2) != (q_bits - 2));

        //And check whether this generates a safe-prime or not
        q = (2*sg_prime) + 1;
    } while (!mpz_probab_prime_p(q.get_mpz_t(), 25));
    
    //Accordingly choose p
    mpz_class temp;
    do {
        temp = rng.get_z_bits(p_bits - q_bits);
        p = temp * q + 1;
    } while (!mpz_probab_prime_p(p.get_mpz_t(), 25));

    //Then determine the generator g, which will define the cyclic group GG
    mpz_class h, exp;
    exp = (p - 1) / q;
    do {
        h = rng.get_z_range(p - 1) + 1;
        mpz_powm(g.get_mpz_t(), h.get_mpz_t(), exp.get_mpz_t(), p.get_mpz_t());
    } while (g == 1);

    //Choose 2 as the generator of the multiplicative sub-group ZZ_q*
    g_q = mpz_class(2);
 
    //Randomly choose a prime number r
    do {
        r = rng.get_z_bits(r_bits);
        mpz_nextprime(r.get_mpz_t(), r.get_mpz_t());
    } while (mpz_sizeinbase(r.get_mpz_t(), 2) != r_bits);

    PrintLog(LOG_LEVEL_TRACE, __FILE__, __LINE__, "Initialized p,q,g and r");

    return;
}

static int SendInitializedParamsToAllServers(){
    int ret = 0;

    //Save parameters to local files
    export_to_file_from_mpz_class(ONE_TIME_MATERIALS_LOCATION_BETA + "p.bin", p);
    export_to_file_from_mpz_class(ONE_TIME_MATERIALS_LOCATION_BETA + "q.bin", q);
    export_to_file_from_mpz_class(ONE_TIME_MATERIALS_LOCATION_BETA + "g.bin", g);
    export_to_file_from_mpz_class(ONE_TIME_MATERIALS_LOCATION_BETA + "g_q.bin", g_q);
    export_to_file_from_mpz_class(ONE_TIME_MATERIALS_LOCATION_BETA + "r.bin", r);
    export_to_file_from_mpz_class(ONE_TIME_MATERIALS_LOCATION_BETA + "pk_E.bin", pk_E);
    export_to_file_from_mpz_class(ONE_TIME_MATERIALS_LOCATION_BETA + "sk_E.bin", sk_E);
    export_to_file_from_mpz_class(ONE_TIME_MATERIALS_LOCATION_BETA + "pk_E_q.bin", pk_E_q);
    export_to_file_from_mpz_class(ONE_TIME_MATERIALS_LOCATION_BETA + "sk_E_q.bin", sk_E_q);

    Serial::SerializeToFile(ONE_TIME_MATERIALS_LOCATION_BETA + "FHEcryptoContext.bin", FHEcryptoContext, SerType::BINARY);
    Serial::SerializeToFile(ONE_TIME_MATERIALS_LOCATION_BETA + "pk_F.bin", pk_F, SerType::BINARY);
    Serial::SerializeToFile(ONE_TIME_MATERIALS_LOCATION_BETA + "sk_F.bin", sk_F, SerType::BINARY);
    Serial::SerializeToFile(ONE_TIME_MATERIALS_LOCATION_BETA + "vectorOnesforElement_ct.bin", vectorOnesforElement_ct, SerType::BINARY);
    Serial::SerializeToFile(ONE_TIME_MATERIALS_LOCATION_BETA + "vectorOnesforTag_ct.bin", vectorOnesforTag_ct, SerType::BINARY);

    PrintLog(LOG_LEVEL_TRACE, __FILE__, __LINE__, "Sending initialized parameters to server alpha");

    //Send parameters to Server Alpha
    (void)sendAll(sock_beta_alpha_con, p.get_str().c_str(), p.get_str().size());
    (void)sendAll(sock_beta_alpha_con, q.get_str().c_str(), q.get_str().size());
    (void)sendAll(sock_beta_alpha_con, g.get_str().c_str(), g.get_str().size());
    (void)sendAll(sock_beta_alpha_con, g_q.get_str().c_str(), g_q.get_str().size());
    (void)sendAll(sock_beta_alpha_con, r.get_str().c_str(), r.get_str().size());
    (void)sendAll(sock_beta_alpha_con, pk_E.get_str().c_str(), pk_E.get_str().size());
    (void)sendAll(sock_beta_alpha_con, pk_E_q.get_str().c_str(), pk_E_q.get_str().size());
    (void)sendAll(sock_beta_alpha_con, Serial::SerializeToString(FHEcryptoContext).c_str(), Serial::SerializeToString(FHEcryptoContext).size());
    (void)sendAll(sock_beta_alpha_con, Serial::SerializeToString(pk_F).c_str(), Serial::SerializeToString(pk_F).size());
    (void)sendAll(sock_beta_alpha_con, Serial::SerializeToString(vectorOnesforElement_ct).c_str(), Serial::SerializeToString(vectorOnesforElement_ct).size());
    (void)sendAll(sock_beta_alpha_con, Serial::SerializeToString(vectorOnesforTag_ct).c_str(), Serial::SerializeToString(vectorOnesforTag_ct).size());

    PrintLog(LOG_LEVEL_TRACE, __FILE__, __LINE__, "Sending initialized parameters to server gamma");

    //Send parameters to Server Gamma
    (void)sendAll(sock_beta_gamma_con, p.get_str().c_str(), p.get_str().size());
    (void)sendAll(sock_beta_gamma_con, q.get_str().c_str(), q.get_str().size());
    (void)sendAll(sock_beta_gamma_con, g.get_str().c_str(), g.get_str().size());
    (void)sendAll(sock_beta_gamma_con, g_q.get_str().c_str(), g_q.get_str().size());
    (void)sendAll(sock_beta_gamma_con, r.get_str().c_str(), r.get_str().size());
    (void)sendAll(sock_beta_gamma_con, pk_E.get_str().c_str(), pk_E.get_str().size());
    (void)sendAll(sock_beta_gamma_con, pk_E_q.get_str().c_str(), pk_E_q.get_str().size());
    (void)sendAll(sock_beta_gamma_con, Serial::SerializeToString(FHEcryptoContext).c_str(), Serial::SerializeToString(FHEcryptoContext).size());
    (void)sendAll(sock_beta_gamma_con, Serial::SerializeToString(pk_F).c_str(), Serial::SerializeToString(pk_F).size());
    (void)sendAll(sock_beta_gamma_con, Serial::SerializeToString(vectorOnesforElement_ct).c_str(), Serial::SerializeToString(vectorOnesforElement_ct).size());
    (void)sendAll(sock_beta_gamma_con, Serial::SerializeToString(vectorOnesforTag_ct).c_str(), Serial::SerializeToString(vectorOnesforTag_ct).size());

    return 0;
}

static int OneTimeInit_beta(){
    int ret = 0;

    //Initialize p, q, g, GG(cyclic group) and r
    Init_parameters(P_BITS, Q_BITS, R_BITS);

    //Initialize El-Gamal key-pair
    std::tie(pk_E, sk_E) = ElGamal_keyGen();
    PrintLog(LOG_LEVEL_TRACE, __FILE__, __LINE__, "Generated El-Gamal key-pair");

    //Initialize El-Gamal_q key-pair
    std::tie(pk_E_q, sk_E_q) = ElGamal_q_keyGen();
    PrintLog(LOG_LEVEL_TRACE, __FILE__, __LINE__, "Generated El-Gamal key-pair in ZZ*_q");

    //Initialize FHE key-pair
    ret = FHE_keyGen();
    
    if (ret != 0) {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Failed to generate FHE key-pair");
        return -1;
    } else {
        PrintLog(LOG_LEVEL_TRACE, __FILE__, __LINE__, "Generated FHE key-pair");
    }

    FHE_EncOfOnes(vectorOnesforElement_ct, vectorOnesforTag_ct);

    //TODO send {p, q, r, g, pk_E, pk_F} to other parties over network in serialized format
    ret = SendInitializedParamsToAllServers();

    //Other parties must store them in their own pir_common.cpp file

    return ret;
}

static int InitSrv_beta(){
    int ret = 0;
    // Initialize random number generator
    // good-ish seed source: std::random_device (combine two samples)
    std::random_device rd;
    unsigned long seed = (static_cast<unsigned long>(rd()) << 1) ^ rd(); // Seed for RNG
    rng.seed(seed); // seed() seeds the gmp_randclass

    //Initialize sockets for communication with server alpha
    ret = InitAcceptingSocket(BETA_LISTENING_TO_ALPHA_PORT, &sock_beta_alpha_srv, &sock_beta_alpha_con);

    if (ret != 0) {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Cannot establish communication with Server Alpha!!");
        ret = -1;
        goto exit;
    }

    ret = InitAcceptingSocket(BETA_LISTENING_TO_GAMMA_PORT, &sock_beta_gamma_srv, &sock_beta_gamma_con);

    if (ret != 0) {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Cannot establish communication with Server Gamma!!");
        ret = -1;
        goto exit;
    }
    
    PrintLog(LOG_LEVEL_INFO, __FILE__, __LINE__, "Server Beta initialization complete");
exit:
    if (ret != 0){
        FinSrv_beta();
    }

    return ret;
}

static int CreateRandomDatabase(){
    int ret = 0;
    mpz_class rand_block_content;
    size_t count;
    plain_db_entry random_entry, read_entry;

    PrintLog(LOG_LEVEL_SPECIAL, __FILE__, __LINE__, "Creating database with random content:"+ DATABASE_LOCATION_BETA);
    pdb.open(pdb_filename, std::ios::in | std::ios::out | std::ios::binary | std::ios::app);

    for (uint64_t i = 0; i < N; ++i) {
        memset(random_entry.element, 0, sizeof(random_entry.element));
        rand_block_content = rng.get_z_bits(PLAINTEXT_PIR_BLOCK_DATA_SIZE);
        mpz_export(random_entry.element, &count, 1, 1, 1, 0, rand_block_content.get_mpz_t());
        insert_pdb_entry(pdb, i, random_entry);
    }

    #if 0
    //Verify 
    for (uint64_t i = 0; i < N; ++i)
    {
        read_pdb_entry(pdb, i, read_entry);

        mpz_t tmp;
        mpz_init(tmp);
        mpz_import(tmp, sizeof(read_entry.element), 1, 1, 1, 0, read_entry.element);
        mpz_class imported_value(tmp);
        mpz_clear(tmp);

        if (imported_value != rand_block_content[i])
        {
            PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Mismatch at location " + std::to_string(i) + " Expected: " + rand_block_content[i].get_str() + "\nRead: " + imported_value.get_str());
        }
    }
    #endif

    pdb.close();

    PrintLog(LOG_LEVEL_SPECIAL, __FILE__, __LINE__, "Database creation complete");

    return ret;
}

/* Its current transcripts are:
13-09-2025 17:41:56:987] [server_beta.cpp:390] TRACE: Number of flushed item is: 199999
[13-09-2025 17:43:31:250] [server_beta.cpp:390] TRACE: Number of flushed item is: 299999
[13-09-2025 17:45:04:080] [server_beta.cpp:390] TRACE: Number of flushed item is: 399999
[13-09-2025 17:46:36:968] [server_beta.cpp:390] TRACE: Number of flushed item is: 499999
[13-09-2025 17:48:10:243] [server_beta.cpp:390] TRACE: Number of flushed item is: 599999
[13-09-2025 17:49:43:623] [server_beta.cpp:390] TRACE: Number of flushed item is: 699999
[13-09-2025 17:51:16:732] [server_beta.cpp:390] TRACE: Number of flushed item is: 799999
[13-09-2025 17:52:50:328] [server_beta.cpp:390] TRACE: Number of flushed item is: 899999
[13-09-2025 17:54:23:701] [server_beta.cpp:390] TRACE: Number of flushed item is: 999999

From it, it looks like it will take more than 30hours, even if used 16-cores simultaneously.
Most expensive operation is exponentiation
*/
/* TODO: This function, currently does not matches with the diagram of the paper. It has to be updated or the diagram has to be modified. */
static int PerEpochOperations_beta(){
    int ret = 0;
    // RNG: mt19937 seeded from random_device
    std::random_device rd;
    std::mt19937 gen(rd());

    PrintLog(LOG_LEVEL_INFO, __FILE__, __LINE__, "Server Beta: Starting PerEpochOperations sequence");

    /* First of all retrieve all the one-time initialized materials from the saved location */
    p = import_from_file_to_mpz_class(ONE_TIME_MATERIALS_LOCATION_BETA + "p.bin");
    q = import_from_file_to_mpz_class(ONE_TIME_MATERIALS_LOCATION_BETA + "q.bin");
    g = import_from_file_to_mpz_class(ONE_TIME_MATERIALS_LOCATION_BETA + "g.bin");
    g_q = import_from_file_to_mpz_class(ONE_TIME_MATERIALS_LOCATION_BETA + "g_q.bin");
    r = import_from_file_to_mpz_class(ONE_TIME_MATERIALS_LOCATION_BETA + "r.bin");
    pk_E = import_from_file_to_mpz_class(ONE_TIME_MATERIALS_LOCATION_BETA + "pk_E.bin");
    sk_E = import_from_file_to_mpz_class(ONE_TIME_MATERIALS_LOCATION_BETA + "sk_E.bin");
    pk_E_q = import_from_file_to_mpz_class(ONE_TIME_MATERIALS_LOCATION_BETA + "pk_E_q.bin");
    sk_E_q = import_from_file_to_mpz_class(ONE_TIME_MATERIALS_LOCATION_BETA + "sk_E_q.bin");
    Serial::DeserializeFromFile(ONE_TIME_MATERIALS_LOCATION_BETA + "FHEcryptoContext.bin", FHEcryptoContext, SerType::BINARY);
    Serial::DeserializeFromFile(ONE_TIME_MATERIALS_LOCATION_BETA + "pk_F.bin", pk_F, SerType::BINARY);
    Serial::DeserializeFromFile(ONE_TIME_MATERIALS_LOCATION_BETA + "sk_F.bin", sk_F, SerType::BINARY);
    Serial::DeserializeFromFile(ONE_TIME_MATERIALS_LOCATION_BETA + "vectorOnesforElement_ct.bin", vectorOnesforElement_ct, SerType::BINARY);
    Serial::DeserializeFromFile(ONE_TIME_MATERIALS_LOCATION_BETA + "vectorOnesforTag_ct.bin", vectorOnesforTag_ct, SerType::BINARY);

    PrintLog(LOG_LEVEL_TRACE, __FILE__, __LINE__, "Server Beta: Loaded one-time initialization materials");

    /* Start sync messages to both the servers, so that everyone is in sync regarding re-initialization process for the epoch */
    (void)sendAll(sock_beta_alpha_con, start_reinit_for_epoch_message.c_str(), start_reinit_for_epoch_message.size());
    (void)sendAll(sock_beta_gamma_con, start_reinit_for_epoch_message.c_str(), start_reinit_for_epoch_message.size());

    pdb.open(pdb_filename, std::ios::in | std::ios::binary | std::ios::app);
    D_K.open(DK_filename, std::ios::in | std::ios::out | std::ios::binary | std::ios::trunc);
    D_alpha.open(D_alpha_filename, std::ios::in | std::ios::out | std::ios::binary | std::ios::trunc);
    D_gamma.open(D_gamma_filename, std::ios::in | std::ios::out | std::ios::binary | std::ios::trunc);

    // 1. Randomly choose Rho in ZZ_((q-1)/2)*
    Rho = rng.get_z_range(((q-1)/2)) + 1;

    /* 2. TODO. And then prepare E_q(Rho) */
    E_q_Rho = ElGamal_q_encrypt(Rho, pk_E_q);

    /* During each request, the client first fetches this from the server_beta, currently storing them in the disk */
    export_to_file_from_mpz_class(PER_EPOCH_MATERIALS_LOCATION_BETA + "Rho.bin", Rho);
    export_to_file_from_mpz_class(PER_EPOCH_MATERIALS_LOCATION_BETA + "E_q_Rho_1.bin", E_q_Rho.first);
    export_to_file_from_mpz_class(PER_EPOCH_MATERIALS_LOCATION_BETA + "E_q_Rho_2.bin", E_q_Rho.second);

    // 3.1 build SS = {1, 2, ..., (N + sqrt_N))}
    uint64_t M = N + sqrt_N;
    std::vector<uint64_t> SS;
    SS.reserve(M);
    for (uint64_t i = 1; i <= M; ++i){
        SS.push_back(i);
    }
    
    if (SS.empty()){
        ret = -1; // nothing to do for N == 0
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Error nothing to do per-epoch, since N is: 0");
        goto exit;
    }

    // 3.2 Clear SetPhi 
    SetPhi.clear();

    /* Create a random shuffling of the list */
    std::shuffle(SS.begin(), SS.end(), gen);

    PrintLog(LOG_LEVEL_TRACE, __FILE__, __LINE__, "Server Beta: Shuffling complete");

    // repeatedly pick a random index in [0, SS.size()-1], print element,
    // then remove it by swapping with the last element and pop_back()

    for (uint64_t iter = 0; iter < M;) {
        //#pragma omp parallel for
        #warning multi-threaded shuffling is not currently working
        for (uint64_t j = 0; (j < NUM_ITEMS_IN_TMP_BUF); ++j)
        {
            size_t send_size = 0;
            mpz_t tmp;
            mpz_class Rho_pow_I;
            mpz_class T_I;
            mpz_class d, d_alpha, d_gamma;
            plain_db_entry read_entry;
            mpz_init(tmp);
            mpz_class d_part, d_alpha_part, d_gamma_part;
            mpz_class mask;
            uint64_t I;
            mpz_class mpz_I;
            unsigned char net_buf_local[(P_BITS/8)];

            // 4. Instead of randomly choose an index I from SS, choose unique index from an already shuffled SS[]
            I = SS[(iter+j)];
            mpz_I  = mpz_class(I);

            /* TODO: Temporarily store the index-location mapping. For the purpose of testing */
#if TEST_SHUFF_DB_FETCH            
            TMP_IDX_LOC_MAP[(I - 1)] = (iter+j); // Since I = 0, is not a valid index
            PrintLog(LOG_LEVEL_TRACE, __FILE__, __LINE__, "Keeping index I: " + std::to_string(I) + " at location: " + std::to_string((iter+j)));
#endif

            // 5. Compute T_I = g^{Rho^I mod p}
            // These two exponentiations takes a long time
            mpz_powm(Rho_pow_I.get_mpz_t(), Rho.get_mpz_t(), mpz_I.get_mpz_t(), q.get_mpz_t());
            mpz_powm(T_I.get_mpz_t(), g.get_mpz_t(), Rho_pow_I.get_mpz_t(), p.get_mpz_t());

            if (I <= N)
            {
                // 6. Compute d = (block_I || I)
                #pragma omp critical
                {
                    /* Index D[I] is located as location (I-1) */
                    // Accessing the disk, randomly takes a long time
                    read_pdb_entry(pdb, (I - 1), read_entry);
                }

                mpz_import(tmp, sizeof(read_entry.element), 1, 1, 1, 0, read_entry.element);
                d = mpz_class(tmp);

                mpz_mul_2exp(d.get_mpz_t(), d.get_mpz_t(), log_N);               // Left shift log_N-bits, so that index can be appended next
                mpz_ior(d.get_mpz_t(), d.get_mpz_t(), mpz_I.get_mpz_t()); // Attach the index at the end
            }
            else
            {
                // 7. Choose d as {0}^{B+log_N} and append T_I to SetPhi
                d = mpz_class(0);
                #pragma omp critical
                {
                    SetPhi.push_back(T_I);
                }
            }
            
            /* 8.1 First create a random number as the secret-share for server_alpha */
            #pragma omp critical
            {
                d_alpha = rng.get_z_bits((PLAINTEXT_PIR_BLOCK_DATA_SIZE + log_N) - 1); /* Since the secret share must be almost half of the original number, make it one bit smaller */
            }

            /* 9.Convert T_I to cuckoo hash key and save that to the buffer */
            mpz_export(net_buf_local, &send_size, 1, 1, 1, 0, T_I.get_mpz_t());
            convert_buf_to_item_type2((const unsigned char*)net_buf_local, (P_BITS/8), TMP_KEY_BUF[(iter+j) % NUM_ITEMS_IN_TMP_BUF]);
            #if 0
            PrintLog(LOG_LEVEL_TRACE, __FILE__, __LINE__, "For item: " + std::to_string(I) + " converted the tag value: " + T_I.get_str() + " to Kuku-key value: ");
            std::cout << "[ ";
            for (const auto& byte : TMP_KEY_BUF[(iter+j) % NUM_ITEMS_IN_TMP_BUF]) {
                // static_cast<int> is CRITICAL. 
                // Without it, cout tries to print the ASCII character.
                std::cout << static_cast<int>(byte) << " "; 
            }
            std::cout << "]" << std::endl;
            #else
            PrintLog(LOG_LEVEL_TRACE, __FILE__, __LINE__, "Iteration: " + std::to_string(iter+j) + " item: " + std::to_string(I));
            #endif
            //(void)sendAll(sock_beta_alpha_con, net_buf, send_size);
            //(void)sendAll(sock_beta_gamma_con, net_buf, send_size);

            /* 10.1 Store the d_alpha share in the local buffer */
            mpz_export(&TMP_D_ALPHA_BUF[(NUM_BYTES_PER_SDB_ELEMENT * ((iter+j) % NUM_ITEMS_IN_TMP_BUF))], &send_size, 1, 1, 1, 0, d_alpha.get_mpz_t());
            //(void)sendAll(sock_beta_alpha_con, net_buf, send_size);

            /* For some numbers, gmp exporting an additional byte. This is a corresponding fix */
            if (send_size != NUM_BYTES_PER_SDB_ELEMENT)
            {
                /* For some reason, for dummy elements exporting d_gamma takes one more byte and that is causing the problem. Hence, first exporting that to a different buffer and then copy the content from there. */
                mpz_export(net_buf_local, &send_size, 1, 1, 1, 0, d_alpha.get_mpz_t());
                memcpy(&TMP_D_ALPHA_BUF[(NUM_BYTES_PER_SDB_ELEMENT * ((iter+j) % NUM_ITEMS_IN_TMP_BUF))], (net_buf_local + 1), NUM_BYTES_PER_SDB_ELEMENT);
            }

            /* 8.2 Create the second share for server_gamma */
            // d_gamma = (d - d_alpha);/* Another share */ But this is creating error while combining homomorphically

            // PrintLog(LOG_LEVEL_TRACE, __FILE__, __LINE__, "For I = " + std::to_string(I) + " value of d: " + d.get_str(16) + " and d_alpha: " + d_alpha.get_str(16));

            d_gamma = mpz_class(0);
            mask = mpz_class((1 << PLAINTEXT_FHE_BLOCK_SIZE) - 1);

            for (unsigned int i = 0; i < TOTAL_NUM_FHE_BLOCKS_PER_ELEMENT; i++)
            {
                /* Extract least significant PLAINTEXT_FHE_BLOCK_SIZE-bits of d and d_alpha */
                d_alpha_part = (d_alpha & mask);
                d_part = (d & mask);

                /* Compute the difference between two parts. And take only PLAINTEXT_FHE_BLOCK_SIZE-bits */
                d_gamma_part = (d_part - d_alpha_part) & mask;

                /* Append the part at the proper location */
                d_gamma = (d_gamma | d_gamma_part);

                mask = mask << PLAINTEXT_FHE_BLOCK_SIZE;
            }

            /* 10.2 Store the d_gamma share in the local buffer */
            mpz_export(&TMP_D_GAMMA_BUF[(NUM_BYTES_PER_SDB_ELEMENT * ((iter+j) % NUM_ITEMS_IN_TMP_BUF))], &send_size, 1, 1, 1, 0, d_gamma.get_mpz_t());

            /* For some numbers, gmp exporting an additional byte. This is a corresponding fix */
            if (send_size != NUM_BYTES_PER_SDB_ELEMENT)
            {
                /* For some reason, for dummy elements exporting d_gamma takes one more byte and that is causing the problem. Hence, first exporting that to a different buffer and then copy the content from there. */
                mpz_export(net_buf_local, &send_size, 1, 1, 1, 0, d_gamma.get_mpz_t());
                memcpy(&TMP_D_GAMMA_BUF[(NUM_BYTES_PER_SDB_ELEMENT * ((iter+j) % NUM_ITEMS_IN_TMP_BUF))], (net_buf_local + 1), NUM_BYTES_PER_SDB_ELEMENT);
            }
        }

        //PrintLog(LOG_LEVEL_TRACE, __FILE__, __LINE__, "d_gamma: " + d_gamma.get_str(16) + " exported size: " + std::to_string(send_size));

        //(void)sendAll(sock_beta_gamma_con, net_buf, send_size);

        // 11. Remove chosen element from SS (order not preserved)
        // This step is not required, since SS[] is already shuffled and we are choosing all the indices only once, one by one

        iter += NUM_ITEMS_IN_TMP_BUF;

        /* Flush into the disk */
        D_K.write((const char *)TMP_KEY_BUF, sizeof(TMP_KEY_BUF));
        D_alpha.write((const char *)TMP_D_ALPHA_BUF, sizeof(TMP_D_ALPHA_BUF));
        D_gamma.write((const char *)TMP_D_GAMMA_BUF, sizeof(TMP_D_GAMMA_BUF));

        /* Clear the buffers */
        memset(TMP_KEY_BUF, 0, ((sizeof(item_type)) * NUM_ITEMS_IN_TMP_BUF));
        memset(TMP_D_ALPHA_BUF, 0, (NUM_BYTES_PER_SDB_ELEMENT * NUM_ITEMS_IN_TMP_BUF));
        memset(TMP_D_GAMMA_BUF, 0, (NUM_BYTES_PER_SDB_ELEMENT * NUM_ITEMS_IN_TMP_BUF));

        PrintLog(LOG_LEVEL_TRACE, __FILE__, __LINE__, "Number of flushed item is: " + std::to_string(iter));
    }

    /* Files are required to close to ensure everything is actually flushed to the disk */
    pdb.close();
    D_alpha.close();
    D_gamma.close();
    D_K.close();

    if (!save_mpz_vector(SetPhi, SetPhi_filename))
    { 
        PrintLog(LOG_LEVEL_TRACE, __FILE__, __LINE__, "Error while exporting the SetPhi");
        goto exit;
    }

    PrintLog(LOG_LEVEL_TRACE, __FILE__, __LINE__, "Completed preparing the database shares and key-database");

    PrintLog(LOG_LEVEL_SPECIAL, __FILE__, __LINE__, "Please manually transfer the shuffled and secret-shared databases to server_alpha and server_gamma.");
    PrintLog(LOG_LEVEL_SPECIAL, __FILE__, __LINE__, "Please transfer the secret-shared database to server_alpha, which is located at:" + D_alpha_filename);
    PrintLog(LOG_LEVEL_SPECIAL, __FILE__, __LINE__, "Please transfer the secret-shared database to server_gama, which is located at:" + D_gamma_filename);
    PrintLog(LOG_LEVEL_SPECIAL, __FILE__, __LINE__, "Please transfer the key-file to both server_alpha and server_gamma, which is located at:" + DK_filename);

    PrintLog(LOG_LEVEL_SPECIAL, __FILE__, __LINE__, "Please press enter to continue..!!");

    std::cin.get();

    /* Send ready message to Server Alpha and Server Gamma */
    (void)sendAll(sock_beta_alpha_con, completed_reinit_for_epoch_message.c_str(), completed_reinit_for_epoch_message.size());
    (void)sendAll(sock_beta_gamma_con, completed_reinit_for_epoch_message.c_str(), completed_reinit_for_epoch_message.size());

    PrintLog(LOG_LEVEL_INFO, __FILE__, __LINE__, "Server Beta: Completed PerEpochOperations for new epoch");
exit:
    //mpz_clear(tmp);
    pdb.close();
    D_alpha.close();
    D_gamma.close();
    D_K.close();
    
    return ret;
}

static int FinSrv_beta(){
    int ret = 0;

    // Close the sockets
    if (sock_beta_alpha_srv != -1) {
        close(sock_beta_alpha_srv);
        sock_beta_alpha_srv = -1;
    }
    if (sock_beta_alpha_con != -1) {
        close(sock_beta_alpha_con);
        sock_beta_alpha_con = -1;
    }

    PrintLog(LOG_LEVEL_SPECIAL, __FILE__, __LINE__, "Finalized Server Beta");

    return ret;
}

static int ShelterTagDetermination_beta(){
    int ret = 0;
    size_t received_sz = 0;
    int ret_recv = 0;
    std::pair<mpz_class, mpz_class> E_q_Rho_pow_I__mul__h_C;
    std::pair<mpz_class, mpz_class> E_g_pow_Rho_pow_I__mul_a_mul_c;
    std::pair<mpz_class, mpz_class> E_g_pow_Rho_pow_I__mul__h_C;
    mpz_class Rho_pow_I__mul__h_C;
    mpz_class g_pow_Rho_pow_I__mul__h_C;
    mpz_class g_pow_Rho_pow_I__mul_a_mul_c;
    mpz_class widehat_T_I;    

    /* Step 2.3.2.1 of the sequence diagram */
    // Receive E_q_Rho_pow_I__mul__h_C.first
    ret_recv = recvAll(sock_beta_client_con, net_buf, sizeof(net_buf), &received_sz);
    if (ret_recv != 0)
    {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Failed to receive E_q_Rho_pow_I__mul__h_C.first from the client");
        return -1;
    }
    /* From now on, starting client request processing */
    PrintLog(LOG_LEVEL_TRACE, __FILE__, __LINE__, "Server Beta starts client request processing from this point");

    E_q_Rho_pow_I__mul__h_C.first = mpz_class(std::string(net_buf, received_sz));

    /* Step 2.3.2.1 of the sequence diagram */
    // Receive E_q_Rho_pow_I__mul__h_C.second
    ret_recv = recvAll(sock_beta_client_con, net_buf, sizeof(net_buf), &received_sz);
    if (ret_recv != 0)
    {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Failed to receive E_q_Rho_pow_I__mul__h_C.second from the client");
        return -1;
    }
    E_q_Rho_pow_I__mul__h_C.second = mpz_class(std::string(net_buf, received_sz));

    /* Step 3 */
    Rho_pow_I__mul__h_C = ElGamal_q_decrypt(E_q_Rho_pow_I__mul__h_C, sk_E_q);

    /* Step 4.1 perform g^{Rho_pow_I__mul__h_C} mod p */
    mpz_powm(g_pow_Rho_pow_I__mul__h_C.get_mpz_t(), g.get_mpz_t(), Rho_pow_I__mul__h_C.get_mpz_t(), p.get_mpz_t());
    /* Step 4.2 Encrypts under ElGamal encryption in GG */
    E_g_pow_Rho_pow_I__mul__h_C = ElGamal_encrypt(g_pow_Rho_pow_I__mul__h_C, pk_E);
    /* Step 4.3 Send both the coponents of the ciphtertext to Server alpha */
    (void)sendAll(sock_beta_alpha_con, E_g_pow_Rho_pow_I__mul__h_C.first.get_str().c_str(), E_g_pow_Rho_pow_I__mul__h_C.first.get_str().size());
    (void)sendAll(sock_beta_alpha_con, E_g_pow_Rho_pow_I__mul__h_C.second.get_str().c_str(), E_g_pow_Rho_pow_I__mul__h_C.second.get_str().size());

    /* Step 10.3.2 Receive the first component of E_g_pow_Rho_pow_I__mul_a_mul_c from Server Gamma */
    ret_recv = recvAll(sock_beta_gamma_con, net_buf, sizeof(net_buf), &received_sz);
    if (ret_recv != 0)
    {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Failed to receive E_g_pow_Rho_pow_I__mul_a_mul_c.first from the Server Gamma");
        return -1;
    }
    E_g_pow_Rho_pow_I__mul_a_mul_c.first = mpz_class(std::string(net_buf, received_sz));

    /* Step 10.4.2 Receive the second component of E_g_pow_Rho_pow_I__mul_a_mul_c from Server Gamma */
    ret_recv = recvAll(sock_beta_gamma_con, net_buf, sizeof(net_buf), &received_sz);
    if (ret_recv != 0)
    {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Failed to receive E_g_pow_Rho_pow_I__mul_a_mul_c.second from the Server Gamma");
        return -1;
    }
    E_g_pow_Rho_pow_I__mul_a_mul_c.second = mpz_class(std::string(net_buf, received_sz));

    /* Step 11.1 Decrypt E_g_pow_Rho_pow_I__mul_a_mul_c  */
    g_pow_Rho_pow_I__mul_a_mul_c = ElGamal_decrypt(E_g_pow_Rho_pow_I__mul_a_mul_c, sk_E);

    /* Step 11.2 Determine the shelter tag, \widehat{T_I}  */
    widehat_T_I = (g_pow_Rho_pow_I__mul_a_mul_c * b) % p;

    /* Step 11.3 Determine \widehat{t_I} */
    widehat_t_I = g_pow_Rho_pow_I__mul_a_mul_c % r;

    /* TODO: This step is only for verification */
#if 0
    mpz_class RhoExpI, gExp_RhoExpI, I, a, c;

    printf("Enter the value of I (base 10): ");
    mpz_inp_str(I.get_mpz_t(), stdin, 10);
    printf("Enter the value of a (base 10): ");
    mpz_inp_str(a.get_mpz_t(), stdin, 10);
    printf("Enter the value of c (base 10): ");
    mpz_inp_str(c.get_mpz_t(), stdin, 10);

    mpz_powm(RhoExpI.get_mpz_t(), Rho.get_mpz_t(), mpz_class(I).get_mpz_t(), q.get_mpz_t()); // Since this is an exponenet
    PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Calculated expected RhoExpI: " + RhoExpI.get_str());

    mpz_powm(gExp_RhoExpI.get_mpz_t(), g.get_mpz_t(), RhoExpI.get_mpz_t(), p.get_mpz_t());
    PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Calculated expected RhoExpI: " + gExp_RhoExpI.get_str());

    mpz_class expected_widehat_T_I = (gExp_RhoExpI * a * b * c) % p;
    if (expected_widehat_T_I != widehat_T_I)
    {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Calculated widehat_T_I does not match: expected " + expected_widehat_T_I.get_str() + ", got " + widehat_T_I.get_str());
    }
    else
    {
        PrintLog(LOG_LEVEL_INFO, __FILE__, __LINE__, "widehat_T_I matched with expected value. Haha..Thank you..:) :) ");
    }
    /* Upto this point, for the verification of shelter_tag_determination flow */
#endif

    return ret;
}

static int ObliviouslySearchShelter_beta() {
    // Set up variables
    Fss fClient, fServer;
    ServerKeyEq K_alpha;
    ServerKeyEq K_gamma;
    size_t serializedFssSize;
    size_t received_sz = 0;
    Ciphertext<DCRTPoly> tmp_ct;
    mpz_class tmp_pt;
    int ret = 0;

    PrintLog(LOG_LEVEL_TRACE, __FILE__, __LINE__, "Generating DPF keys");

    // Step 1.1 Initialize client, use 64 bits in domain as example
    initializeClient(&fClient, R_BITS, 2); // If bit length is not set properly, then incorrect answer will be returned

    // Step 1.2 Generate keys for equality FSS test
    generateTreeEq(&fClient, &K_alpha, &K_gamma, widehat_t_I, 1);//So that the point function will evaluate as 1 at location i, and zero elsewhere

    // Step 1.3 Initialize server structure
    initializeServer(&fServer, &fClient);

    /* 2.a.1 Genrate FSS-key parts and send them to the server_alpha and server_gamma */
    serializedFssSize = serializeFssAndServerKeyEq(fServer, K_alpha, net_buf, sizeof(net_buf));

    (void)sendAll(sock_beta_alpha_con, net_buf, serializedFssSize);

    serializedFssSize = serializeFssAndServerKeyEq(fServer, K_gamma, net_buf, sizeof(net_buf));

    (void)sendAll(sock_beta_gamma_con, net_buf, serializedFssSize);

#if TEST_SHELTER_FOUND
    /* This line is only for testing purpose  */
    PrintLog(LOG_LEVEL_SPECIAL, __FILE__, __LINE__, "Fore verification only. Enter this search tag in Server_Alpha and Server_Beta prompt to make DPF search successful: " + widehat_t_I.get_str());
#endif

    /************************ Refresh fnd_ct_element ciphertext ***********************/

    /* Receive the ciphertext fnd_ct_element */
    ret = recvAll(sock_beta_alpha_con, net_buf, sizeof(net_buf), &received_sz);
    if (ret != 0)
    {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Failed to receive fnd_ct_element from Server Alpha");
        return -1;
    }
    Serial::DeserializeFromString(tmp_ct, std::string(net_buf, received_sz));

    /* Decrypt then re-encrypt and send */
    FHE_Dec_SDBElement(tmp_ct, tmp_pt);
    tmp_ct = FHE_Enc_SDBElement(tmp_pt);

    (void)sendAll(sock_beta_alpha_con, Serial::SerializeToString(tmp_ct).c_str(), Serial::SerializeToString(tmp_ct).size());

    /************************ Refresh fnd_ct_tag ciphertext ***********************/

    /* Receive the ciphertext fnd_ct_tag */
    ret = recvAll(sock_beta_alpha_con, net_buf, sizeof(net_buf), &received_sz);
    if (ret != 0)
    {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Failed to receive fnd_ct_tag from Server Alpha");
        return -1;
    }
    Serial::DeserializeFromString(tmp_ct, std::string(net_buf, received_sz));
    
    /* Decrypt then re-encrypt and send */
    FHE_Dec_Tag(tmp_ct, tmp_pt);
    tmp_ct = FHE_Enc_Tag(tmp_pt);

    (void)sendAll(sock_beta_alpha_con, Serial::SerializeToString(tmp_ct).c_str(), Serial::SerializeToString(tmp_ct).size());

    return ret;
}

static int ShelterUpdate_beta(){
    int ret = 0;
    size_t received_sz = 0;
    int ret_recv = 0;
    std::pair<mpz_class, mpz_class> E_T_star_a_dashed_c_dashed_h_gamma0;
    mpz_class T_star_a_dashed_c_dashed_h_gamma0, T_star_a_dashed_b_dashed_c_dashed_h_gamma0;
    
    /* 1.b. Randomly select b_dashed */
    mpz_class b_dashed = ElGamal_randomGroupElement();

    /* 3.4.1 Receive first component of E(T_*.a'.c'.h_gamma0) */
    ret = recvAll(sock_beta_gamma_con, net_buf, sizeof(net_buf), &received_sz);
    if (ret != 0) {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Failed to receive first component of E(T_*.a'.c'.h_gamma0)) from Server Gamma");
        close(sock_beta_gamma_con);
        goto exit;
    }
    E_T_star_a_dashed_c_dashed_h_gamma0.first = mpz_class(std::string(net_buf, received_sz));

    /* 3.4.2 Receive second component of E(T_*.a'.c'.h_gamma0) */
    ret = recvAll(sock_beta_gamma_con, net_buf, sizeof(net_buf), &received_sz);
    if (ret != 0) {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Failed to receive second component of E(T_*.a'.c'.h_gamma0)) from Server Gamma");
        close(sock_beta_gamma_con);
        goto exit;
    }
    E_T_star_a_dashed_c_dashed_h_gamma0.second = mpz_class(std::string(net_buf, received_sz));

    /* 4.1 Decrypt E(T_*.a'.c'.h_gamma0) to T_*.a'.c'.h_gamma0 */
    T_star_a_dashed_c_dashed_h_gamma0 = ElGamal_decrypt(E_T_star_a_dashed_c_dashed_h_gamma0, sk_E);

    /* 4.2 Compute T_*.a'.b'.c'.h_gamma0 */
    T_star_a_dashed_b_dashed_c_dashed_h_gamma0 = (T_star_a_dashed_c_dashed_h_gamma0 * b_dashed) % p;

    /* 4.3 Return T_*.a'.b'.c'.h_gamma0 to the server gamma */
    (void)sendAll(sock_beta_gamma_con, T_star_a_dashed_b_dashed_c_dashed_h_gamma0.get_str().c_str(), T_star_a_dashed_b_dashed_c_dashed_h_gamma0.get_str().size());

exit:
    return ret;
}

static int ProcessClientRequest_beta(){
    int ret = -1;
    struct sockaddr_in address;

    PrintLog(LOG_LEVEL_INFO, __FILE__, __LINE__, "Server Beta: Starting Processing client request");

    /* First of all retrieve all the one-time initialized materials from the saved location */
    p = import_from_file_to_mpz_class(ONE_TIME_MATERIALS_LOCATION_BETA + "p.bin");
    q = import_from_file_to_mpz_class(ONE_TIME_MATERIALS_LOCATION_BETA + "q.bin");
    g = import_from_file_to_mpz_class(ONE_TIME_MATERIALS_LOCATION_BETA + "g.bin");
    g_q = import_from_file_to_mpz_class(ONE_TIME_MATERIALS_LOCATION_BETA + "g_q.bin");
    r = import_from_file_to_mpz_class(ONE_TIME_MATERIALS_LOCATION_BETA + "r.bin");
    pk_E = import_from_file_to_mpz_class(ONE_TIME_MATERIALS_LOCATION_BETA + "pk_E.bin");
    sk_E = import_from_file_to_mpz_class(ONE_TIME_MATERIALS_LOCATION_BETA + "sk_E.bin");
    pk_E_q = import_from_file_to_mpz_class(ONE_TIME_MATERIALS_LOCATION_BETA + "pk_E_q.bin");
    sk_E_q = import_from_file_to_mpz_class(ONE_TIME_MATERIALS_LOCATION_BETA + "sk_E_q.bin");
    E_q_Rho.first = import_from_file_to_mpz_class(PER_EPOCH_MATERIALS_LOCATION_BETA + "E_q_Rho_1.bin");
    E_q_Rho.second = import_from_file_to_mpz_class(PER_EPOCH_MATERIALS_LOCATION_BETA + "E_q_Rho_2.bin");
    Rho = import_from_file_to_mpz_class(PER_EPOCH_MATERIALS_LOCATION_BETA + "Rho.bin");

    Serial::DeserializeFromFile(ONE_TIME_MATERIALS_LOCATION_BETA + "FHEcryptoContext.bin", FHEcryptoContext, SerType::BINARY);
    Serial::DeserializeFromFile(ONE_TIME_MATERIALS_LOCATION_BETA + "pk_F.bin", pk_F, SerType::BINARY);
    Serial::DeserializeFromFile(ONE_TIME_MATERIALS_LOCATION_BETA + "sk_F.bin", sk_F, SerType::BINARY);
    Serial::DeserializeFromFile(ONE_TIME_MATERIALS_LOCATION_BETA + "vectorOnesforElement_ct.bin", vectorOnesforElement_ct, SerType::BINARY);
    Serial::DeserializeFromFile(ONE_TIME_MATERIALS_LOCATION_BETA + "vectorOnesforTag_ct.bin", vectorOnesforTag_ct, SerType::BINARY);

    if (!load_mpz_vector(SetPhi, SetPhi_filename))
    {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Server Beta: Error while loding the SetPhi");
        goto exit;
    }

    PrintLog(LOG_LEVEL_TRACE, __FILE__, __LINE__, "Server Beta: Loaded one-time initialization materials");

    //Always initialize them
    K = 0;
    b = 1;/* Initialize with 1, so that, during the first iteration b = b'*(b^-1) = b'*1 = b' */

    /* For debugging only */
#if 1
    /* Populate the shelter, with random elements */
    for(size_t k = 0; k < sqrt_N; k++) {
        //if (!std::filesystem::exists(SHELTER_STORING_LOCATION + "sh[" + std::to_string(k) + "].ct")) {
        if (1) {
            // Generate random block_content of PLAINTEXT_PIR_BLOCK_DATA_SIZE bits of random | k as the block index
            mpz_class data = rng.get_z_bits(PLAINTEXT_PIR_BLOCK_DATA_SIZE);
            mpz_class index = rng.get_z_bits(log_N);
            Ciphertext<DCRTPoly> tmp_ct = FHE_Enc_SDBElement(( data << log_N) | index);

            /* Printing, so that afer the DPF search it can be verified */
            if (k == 0){
                PrintLog(LOG_LEVEL_TRACE, __FILE__, __LINE__, "The value of data, which should be found during DPF search experiment is: " + data.get_str()+ " index: " + index.get_str());
            }

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

            if (k == 0){
                PrintLog(LOG_LEVEL_TRACE, __FILE__, __LINE__, "Size of the FHE-ciphertext: "+ std::to_string(std::filesystem::file_size(SHELTER_STORING_LOCATION + "sh[" + std::to_string(k) + "].ct")));
            }
        }

        mpz_class T_hat = ElGamal_randomGroupElement();
        export_to_file_from_mpz_class(SHELTER_STORING_LOCATION + "sh[" + std::to_string(k) + "].T_hat", T_hat);
        export_to_file_from_mpz_class(SHELTER_STORING_LOCATION + "sh[" + std::to_string(k) + "].t_hat", (T_hat % r));
    }
    
    /* Load the shelter into the RAM */
    for(size_t k = 0; k < sqrt_N; k++) {
        sh[k].element_FHE_ct = import_from_file_to_mpz_class(SHELTER_STORING_LOCATION + "sh[" + std::to_string(k) + "].ct");

        /* Generate the tags and keep them in the variable, which will be used for DPF search */
        sh[k].tag = import_from_file_to_mpz_class(SHELTER_STORING_LOCATION + "sh[" + std::to_string(k) + "].T_hat");
        sh[k].tag_short = import_from_file_to_mpz_class(SHELTER_STORING_LOCATION + "sh[" + std::to_string(k) + "].t_hat"); // Create a random short tag
    }
#endif    

    while (K < sqrt_N){        
        PrintLog(LOG_LEVEL_INFO, __FILE__, __LINE__, "Waiting for processing the PIR request number: "+ std::to_string(K+1) +" from the client..!!");

        ret = InitAcceptingSocket(BETA_LISTENING_TO_CLIENT_PORT, &sock_beta_client_srv, &sock_beta_client_con);

        if (ret != 0)
        {
            PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Cannot establish communication with the client!!");
            ret = -1;
            goto exit;
        }

        PrintLog(LOG_LEVEL_TRACE, __FILE__, __LINE__, "Sending all the required, publicly known, parameters to the client");

        // Send the required parameters to the client
        (void)sendAll(sock_beta_client_con, p.get_str().c_str(), p.get_str().size());
        (void)sendAll(sock_beta_client_con, q.get_str().c_str(), q.get_str().size());
        (void)sendAll(sock_beta_client_con, g.get_str().c_str(), g.get_str().size());
        (void)sendAll(sock_beta_client_con, g_q.get_str().c_str(), g_q.get_str().size());
        (void)sendAll(sock_beta_client_con, r.get_str().c_str(), r.get_str().size());
        (void)sendAll(sock_beta_client_con, pk_E.get_str().c_str(), pk_E.get_str().size());
        (void)sendAll(sock_beta_client_con, pk_E_q.get_str().c_str(), pk_E_q.get_str().size());
        (void)sendAll(sock_beta_client_con, Serial::SerializeToString(FHEcryptoContext).c_str(), Serial::SerializeToString(FHEcryptoContext).size());
        (void)sendAll(sock_beta_client_con, Serial::SerializeToString(pk_F).c_str(), Serial::SerializeToString(pk_F).size());
        (void)sendAll(sock_beta_client_con, E_q_Rho.first.get_str().c_str(), E_q_Rho.first.get_str().size());
        (void)sendAll(sock_beta_client_con, E_q_Rho.second.get_str().c_str(), E_q_Rho.second.get_str().size());

        ret = ShelterTagDetermination_beta();
        if (ret != 0)
        {
            PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Problem while determining the shelter tag..!!");
            ret = -1;
            goto exit;
        }

        /* For the first request, the shelter is not required to be searched */
        if (K > 0)
        {
            ret = ObliviouslySearchShelter_beta();
            if (ret != 0)
            {
                PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Problem during the shelter search operation..!!");
                ret = -1;
                goto exit;
            }
        }

        ret = SelShuffDBSearchTag_beta();

        if (ret != 0){
            PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Problem during the selecting shuffled database tag..!!");
            ret = -1;
            goto exit;
        }

        #warning Temporary from here
        size_t received_sz = 0;
        mpz_class dec_element_content, dec_element_index, dec_content_and_index, dec_tmp;
        Ciphertext<DCRTPoly> tmp;
        (void)sendAll(sock_beta_gamma_con, Serial::SerializeToString(sk_F).c_str(), Serial::SerializeToString(sk_F).size());

        ret = recvAll(sock_beta_gamma_con, net_buf, sizeof(net_buf), &received_sz);
        if (ret != 0)
        {
            PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Failed to receive FHE Ciphertext requested_element_ct from Server Gamma");
            goto exit;
        }
        Serial::DeserializeFromString(tmp, std::string(net_buf, received_sz));
        FHE_Dec_SDBElement(tmp, dec_tmp);
        PrintLog(LOG_LEVEL_TRACE, __FILE__, __LINE__, "Receivded SR_D: " + dec_tmp.get_str(16));

        ret = recvAll(sock_beta_gamma_con, net_buf, sizeof(net_buf), &received_sz);
        if (ret != 0)
        {
            PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Failed to receive FHE Ciphertext fnd_ct_element from Server Gamma");
            goto exit;
        }
        Serial::DeserializeFromString(tmp, std::string(net_buf, received_sz));
        FHE_Dec_SDBElement(tmp, dec_tmp);
        PrintLog(LOG_LEVEL_TRACE, __FILE__, __LINE__, "Receivded fnd_element: " + dec_tmp.get_str(16));

        ret = recvAll(sock_beta_gamma_con, net_buf, sizeof(net_buf), &received_sz);
        if (ret != 0)
        {
            PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Failed to receive FHE Ciphertext requested_element_ct from Server Gamma");
            goto exit;
        }
        Serial::DeserializeFromString(tmp, std::string(net_buf, received_sz));
        FHE_Dec_SDBElement(tmp, dec_tmp);
        PrintLog(LOG_LEVEL_TRACE, __FILE__, __LINE__, "Receivded SR_sh: " + dec_tmp.get_str(16));

        ret = recvAll(sock_beta_gamma_con, net_buf, sizeof(net_buf), &received_sz);
        if (ret != 0)
        {
            PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Failed to receive FHE Ciphertext requested_element_ct from Server Alpha");
            goto exit;
        }
        Serial::DeserializeFromString(requested_element_ct, std::string(net_buf, received_sz));
        FHE_Dec_SDBElement(requested_element_ct, dec_content_and_index);
        dec_element_content = (dec_content_and_index >> log_N);
        dec_element_index = (dec_content_and_index & ((1U << log_N) - 1U));
        PrintLog(LOG_LEVEL_TRACE, __FILE__, __LINE__, "Selected element is: " + dec_content_and_index.get_str(16));
        PrintLog(LOG_LEVEL_TRACE, __FILE__, __LINE__, "Fetcheched index is: " + dec_element_index.get_str());
        #warning Temporary upto here

        ShelterUpdate_beta();

        /* Close the connection with existing client */
        close(sock_beta_client_srv);
        close(sock_beta_client_con);

        K++;
    }

    PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Current epoch is completed. Please re-perform the per-epoch initialization");

exit:
    if(ret != 0){
        close(sock_beta_client_srv);
        close(sock_beta_client_con);
    }

    return ret;
}

static int SelShuffDBSearchTag_beta(){
    int ret = -1;
    size_t received_sz = 0;

    PrintLog(LOG_LEVEL_TRACE, __FILE__, __LINE__, "Server Beta: Starting SelShuffDBSearchTag sequence");

    /* 1.b.1 Select random h_{\beta 0} */
    mpz_class h_beta0 = ElGamal_randomGroupElement();
    /* Also select its inverse */
    mpz_class h_beta0_1;
    mpz_invert(h_beta0_1.get_mpz_t(), h_beta0.get_mpz_t(), p.get_mpz_t());

    // 2.1.b Receive first component of E(T_I.h_{\alpha 1})
    ret = recvAll(sock_beta_alpha_con, net_buf, sizeof(net_buf), &received_sz);
    if (ret != 0) {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Server Beta: Failed to receive first component of E(T_I.h_{\\alpha 1} from Server Alpha");
        close(sock_beta_alpha_con);
        return ret;
    }
    mpz_class E_T_I_h_alpha1_1 = mpz_class(std::string(net_buf, received_sz));

    //2.2.b Receive second component of E(T_I.h_{\alpha 1})
    ret = recvAll(sock_beta_alpha_con, net_buf, sizeof(net_buf), &received_sz);
    if (ret != 0) {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Server Beta: Failed to receive second component of E(T_I.h_{\\alpha 1} from Server Alpha");
        close(sock_beta_alpha_con);
        return ret;
    }

    //3. Decrypt to T_I.h_{\\alpha 1}
    mpz_class E_T_I_h_alpha1_2 = mpz_class(std::string(net_buf, received_sz));
    std::pair<mpz_class, mpz_class> E_T_I_h_alpha1_local = std::make_pair(E_T_I_h_alpha1_1, E_T_I_h_alpha1_2);
    mpz_class T_I_h_alpha1 = ElGamal_decrypt(E_T_I_h_alpha1_local, sk_E);

    //PrintLog(LOG_LEVEL_INFO, __FILE__, __LINE__, "Decrypted T_I.h_{\\alpha 1}: " + T_I_h_alpha1.get_str());

    /* 4.b.1 Compute T_I.h_{\alpha 1}h_{\beta 0} */
    mpz_class T_I_h_alpha1_h_beta0 = (T_I_h_alpha1*h_beta0) % p;
    /* 4.b.2 Send T_I.h_{\alpha 1}h_{\beta 0} to server alpha */
    (void)sendAll(sock_beta_alpha_con, T_I_h_alpha1_h_beta0.get_str().c_str(), T_I_h_alpha1_h_beta0.get_str().size());

    /* 7. Select a new dummy tag */
    mpz_class T_phi;
    if (!SetPhi.empty()) {
        // move the last element into `T_phi` and remove it from the vector
        T_phi = std::move(SetPhi.back());
        SetPhi.pop_back();
        PrintLog(LOG_LEVEL_TRACE, __FILE__, __LINE__, "Selected new dummy tag is: " + T_phi.get_str());
    }
    else
    {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "SetPhi is empty");
        return -1;
    }

    /* 8.b.1 Compute T_phi.h_{\beta 0} */
    mpz_class T_phi_h_beta0 = (T_phi*h_beta0) % p;
    /* 8.b.2 Send to server gamma */
    (void)sendAll(sock_beta_gamma_con, T_phi_h_beta0.get_str().c_str(), T_phi_h_beta0.get_str().size());

    //11.b.1 Receive FHE-ciphertext of (T_star.h_{\alpha 2}.h_{\beta 0})
    ret = recvAll(sock_beta_alpha_con, net_buf, sizeof(net_buf), &received_sz);
    if (ret != 0) {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Server Beta: Failed to receive FHE ciphertext of T_star.h_{\\alpha 2}.h_{\\beta 0} from Server Alpha");
        close(sock_beta_alpha_con);
        return ret;
    }
    Ciphertext<DCRTPoly> FHE_ct_T_star_h_alpha2_h_beta0;
    Serial::DeserializeFromString(FHE_ct_T_star_h_alpha2_h_beta0, std::string(net_buf, received_sz));

    /* 12.1 Decrypt the FHE-ciphertext */
    mpz_class T_star_h_alpha2_h_beta0;
    FHE_Dec_Tag(FHE_ct_T_star_h_alpha2_h_beta0, T_star_h_alpha2_h_beta0);

    /* 12.2Remove h_{\beta 0} */
    mpz_class T_star_h_alpha2 = (T_star_h_alpha2_h_beta0*h_beta0_1) % p;
    
    /* 13.a Send to server alpha */
    (void)sendAll(sock_beta_alpha_con, T_star_h_alpha2.get_str().c_str(), T_star_h_alpha2.get_str().size());

    return ret;
}

static void TestSelShuffDBSearchTag_beta(){
    int ret = -1;

    ret = SelShuffDBSearchTag_beta();

    return;
}


static void TestBlindedExponentiation() {
    int ret;

    PrintLog(LOG_LEVEL_INFO, __FILE__, __LINE__, "El-Gamal parameters p: " + p.get_str() + " q: " + q.get_str()+ " g: " + g.get_str()+ " pk_E: " + pk_E.get_str()+ " sk_E: " + sk_E.get_str());

    Rho = rng.get_z_range(((q-1)/2)) + 1; // Randomly choose Rho in ZZ_((q-1)/2)*

    mpz_class h = rng.get_z_range(q-1) + 1;//To ensure that the element belongs to ZZ_q*
    mpz_class h_1;
    mpz_invert(h_1.get_mpz_t(), h.get_mpz_t(), q.get_mpz_t());//Compute h^-1 in ZZ_q*

    mpz_class alpha = rng.get_z_range(q-1) + 1;
    mpz_class alpha_1;
    mpz_invert(alpha_1.get_mpz_t(), alpha.get_mpz_t(), q.get_mpz_t());

    //Confirm that the inverses really work
    assert((h * h_1) % q == 1);
    assert((alpha * alpha_1) % q == 1);

    mpz_class I = rng.get_z_range(N-1) + 1;//I will be in the range of [1, N]
    PrintLog(LOG_LEVEL_INFO, __FILE__, __LINE__, "Chosen plain text messages are Rho: " + Rho.get_str() + " I: " + I.get_str() + " h: " + h.get_str() + " h_1: " + h_1.get_str()+ " alpha: " + alpha.get_str() + " alpha_1: " + alpha_1.get_str() );
    mpz_class RhoExpI;
    //mpz_powm(RhoExpI.get_mpz_t(), Rho.get_mpz_t(), I.get_mpz_t(), p.get_mpz_t()); previous
    //mpz_class RhoExpI_h = (RhoExpI * h) % p; previous
    mpz_powm(RhoExpI.get_mpz_t(), Rho.get_mpz_t(), I.get_mpz_t(), q.get_mpz_t());//Since this is an exponenet
    mpz_class RhoExpI_h = (RhoExpI * h) % q;

    mpz_class gExp_RhoExpI;
    mpz_powm(gExp_RhoExpI.get_mpz_t(), g.get_mpz_t(), RhoExpI.get_mpz_t(), p.get_mpz_t());
    mpz_class gExp_RhoExpI_h;
    mpz_powm(gExp_RhoExpI_h.get_mpz_t(), g.get_mpz_t(), RhoExpI_h.get_mpz_t(), p.get_mpz_t());

    //Reduce exponents to mod q
    mpz_class RhoExpI_mod_q = RhoExpI % q;
    mpz_class gExp_RhoExpI_mod_q;
    mpz_powm(gExp_RhoExpI_mod_q.get_mpz_t(), g.get_mpz_t(), RhoExpI_mod_q.get_mpz_t(), p.get_mpz_t());
    mpz_class RhoExpI_h__mod_q = RhoExpI_h % q;
    mpz_class gExp_RhoExpI_h__mod_q;
    mpz_powm(gExp_RhoExpI_h__mod_q.get_mpz_t(), g.get_mpz_t(), RhoExpI_h__mod_q.get_mpz_t(), p.get_mpz_t());

    std::pair<mpz_class, mpz_class> E_Rho = ElGamal_encrypt(Rho, pk_E);//Compute E(\rho)
    std::pair<mpz_class, mpz_class> E_h = ElGamal_encrypt(h, pk_E);//Compute E(h)
    std::pair<mpz_class, mpz_class> E_RhoExpI = ElGamal_exp_ct(E_Rho, I, pk_E);//Compute E(\rho^I)
    std::pair<mpz_class, mpz_class> E_RhoExpImulh = ElGamal_mult_ct(E_RhoExpI, E_h);//Compute E(\rho^I.h)

    mpz_class decrypted_RhoExpI = ElGamal_decrypt(E_RhoExpI, sk_E);
    decrypted_RhoExpI = decrypted_RhoExpI % q;//Since this should be in exponent
    mpz_class decrypted_RhoExpImulh = ElGamal_decrypt(E_RhoExpImulh, sk_E);
    decrypted_RhoExpImulh = decrypted_RhoExpImulh % q;//Since this should be in exponent

    std::pair<mpz_class, mpz_class> E_gExp_RhoExpI_h__mod_q = ElGamal_encrypt(gExp_RhoExpI_h__mod_q, pk_E);//Compute E(g^{\rho^{I}.h mod q})
    std::pair<mpz_class, mpz_class> E_gExp_RhoExpI_hh_1_mod_q = ElGamal_exp_ct(E_gExp_RhoExpI_h__mod_q, h_1, pk_E);//Compute E(g^{\rho^{I}.h.h_1 mod q})
    std::pair<mpz_class, mpz_class> E_gExp_RhoExpI_mulh_alpha_mod_q = ElGamal_exp_ct(E_gExp_RhoExpI_h__mod_q, alpha, pk_E);//Compute E(g^{\rho^{I}.h.alpha})
    std::pair<mpz_class, mpz_class> E_gExp_RhoExpI_mulh_alpha_alpha_1_mod_q = ElGamal_exp_ct(E_gExp_RhoExpI_mulh_alpha_mod_q, alpha_1, pk_E);//Compute E(g^{\rho^{I}.h.alpha.h^{-1}})
    std::pair<mpz_class, mpz_class> E_gExp_RhoExpI_mulh_mulalpha_alpha_1_h_1_mod_q = ElGamal_exp_ct(E_gExp_RhoExpI_mulh_alpha_alpha_1_mod_q, h_1, pk_E);//Compute E(g^{\rho^{I}.h.alpha.h^{-1}})

    mpz_class decrypted_h = ElGamal_decrypt(E_h, sk_E);
    mpz_class decrypted_Rho = ElGamal_decrypt(E_Rho, sk_E);
    mpz_class decrypted_gExp_RhoExpI_hh_1_mod_q = ElGamal_decrypt(E_gExp_RhoExpI_hh_1_mod_q, sk_E);
    mpz_class decrypted_gExp_RhoExpI_mulh_alpha_alpha_1_mod_q = ElGamal_decrypt(E_gExp_RhoExpI_mulh_alpha_alpha_1_mod_q, sk_E);
    mpz_class decrypted_gExp_RhoExpI_mulh_mulalpha_alpha_1_mod_q = ElGamal_decrypt(E_gExp_RhoExpI_mulh_alpha_alpha_1_mod_q, sk_E);
    mpz_class decrypted_gExp_RhoExpI_mulh_mulalpha_alpha_1_h_1_mod_q = ElGamal_decrypt(E_gExp_RhoExpI_mulh_mulalpha_alpha_1_h_1_mod_q, sk_E);

    if (Rho != decrypted_Rho) {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Decryption failed for Rho: expected " + Rho.get_str() + ", got " + decrypted_Rho.get_str());
        //return;
    } else {
        PrintLog(LOG_LEVEL_INFO, __FILE__, __LINE__, "Rho matched with decrypted_Rho ");
    }

    if (h != decrypted_h) {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Decryption failed for h: expected " + h.get_str() + ", got " + decrypted_h.get_str());
        //return;
    } else {
        PrintLog(LOG_LEVEL_INFO, __FILE__, __LINE__, "h matched with decrypted_h ");
    }

    if (RhoExpI != decrypted_RhoExpI) {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Decryption failed for RhoExpI: expected " + RhoExpI.get_str() + ", got " + decrypted_RhoExpI.get_str());
        //return;
    } else {
        PrintLog(LOG_LEVEL_INFO, __FILE__, __LINE__, "RhoExpI matched with decrypted_RhoExpI ");
    }

    if (RhoExpI_h__mod_q != decrypted_RhoExpImulh) {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Decryption failed for RhoExpI_h: expected " + RhoExpI_h__mod_q.get_str() + ", got " + decrypted_RhoExpImulh.get_str());
        //return;
    } else {
        PrintLog(LOG_LEVEL_INFO, __FILE__, __LINE__, "RhoExpI_h matched with decrypted_RhoExpImulh ");
    }

    if (gExp_RhoExpI != gExp_RhoExpI_mod_q) {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "g^RhoExpI does not match with g^RhoExpI mod q: expected " + gExp_RhoExpI.get_str() + ", got " + gExp_RhoExpI_mod_q.get_str());
        //return;
    } else {
        PrintLog(LOG_LEVEL_INFO, __FILE__, __LINE__, "g^RhoExpI matches with g^RhoExpI mod q"); 
    }

    if (gExp_RhoExpI_h != gExp_RhoExpI_h__mod_q) {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "g^RhoExpI_h does not match with g^RhoExpI_h mod q: expected " + gExp_RhoExpI_h.get_str() + ", got " + gExp_RhoExpI_h__mod_q.get_str());
        //return;
    } else {
        PrintLog(LOG_LEVEL_INFO, __FILE__, __LINE__, "g^RhoExpI_h matches with g^RhoExpI_h mod q");
    }

    if (decrypted_gExp_RhoExpI_hh_1_mod_q != gExp_RhoExpI_mod_q) {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "g^{Rho^{I}} mod q does not match with g^{Rho^{I}.h.h^{-1}} mod q: expected " + gExp_RhoExpI.get_str() + ", got " + decrypted_gExp_RhoExpI_hh_1_mod_q.get_str());
        //return;
    } else {
        PrintLog(LOG_LEVEL_INFO, __FILE__, __LINE__, "g^{Rho^{I}} mod q matches with g^{Rho^{I}.h.h^{-1}} mod q..!!");
    }

    if (decrypted_gExp_RhoExpI_mulh_alpha_alpha_1_mod_q != gExp_RhoExpI_h__mod_q) {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "g^{Rho^{I}.h} mod q does not match with g^{Rho^{I}.h.alpha.alpha^{-1}} mod q: expected " + gExp_RhoExpI_h.get_str() + ", got " + decrypted_gExp_RhoExpI_mulh_alpha_alpha_1_mod_q.get_str());
        //return;
    } else {
        PrintLog(LOG_LEVEL_INFO, __FILE__, __LINE__, "g^{Rho^{I}.h} mod q matches with g^{Rho^{I}.h.alpha.alpha^{-1}} mod q..!!");
    }
#if 0
    mpz_class reduced_RhoExpImulh = decrypted_RhoExpImulh % q;

    mpz_class Exp_rhoExpImulh, Exp_RedrhoExpImulh, h_inverse, final_result, reduced_final_result;
    mpz_powm(Exp_rhoExpImulh.get_mpz_t(), g.get_mpz_t(), decrypted_rhoExpImulh.get_mpz_t(), p.get_mpz_t());
    mpz_powm(Exp_RedrhoExpImulh.get_mpz_t(), g.get_mpz_t(), reduced_rhoExpImulh.get_mpz_t(), p.get_mpz_t());
    mpz_invert(h_inverse.get_mpz_t(), h.get_mpz_t(), p.get_mpz_t());
    mpz_class reduced_h_inverse = h_inverse % q;
    mpz_powm(final_result.get_mpz_t(), Exp_rhoExpImulh.get_mpz_t(), h_inverse.get_mpz_t(), p.get_mpz_t());
    mpz_powm(reduced_final_result.get_mpz_t(), Exp_RedrhoExpImulh.get_mpz_t(), reduced_h_inverse.get_mpz_t(), p.get_mpz_t());

    PrintLog(LOG_LEVEL_INFO, __FILE__, __LINE__, "Decrypted messages are decrypted_rhoExpI: " + decrypted_rhoExpI.get_str() + ", decrypted_rhoExpImulh: " + decrypted_rhoExpImulh.get_str() + ", Exp_rhoExpImulh: " + Exp_rhoExpImulh.get_str() + ", Exp_RedrhoExpImulh: " + Exp_RedrhoExpImulh.get_str());
    PrintLog(LOG_LEVEL_INFO, __FILE__, __LINE__, "h_inverse: " + h_inverse.get_str() + ", reduced_h_inverse: " + reduced_h_inverse.get_str() + ", final_result: " + final_result.get_str() + ", reduced_final_result: " + reduced_final_result.get_str());
#endif
}

static void TestBlindedExponentiation1() {
    int ret;

    PrintLog(LOG_LEVEL_INFO, __FILE__, __LINE__, "El-Gamal parameters p: " + p.get_str() + " q: " + q.get_str()+ " g: " + g.get_str()+ " pk_E: " + pk_E.get_str()+ " sk_E: " + sk_E.get_str());

    Rho = rng.get_z_range(((q-1)/2)) + 1; // Randomly choose Rho in ZZ_((q-1)/2)*

    mpz_class h = rng.get_z_range(q-1) + 1;//To ensure that the element belongs to ZZ_q*
    mpz_class h_1;
    mpz_invert(h_1.get_mpz_t(), h.get_mpz_t(), q.get_mpz_t());//Compute h^-1 in ZZ_q*

    mpz_class alpha = rng.get_z_range(q-1) + 1;
    mpz_class alpha_1;
    mpz_invert(alpha_1.get_mpz_t(), alpha.get_mpz_t(), q.get_mpz_t());

    //Confirm that the inverses really work
    assert((h * h_1) % q == 1);
    assert((alpha * alpha_1) % q == 1);

    mpz_class I = rng.get_z_range(
        N-1) + 1;//I will be in the range of [1, N]
    PrintLog(LOG_LEVEL_INFO, __FILE__, __LINE__, "Chosen plain text messages are Rho: " + Rho.get_str() + " I: " + I.get_str() + " h: " + h.get_str() + " h_1: " + h_1.get_str()+ " alpha: " + alpha.get_str() + " alpha_1: " + alpha_1.get_str() );
    mpz_class RhoExpI;
    //mpz_powm(RhoExpI.get_mpz_t(), Rho.get_mpz_t(), I.get_mpz_t(), p.get_mpz_t()); previous
    //mpz_class RhoExpI_h = (RhoExpI * h) % p; previous
    mpz_powm(RhoExpI.get_mpz_t(), Rho.get_mpz_t(), I.get_mpz_t(), q.get_mpz_t());//Since this is an exponenet
    //mpz_class RhoExpI_h = (RhoExpI * h) % q;

    mpz_class gExp_RhoExpI;
    mpz_powm(gExp_RhoExpI.get_mpz_t(), g.get_mpz_t(), RhoExpI.get_mpz_t(), p.get_mpz_t());

    //Server beta broadcasts E_Rho
    std::pair<mpz_class, mpz_class> E_Rho = ElGamal_encrypt(Rho, pk_E);//Compute E(\rho)

    //Client computes E_RhoExpI
    std::pair<mpz_class, mpz_class> E_RhoExpI = ElGamal_exp_ct(E_Rho, I, pk_E);//Compute E(\rho^I)
    //Client computes E_RhoExpI_h and sends to server beta
    std::pair<mpz_class, mpz_class> E_h = ElGamal_encrypt(h, pk_E);//Compute E(h)
    std::pair<mpz_class, mpz_class> E_RhoExpI_h = ElGamal_mult_ct(E_RhoExpI, E_h);//Compute E(\rho^I.h)

    //Server beta decrypts and reduces to mod q
    mpz_class decrypted_RhoExpI_h = ElGamal_decrypt(E_RhoExpI_h, sk_E);
    mpz_class decrypted_RhoExpI_h__modq = decrypted_RhoExpI_h % q;//Reduce to mod q
    //Server beta raises to g
    mpz_class gExp_RhoExpI_h__modq;
    mpz_powm(gExp_RhoExpI_h__modq.get_mpz_t(), g.get_mpz_t(), decrypted_RhoExpI_h__modq.get_mpz_t(), p.get_mpz_t());
    std::pair<mpz_class, mpz_class> E_gExp_RhoExpI_h__modq = ElGamal_encrypt(gExp_RhoExpI_h__modq, pk_E);//Sends corresponding ciphertext

    //Server alpha removes h by raising it to h^{-1}
    std::pair<mpz_class, mpz_class> E_gExp_RhoExpI_hh_1__modq = ElGamal_exp_ct(E_gExp_RhoExpI_h__modq, h_1, pk_E);

    //If decrypted must match with gExp_RhoExpI
    mpz_class decrypted_gExp_RhoExpI_hh_1__modq = ElGamal_decrypt(E_gExp_RhoExpI_hh_1__modq, sk_E);

    assert(decrypted_gExp_RhoExpI_hh_1__modq == gExp_RhoExpI);

    if (decrypted_gExp_RhoExpI_hh_1__modq == gExp_RhoExpI) {
        PrintLog(LOG_LEVEL_INFO, __FILE__, __LINE__, "Server Beta: Decryption successful and matches with gExp_RhoExpI");
    } else {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Server Beta: Decryption failed or does not match with gExp_RhoExpI. Expected :" + gExp_RhoExpI.get_str() + " but got: " + decrypted_gExp_RhoExpI_hh_1__modq.get_str());
    }
}

static void TestBlindedExponentiation2() {
    int ret;

    PrintLog(LOG_LEVEL_INFO, __FILE__, __LINE__, "El-Gamal parameters p: " + p.get_str() + " q: " + q.get_str()+ " g: " + g.get_str()+ " pk_E: " + pk_E.get_str()+ " sk_E: " + sk_E.get_str());
    //Choose random message and random exponent
    mpz_class m1 = ElGamal_randomGroupElement();
    mpz_class m2 = ElGamal_randomGroupElement();
    mpz_class exp = rng.get_z_range(q);
    PrintLog(LOG_LEVEL_INFO, __FILE__, __LINE__, "Chosen plaintext messages are m1: " + m1.get_str() + " m2: " + m2.get_str() + " exp: " + exp.get_str());

    mpz_class m3 = (m1*m2)%p;
    mpz_class m4;
    mpz_powm(m4.get_mpz_t(), m1.get_mpz_t(), exp.get_mpz_t(), p.get_mpz_t());

    mpz_class m5 = (m4*m2)%p;//(m1^exp)*m2

    std::pair<mpz_class, mpz_class> E_m1 = ElGamal_encrypt(m1, pk_E);
    std::pair<mpz_class, mpz_class> E_m2 = ElGamal_encrypt(m2, pk_E);
    std::pair<mpz_class, mpz_class> E_m3 = ElGamal_mult_ct(E_m1, E_m2);
    std::pair<mpz_class, mpz_class> E_m4 = ElGamal_exp_ct(E_m1, exp, pk_E);
    std::pair<mpz_class, mpz_class> E_m5 = ElGamal_mult_ct(E_m4, E_m2);

    mpz_class decrypted_m1 = ElGamal_decrypt(E_m1, sk_E);
    mpz_class decrypted_m3 = ElGamal_decrypt(E_m3, sk_E);
    mpz_class decrypted_m4 = ElGamal_decrypt(E_m4, sk_E);
    mpz_class decrypted_m5 = ElGamal_decrypt(E_m5, sk_E);

    if (decrypted_m1 == m1) {
        PrintLog(LOG_LEVEL_INFO, __FILE__, __LINE__, "El-Gamal encryption works");
    } else {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "El-Gamal encryption is not working. Expected: " + m1.get_str() + " but got: " + decrypted_m1.get_str());
    }

    if (decrypted_m3 == m3) {
        PrintLog(LOG_LEVEL_INFO, __FILE__, __LINE__, "El-Gamal multiplication works");
    } else {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "El-Gamal multiplication is not working. Expected: " + m3.get_str() + " but got: " + decrypted_m3.get_str());
    }

    if (decrypted_m4 == m4) {
        PrintLog(LOG_LEVEL_INFO, __FILE__, __LINE__, "El-Gamal exponentiation works");
    } else {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "El-Gamal exponentiation is not working. Expected: " + m4.get_str() + " but got: " + decrypted_m4.get_str());
    }

    if (decrypted_m5 == m5) {
        PrintLog(LOG_LEVEL_INFO, __FILE__, __LINE__, "El-Gamal multiplication after exponentiation works");
    } else {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "El-Gamal multiplication after exponentiation is not working. Expected: " + m5.get_str() + " but got: " + decrypted_m5.get_str());
    }

    PrintLog(LOG_LEVEL_INFO, __FILE__, __LINE__, "El-Gamal in ZZ*_q parameters g_q: " + g_q.get_str()+ " pk_E_q: " + pk_E_q.get_str()+ " sk_E_q: " + sk_E_q.get_str());

    Rho = rng.get_z_range(q-1)+1;//i.e., within ZZ_q*
    mpz_class h = rng.get_z_range(q-1)+1;//i.e., within ZZ_q*
    mpz_class h_1;
    mpz_invert(h_1.get_mpz_t(), h.get_mpz_t(), q.get_mpz_t());

    mpz_class alpha = rng.get_z_range(q-1)+1;//i.e., within ZZ_q*
    mpz_class alpha_1;
    mpz_invert(alpha_1.get_mpz_t(), alpha.get_mpz_t(), q.get_mpz_t());

    mpz_class I = rng.get_z_range(q-1)+1;//i.e., within ZZ_q*
    PrintLog(LOG_LEVEL_INFO, __FILE__, __LINE__, "Chosen plaintext messages are Rho: " + Rho.get_str() + " h: " + h.get_str()+ " h_1: " + h_1.get_str() + " I: " + I.get_str() + " alpha: " + alpha.get_str() + " alpha_1: " + alpha_1.get_str());

    mpz_class Rho_h = (Rho*h)%q;
    mpz_class Rho_pow_I;
    mpz_powm(Rho_pow_I.get_mpz_t(), Rho.get_mpz_t(), I.get_mpz_t(), q.get_mpz_t());
    mpz_class g_pow_Rho_pow_I;
    mpz_powm(g_pow_Rho_pow_I.get_mpz_t(), g.get_mpz_t(), Rho_pow_I.get_mpz_t(), p.get_mpz_t());

    mpz_class Rho_pow_I__h = (Rho_pow_I*h)%q;//(Rho^I)*h

    std::pair<mpz_class, mpz_class> E_Rho = ElGamal_q_encrypt(Rho, pk_E_q);
    std::pair<mpz_class, mpz_class> E_h = ElGamal_q_encrypt(h, pk_E_q);
    std::pair<mpz_class, mpz_class> E_Rho_h = ElGamal_q_mult_ct(E_Rho, E_h);
    std::pair<mpz_class, mpz_class> E_Rho_pow_I = ElGamal_q_exp_ct(E_Rho, I, pk_E_q);
    std::pair<mpz_class, mpz_class> E_Rho_pow_I__h = ElGamal_q_mult_ct(E_Rho_pow_I, E_h);

    mpz_class decrypted_Rho = ElGamal_q_decrypt(E_Rho, sk_E_q);
    mpz_class decrypted_Rho_h = ElGamal_q_decrypt(E_Rho_h, sk_E_q);
    mpz_class decrypted_Rho_pow_I = ElGamal_q_decrypt(E_Rho_pow_I, sk_E_q);
    mpz_class decrypted_Rho_pow_I__h = ElGamal_q_decrypt(E_Rho_pow_I__h, sk_E_q);
   
    //Server beta decrypts and perform g^{decrypted_Rho_pow_I__h} mod p
    mpz_class g_pow_Rho_pow_I__h;
    mpz_powm(g_pow_Rho_pow_I__h.get_mpz_t(), g.get_mpz_t(), decrypted_Rho_pow_I__h.get_mpz_t(), p.get_mpz_t());
    //Encrypts and sends that under ElGamal encryption in GG
    std::pair<mpz_class, mpz_class> E_g_pow_Rho_pow_I__h = ElGamal_encrypt(g_pow_Rho_pow_I__h, pk_E);

    //Server alpha semi-homomorphically raises that to alpha under ElGamal encryption in GG
    std::pair<mpz_class, mpz_class> E_g_pow_Rho_pow_I__h_alpha = ElGamal_exp_ct(E_g_pow_Rho_pow_I__h, alpha, pk_E);

    //Client semi-homomorphically raises that to h_1 under ElGamal encryption in GG
    std::pair<mpz_class, mpz_class> E_g_pow_Rho_pow_I__h_alpha_h_1 = ElGamal_exp_ct(E_g_pow_Rho_pow_I__h_alpha, h_1, pk_E);

    //Server semi-homomorphically raises that to alpha_1 under ElGamal encryption in GG
    std::pair<mpz_class, mpz_class> E_g_pow_Rho_pow_I__h_alpha_h_1_alpha_1 = ElGamal_exp_ct(E_g_pow_Rho_pow_I__h_alpha_h_1, alpha_1, pk_E);

    //Server beta decrypts that under ElGamal encryption in GG
    mpz_class decrypted_g_pow_Rho_pow_I__h_alpha_h_1_alpha_1 = ElGamal_decrypt(E_g_pow_Rho_pow_I__h_alpha_h_1_alpha_1, sk_E);

    if (decrypted_Rho == Rho) {
        PrintLog(LOG_LEVEL_INFO, __FILE__, __LINE__, "El-Gamal encryption works in ZZ*_q");
    } else {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "El-Gamal encryption is not working in ZZ*_q. Expected: " + Rho.get_str() + " but got: " + decrypted_Rho.get_str());
    }

    if (decrypted_Rho_h == Rho_h) {
        PrintLog(LOG_LEVEL_INFO, __FILE__, __LINE__, "El-Gamal multiplication works in ZZ*_q");
    } else {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "El-Gamal multiplication is not working in ZZ*_q. Expected: " + Rho_h.get_str() + " but got: " + decrypted_Rho_h.get_str());
    }

    if (decrypted_Rho_pow_I == Rho_pow_I) {
        PrintLog(LOG_LEVEL_INFO, __FILE__, __LINE__, "El-Gamal exponentiation works in ZZ*_q");
    } else {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "El-Gamal exponentiation is not working in ZZ*_q. Expected: " + Rho_pow_I.get_str() + " but got: " + decrypted_Rho_pow_I.get_str());
    }

    if (decrypted_Rho_pow_I__h == Rho_pow_I__h) {
        PrintLog(LOG_LEVEL_INFO, __FILE__, __LINE__, "El-Gamal multiplication after exponentiation works in ZZ*_q");
    } else {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "El-Gamal multiplication after exponentiation is not working in ZZ*_q. Expected: " + Rho_pow_I__h.get_str() + " but got: " + decrypted_Rho_pow_I__h.get_str());
    }

    if (decrypted_g_pow_Rho_pow_I__h_alpha_h_1_alpha_1 == g_pow_Rho_pow_I) {
        PrintLog(LOG_LEVEL_INFO, __FILE__, __LINE__, "Finding tag homomorphically works..!!..Ha ha..Thank you..:) :) :) :) :)");
    } else {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Finding tag homomorphically is not working :( Expected: " + g_pow_Rho_pow_I.get_str() + " but got: " + decrypted_g_pow_Rho_pow_I__h_alpha_h_1_alpha_1.get_str());
    }
}

// Test function for FHE_Enc_SDBElement and FHE_Dec_SDBElement
static void Test_FHE_DBElement() {
    /* First of all retrieve all the one-time initialized materials from the saved location */
    p = import_from_file_to_mpz_class(ONE_TIME_MATERIALS_LOCATION_BETA + "p.bin");
    q = import_from_file_to_mpz_class(ONE_TIME_MATERIALS_LOCATION_BETA + "q.bin");
    g = import_from_file_to_mpz_class(ONE_TIME_MATERIALS_LOCATION_BETA + "g.bin");
    g_q = import_from_file_to_mpz_class(ONE_TIME_MATERIALS_LOCATION_BETA + "g_q.bin");
    r = import_from_file_to_mpz_class(ONE_TIME_MATERIALS_LOCATION_BETA + "r.bin");
    pk_E = import_from_file_to_mpz_class(ONE_TIME_MATERIALS_LOCATION_BETA + "pk_E.bin");
    sk_E = import_from_file_to_mpz_class(ONE_TIME_MATERIALS_LOCATION_BETA + "sk_E.bin");
    pk_E_q = import_from_file_to_mpz_class(ONE_TIME_MATERIALS_LOCATION_BETA + "pk_E_q.bin");
    sk_E_q = import_from_file_to_mpz_class(ONE_TIME_MATERIALS_LOCATION_BETA + "sk_E_q.bin");
    Serial::DeserializeFromFile(ONE_TIME_MATERIALS_LOCATION_BETA + "FHEcryptoContext.bin", FHEcryptoContext, SerType::BINARY);
    Serial::DeserializeFromFile(ONE_TIME_MATERIALS_LOCATION_BETA + "pk_F.bin", pk_F, SerType::BINARY);
    Serial::DeserializeFromFile(ONE_TIME_MATERIALS_LOCATION_BETA + "sk_F.bin", sk_F, SerType::BINARY);
    Serial::DeserializeFromFile(ONE_TIME_MATERIALS_LOCATION_BETA + "vectorOnesforElement_ct.bin", vectorOnesforElement_ct, SerType::BINARY);
    Serial::DeserializeFromFile(ONE_TIME_MATERIALS_LOCATION_BETA + "vectorOnesforTag_ct.bin", vectorOnesforTag_ct, SerType::BINARY);

    PrintLog(LOG_LEVEL_TRACE, __FILE__, __LINE__, "Server Beta: Loaded one-time initialization materials");

    // Generate random block_content of PLAINTEXT_PIR_BLOCK_DATA_SIZE bits
    mpz_class block_1_content = rng.get_z_bits(PLAINTEXT_PIR_BLOCK_DATA_SIZE);
    // Generate random block_index of log_N bits
    mpz_class block_1_index = rng.get_z_bits(log_N);
    // Generate random block_content of PLAINTEXT_PIR_BLOCK_DATA_SIZE bits
    mpz_class block_2_content = rng.get_z_bits(PLAINTEXT_PIR_BLOCK_DATA_SIZE);
    // Generate random block_index of log_N bits
    mpz_class block_2_index = rng.get_z_bits(log_N);
    // Generate random tag of P_BITS bits
    mpz_class tag_1 = rng.get_z_bits(P_BITS);
    // Generate another random tag of P_BITS bits
    mpz_class tag_2 = rng.get_z_bits(P_BITS);

    //PrintLog(LOG_LEVEL_INFO, __FILE__, __LINE__, "Chosen block 1 content: " + block_1_content.get_str() + "\nblock 1 index: " + block_1_index.get_str() + "\nblock 2 content: " + block_2_content.get_str() + "\nblock 2 index: " + block_2_index.get_str() + "\ntag 1: " + tag_1.get_str()+ "\ntag 2: " + tag_2.get_str());

    // Encrypt
    Ciphertext<DCRTPoly> ct_element_1 = FHE_Enc_SDBElement((block_1_content << log_N) | block_1_index);
    Ciphertext<DCRTPoly> ct_element_2 = FHE_Enc_SDBElement((block_2_content << log_N) | block_2_index);
    Ciphertext<DCRTPoly> ct_tag_1 = FHE_Enc_Tag(tag_1);
    Ciphertext<DCRTPoly> ct_tag_2 = FHE_Enc_Tag(tag_2);
    
    Ciphertext<DCRTPoly> selectElementBits_ct = FHE_Enc_SDBElement(mpz_class(0));
    Ciphertext<DCRTPoly> selectTagBits_ct = FHE_Enc_Tag(mpz_class(0));
    // Also tested with 1, which is 0b...000000000000001000000000000001
    //Ciphertext<DCRTPoly> selectElementBits_ct = vectorOnesforElement_ct;
    //Ciphertext<DCRTPoly> selectTagBits_ct = vectorOnesforTag_ct;
    //Ciphertext<DCRTPoly> selectTagBits_ct = vectorOnesforTag_ct;


    // Decrypt
    mpz_class dec_block_1_content, dec_block_1_index, dec_block_2_content, dec_block_2_index, dec_tag_1, dec_tag_2, dec_fnd, dec_selected_tag, dec_selected_content, dec_selected_index, dec_content_and_index;
    FHE_Dec_SDBElement(ct_element_1, dec_content_and_index);
    dec_block_1_content = (dec_content_and_index >> log_N);
    dec_block_1_index = (dec_content_and_index & ((1U << log_N) - 1U)); 

    FHE_Dec_SDBElement(ct_element_2, dec_content_and_index);
    dec_block_2_content = (dec_content_and_index >> log_N);
    dec_block_2_index = (dec_content_and_index & ((1U << log_N) - 1U)); 
    
    FHE_Dec_Tag(ct_tag_1, dec_tag_1);
    FHE_Dec_Tag(ct_tag_2, dec_tag_2);
    FHE_Dec_Tag(selectTagBits_ct, dec_fnd);
    Ciphertext<DCRTPoly> ct_selected_element = FHE_SelectElement(selectElementBits_ct, ct_element_1, ct_element_2);
    Ciphertext<DCRTPoly> ct_selected_tag = FHE_SelectTag(selectTagBits_ct, ct_tag_1, ct_tag_2);
    FHE_Dec_Tag(ct_selected_tag, dec_selected_tag);
    FHE_Dec_SDBElement(ct_selected_element, dec_content_and_index);
    dec_selected_content = (dec_content_and_index >> log_N);
    dec_selected_index = (dec_content_and_index & ((1U << log_N) - 1U)); 

    // Verify
    if (block_1_content != dec_block_1_content) {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "block_content mismatch!\nExpected: " + block_1_content.get_str() + "\nObtained: " + dec_block_1_content.get_str());
    } else {
        PrintLog(LOG_LEVEL_INFO, __FILE__, __LINE__, "block_content matches.");
    }

    if (block_1_index != dec_block_1_index) {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "block_index mismatch!\nExpected: " + block_1_index.get_str() + "\nObtained: " + dec_block_1_index.get_str());
    } else {
        PrintLog(LOG_LEVEL_INFO, __FILE__, __LINE__, "block_index matches.");
    }

    if (tag_1 != dec_tag_1) {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "tag_1 mismatch!\nExpected: " + tag_1.get_str() + "\nObtained: " + dec_tag_1.get_str());
    } else {
        PrintLog(LOG_LEVEL_INFO, __FILE__, __LINE__, "tag_1 matches.");
    }

    if (tag_2 != dec_tag_2) {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "tag_2 mismatch!\nExpected: " + tag_2.get_str() + "\nObtained: " + dec_tag_2.get_str());
    } else {
        PrintLog(LOG_LEVEL_INFO, __FILE__, __LINE__, "tag_2 matches.");
    }

    if (dec_fnd != 0) {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "fnd mismatch!\nExpected: 0 but got: " + dec_fnd.get_str());
    } else {
        PrintLog(LOG_LEVEL_INFO, __FILE__, __LINE__, "fnd matches");
    }

    if (dec_selected_tag == tag_1) {
        PrintLog(LOG_LEVEL_INFO, __FILE__, __LINE__, "Homomorphic selection works for the tags..!!..Ha ha..Thank you..:) :) :) :) :)");
    } else {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Homomorphic selection is not working for the tags:( Expected: " + tag_1.get_str() + " but got: " + dec_selected_tag.get_str());
    }

    if (dec_selected_content == block_1_content) {
        PrintLog(LOG_LEVEL_INFO, __FILE__, __LINE__, "Homomorphic selection works for the block content..!!..Ha ha..Thank you..:) :) :) :) :)");
    } else {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Homomorphic selection is not working for the block content :( Expected: " + block_1_content.get_str() + " but got: " + dec_selected_content.get_str());
    }

    if (dec_selected_index == block_1_index) {
        PrintLog(LOG_LEVEL_INFO, __FILE__, __LINE__, "Homomorphic selection works for the block index..!!..Ha ha..Thank you..:) :) :) :) :)");
    } else {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Homomorphic selection is not working for the block index :( Expected: " + block_1_index.get_str() + " but got: " + dec_selected_index.get_str());
    }
}


static void TestPKEOperations_beta() {
    int ret = -1;
    size_t received_sz = 0;
    int ret_recv = 0;
    mpz_class m1_local, m2_local, m3_local, m4_local, c11_local, c12_local, c21_local, c22_local, c31_local, c32_local, c41_local, c42_local, tag_local;
    Ciphertext<DCRTPoly> ct_tag_local;

    // Receive from server alpha

    // Receive m1
    ret_recv = recvAll(sock_beta_alpha_con, net_buf, sizeof(net_buf), &received_sz);
    if (ret_recv != 0) {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Server Beta: Failed to receive m1 from Server Alpha");
        close(sock_beta_alpha_con);
        return;
    }
    m1_local = mpz_class(std::string(net_buf, received_sz));

    // Receive m2
    ret_recv = recvAll(sock_beta_alpha_con, net_buf, sizeof(net_buf), &received_sz);
    if (ret_recv != 0) {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Server Beta: Failed to receive m2 from Server Alpha");
        close(sock_beta_alpha_con);
        return;
    }
    m2_local = mpz_class(std::string(net_buf, received_sz));

    // Receive m3
    ret_recv = recvAll(sock_beta_alpha_con, net_buf, sizeof(net_buf), &received_sz);
    if (ret_recv != 0) {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Server Beta: Failed to receive m3 from Server Alpha");
        close(sock_beta_alpha_con);
        return;
    }
    m3_local = mpz_class(std::string(net_buf, received_sz));

    // Receive m4
    ret_recv = recvAll(sock_beta_alpha_con, net_buf, sizeof(net_buf), &received_sz);
    if (ret_recv != 0) {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Server Beta: Failed to receive m4 from Server Alpha");
        close(sock_beta_alpha_con);
        return;
    }
    m4_local = mpz_class(std::string(net_buf, received_sz));

    // Receive c11
    ret_recv = recvAll(sock_beta_alpha_con, net_buf, sizeof(net_buf), &received_sz);
    if (ret_recv != 0) {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Server Beta: Failed to receive c11 from Server Alpha");
        close(sock_beta_alpha_con);
        return;
    }
    c11_local = mpz_class(std::string(net_buf, received_sz));

    // Receive c12
    ret_recv = recvAll(sock_beta_alpha_con, net_buf, sizeof(net_buf), &received_sz);
    if (ret_recv != 0) {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Server Beta: Failed to receive c12 from Server Alpha");
        close(sock_beta_alpha_con);
        return;
    }
    c12_local = mpz_class(std::string(net_buf, received_sz));

    // Receive c21
    ret_recv = recvAll(sock_beta_alpha_con, net_buf, sizeof(net_buf), &received_sz);
    if (ret_recv != 0) {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Server Beta: Failed to receive c21 from Server Alpha");
        close(sock_beta_alpha_con);
        return;
    }
    c21_local = mpz_class(std::string(net_buf, received_sz));

    // Receive c22
    ret_recv = recvAll(sock_beta_alpha_con, net_buf, sizeof(net_buf), &received_sz);
    if (ret_recv != 0) {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Server Beta: Failed to receive c22 from Server Alpha");
        close(sock_beta_alpha_con);
        return;
    }
    c22_local = mpz_class(std::string(net_buf, received_sz));

    // Receive c31
    ret_recv = recvAll(sock_beta_alpha_con, net_buf, sizeof(net_buf), &received_sz);
    if (ret_recv != 0) {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Server Beta: Failed to receive c31 from Server Alpha");
        close(sock_beta_alpha_con);
        return;
    }
    c31_local = mpz_class(std::string(net_buf, received_sz));

    // Receive c32
    ret_recv = recvAll(sock_beta_alpha_con, net_buf, sizeof(net_buf), &received_sz);
    if (ret_recv != 0) {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Server Beta: Failed to receive c32 from Server Alpha");
        close(sock_beta_alpha_con);
        return;
    }
    c32_local = mpz_class(std::string(net_buf, received_sz));

    // Receive c41
    ret_recv = recvAll(sock_beta_alpha_con, net_buf, sizeof(net_buf), &received_sz);
    if (ret_recv != 0) {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Server Beta: Failed to receive c41 from Server Alpha");
        close(sock_beta_alpha_con);
        return;
    }
    c41_local = mpz_class(std::string(net_buf, received_sz));

    // Receive c42
    ret_recv = recvAll(sock_beta_alpha_con, net_buf, sizeof(net_buf), &received_sz);
    if (ret_recv != 0) {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Server Beta: Failed to receive c42 from Server Alpha");
        close(sock_beta_alpha_con);
        return;
    }
    c42_local = mpz_class(std::string(net_buf, received_sz));

    // Receive tag_local
    ret_recv = recvAll(sock_beta_alpha_con, net_buf, sizeof(net_buf), &received_sz);
    if (ret_recv != 0) {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Server Beta: Failed to receive tag_local from Server Alpha");
        close(sock_beta_alpha_con);
        return;
    }
    tag_local = mpz_class(std::string(net_buf, received_sz));

    // Receive ct_tag_local
    ret_recv = recvAll(sock_beta_alpha_con, net_buf, sizeof(net_buf), &received_sz);
    if (ret_recv != 0)
    {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Failed to receive ct_tag_local from Server Beta");
        close(sock_beta_alpha_con);
        return;
    }
    Serial::DeserializeFromString(ct_tag_local, std::string(net_buf, received_sz));

    // Decrypt and check
    mpz_class decrypted_m1 = ElGamal_decrypt({c11_local, c12_local}, sk_E);
    mpz_class decrypted_m2 = ElGamal_decrypt({c21_local, c22_local}, sk_E);
    mpz_class decrypted_m3 = ElGamal_decrypt({c31_local, c32_local}, sk_E);
    mpz_class decrypted_m4 = ElGamal_decrypt({c41_local, c42_local}, sk_E);

    if (decrypted_m1 != m1_local) {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Server Beta: Decrypted m1: " + decrypted_m1.get_str() + " does not match with expected value: " + m1_local.get_str() + " !!");
    } else {
        PrintLog(LOG_LEVEL_INFO, __FILE__, __LINE__, "Server Beta: Decrypted m1 matches with expected value");
    }

    if (decrypted_m2 != m2_local) {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Server Beta: Decrypted m2: " + decrypted_m2.get_str() + " does not match with expected value: " + m2_local.get_str() + " !!");
    } else {
        PrintLog(LOG_LEVEL_INFO, __FILE__, __LINE__, "Server Beta: Decrypted m2 matches with expected value");
    }

    if (decrypted_m3 != m3_local) {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Server Beta: Decrypted m3: " + decrypted_m3.get_str() + " does not match with expected value: " + m3_local.get_str() + " !!");
    } else {
        PrintLog(LOG_LEVEL_INFO, __FILE__, __LINE__, "Server Beta: Decrypted m3 matches with expected value");
    }

    if (decrypted_m4 != m4_local) {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Server Beta: Decrypted m4: " + decrypted_m4.get_str() + " does not match with expected value: " + m4_local.get_str() + " !!");
    } else {
        PrintLog(LOG_LEVEL_INFO, __FILE__, __LINE__, "Server Beta: Decrypted m4 matches with expected value");
    }

    mpz_class decrypted_tag;
    FHE_Dec_Tag(ct_tag_local, decrypted_tag);

    if (decrypted_tag != tag_local) {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Server Beta: Decrypted tag: " + decrypted_tag.get_str() + " does not match with expected value: " + tag_local.get_str() + " !!");
    } else {
        PrintLog(LOG_LEVEL_INFO, __FILE__, __LINE__, "Server Beta: Decrypted tag matches with expected value");
    }



    // Receive from server gamma
    // Receive m1
    ret_recv = recvAll(sock_beta_gamma_con, net_buf, sizeof(net_buf), &received_sz);
    if (ret_recv != 0) {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Server Beta: Failed to receive m1 from Server gamma");
        close(sock_beta_gamma_con);
        return;
    }
    m1_local = mpz_class(std::string(net_buf, received_sz));

    // Receive m2
    ret_recv = recvAll(sock_beta_gamma_con, net_buf, sizeof(net_buf), &received_sz);
    if (ret_recv != 0) {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Server Beta: Failed to receive m2 from Server gamma");
        close(sock_beta_gamma_con);
        return;
    }
    m2_local = mpz_class(std::string(net_buf, received_sz));

    // Receive m3
    ret_recv = recvAll(sock_beta_gamma_con, net_buf, sizeof(net_buf), &received_sz);
    if (ret_recv != 0) {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Server Beta: Failed to receive m3 from Server gamma");
        close(sock_beta_gamma_con);
        return;
    }
    m3_local = mpz_class(std::string(net_buf, received_sz));

    // Receive m4
    ret_recv = recvAll(sock_beta_gamma_con, net_buf, sizeof(net_buf), &received_sz);
    if (ret_recv != 0) {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Server Beta: Failed to receive m4 from Server gamma");
        close(sock_beta_gamma_con);
        return;
    }
    m4_local = mpz_class(std::string(net_buf, received_sz));

    // Receive c11
    ret_recv = recvAll(sock_beta_gamma_con, net_buf, sizeof(net_buf), &received_sz);
    if (ret_recv != 0) {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Server Beta: Failed to receive c11 from Server gamma");
        close(sock_beta_gamma_con);
        return;
    }
    c11_local = mpz_class(std::string(net_buf, received_sz));

    // Receive c12
    ret_recv = recvAll(sock_beta_gamma_con, net_buf, sizeof(net_buf), &received_sz);
    if (ret_recv != 0) {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Server Beta: Failed to receive c12 from Server gamma");
        close(sock_beta_gamma_con);
        return;
    }
    c12_local = mpz_class(std::string(net_buf, received_sz));

    // Receive c21
    ret_recv = recvAll(sock_beta_gamma_con, net_buf, sizeof(net_buf), &received_sz);
    if (ret_recv != 0) {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Server Beta: Failed to receive c21 from Server gamma");
        close(sock_beta_gamma_con);
        return;
    }
    c21_local = mpz_class(std::string(net_buf, received_sz));

    // Receive c22
    ret_recv = recvAll(sock_beta_gamma_con, net_buf, sizeof(net_buf), &received_sz);
    if (ret_recv != 0) {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Server Beta: Failed to receive c22 from Server gamma");
        close(sock_beta_gamma_con);
        return;
    }
    c22_local = mpz_class(std::string(net_buf, received_sz));

    // Receive c31
    ret_recv = recvAll(sock_beta_gamma_con, net_buf, sizeof(net_buf), &received_sz);
    if (ret_recv != 0) {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Server Beta: Failed to receive c31 from Server gamma");
        close(sock_beta_gamma_con);
        return;
    }
    c31_local = mpz_class(std::string(net_buf, received_sz));

    // Receive c32
    ret_recv = recvAll(sock_beta_gamma_con, net_buf, sizeof(net_buf), &received_sz);
    if (ret_recv != 0) {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Server Beta: Failed to receive c32 from Server gamma");
        close(sock_beta_gamma_con);
        return;
    }
    c32_local = mpz_class(std::string(net_buf, received_sz));

    // Receive c41
    ret_recv = recvAll(sock_beta_gamma_con, net_buf, sizeof(net_buf), &received_sz);
    if (ret_recv != 0) {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Server Beta: Failed to receive c41 from Server gamma");
        close(sock_beta_gamma_con);
        return;
    }
    c41_local = mpz_class(std::string(net_buf, received_sz));

    // Receive c42
    ret_recv = recvAll(sock_beta_gamma_con, net_buf, sizeof(net_buf), &received_sz);
    if (ret_recv != 0) {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Server Beta: Failed to receive c42 from Server gamma");
        close(sock_beta_gamma_con);
        return;
    }
    c42_local = mpz_class(std::string(net_buf, received_sz));

    // Receive tag_local
    ret_recv = recvAll(sock_beta_gamma_con, net_buf, sizeof(net_buf), &received_sz);
    if (ret_recv != 0) {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Server Beta: Failed to receive tag_local from Server gamma");
        close(sock_beta_gamma_con);
        return;
    }
    tag_local = mpz_class(std::string(net_buf, received_sz));

    // Receive ct_tag_local
    ret_recv = recvAll(sock_beta_gamma_con, net_buf, sizeof(net_buf), &received_sz);
    if (ret_recv != 0)
    {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Failed to receive ct_tag_local from Server Beta");
        close(sock_beta_gamma_con);
        return;
    }
    Serial::DeserializeFromString(ct_tag_local, std::string(net_buf, received_sz));

    // Decrypt and check
    decrypted_m1 = ElGamal_decrypt({c11_local, c12_local}, sk_E);
    decrypted_m2 = ElGamal_decrypt({c21_local, c22_local}, sk_E);
    decrypted_m3 = ElGamal_decrypt({c31_local, c32_local}, sk_E);
    decrypted_m4 = ElGamal_decrypt({c41_local, c42_local}, sk_E);

    if (decrypted_m1 != m1_local) {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Server Beta: Decrypted m1: " + decrypted_m1.get_str() + " does not match with expected value: " + m1_local.get_str() + " !!");
    } else {
        PrintLog(LOG_LEVEL_INFO, __FILE__, __LINE__, "Server Beta: Decrypted m1 matches with expected value");
    }

    if (decrypted_m2 != m2_local) {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Server Beta: Decrypted m2: " + decrypted_m2.get_str() + " does not match with expected value: " + m2_local.get_str() + " !!");
    } else {
        PrintLog(LOG_LEVEL_INFO, __FILE__, __LINE__, "Server Beta: Decrypted m2 matches with expected value");
    }

    if (decrypted_m3 != m3_local) {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Server Beta: Decrypted m3: " + decrypted_m3.get_str() + " does not match with expected value: " + m3_local.get_str() + " !!");
    } else {
        PrintLog(LOG_LEVEL_INFO, __FILE__, __LINE__, "Server Beta: Decrypted m3 matches with expected value");
    }

    if (decrypted_m4 != m4_local) {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Server Beta: Decrypted m4: " + decrypted_m4.get_str() + " does not match with expected value: " + m4_local.get_str() + " !!");
    } else {
        PrintLog(LOG_LEVEL_INFO, __FILE__, __LINE__, "Server Beta: Decrypted m4 matches with expected value");
    }

    FHE_Dec_Tag(ct_tag_local, decrypted_tag);

    if (decrypted_tag != tag_local) {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Server Beta: Decrypted tag: " + decrypted_tag.get_str() + " does not match with expected value: " + tag_local.get_str() + " !!");
    } else {
        PrintLog(LOG_LEVEL_INFO, __FILE__, __LINE__, "Server Beta: Decrypted tag matches with expected value");
    }

    return;
}

static void TestSrv_beta()
{
    //TestPKEOperations_beta();
    //TestSelShuffDBSearchTag_beta();
    //TestShelterDPFSearch_beta();
    //TestClientProcessing_beta();
    TestShuffDBFetch_beta();
    //Test_FHE_DBElement();
}

#if TEST_SHUFF_DB_FETCH
#warning this function will not work, due to the absense of this large global array
static void TestShuffDBFetch_beta(){
    int ret = -1;
    u_int test_failure = 0;
    plain_db_entry read_entry;
    mpz_t tmp;
    mpz_init(tmp);

    #define NUM_FETCH_TEST_CNT 4
    /* First perform the per-epoch initialization */
    PerEpochOperations_beta();

    // RNG: mt19937 seeded from random_device
    std::random_device rd;
    std::mt19937 gen(rd());

    /* I: Choose a random index to be fetched from the shuffled databases */
    std::uniform_int_distribution<std::uint64_t> dist(1, (N+sqrt_N)); // Indices are from 1 to N+sqrt(N)
    std::uint64_t I;

    pdb.open(pdb_filename, std::ios::in | std::ios::out | std::ios::binary | std::ios::app);
    D_alpha.open(D_alpha_filename, std::ios::in | std::ios::binary | std::ios::app);
    D_gamma.open(D_gamma_filename, std::ios::in | std::ios::binary | std::ios::app);
   
    /* Figureout its corresponding location in the shuffled database */
    uint64_t loc = 0;

    for (unsigned int fetch_trial; fetch_trial < NUM_FETCH_TEST_CNT; fetch_trial++){
        I = dist(gen);


        loc = TMP_IDX_LOC_MAP[(I - 1)]; // Location for item I is stored at index (i-1)

        /* Fetch the share from that location of D_alpha */
        D_alpha.seekg((loc * NUM_BYTES_PER_SDB_ELEMENT), std::ios::beg);
        D_alpha.read(reinterpret_cast<char *>(net_buf), NUM_BYTES_PER_SDB_ELEMENT);
        mpz_import(tmp, NUM_BYTES_PER_SDB_ELEMENT, 1, 1, 1, 0, net_buf);
        mpz_class d_alpha = mpz_class(tmp);

        /* Fetch the share from that location of D_gamma */
        D_gamma.seekg((loc * NUM_BYTES_PER_SDB_ELEMENT), std::ios::beg);
        D_gamma.read(reinterpret_cast<char *>(net_buf), NUM_BYTES_PER_SDB_ELEMENT);
        mpz_import(tmp, NUM_BYTES_PER_SDB_ELEMENT, 1, 1, 1, 0, net_buf);
        mpz_class d_gamma = mpz_class(tmp);

        if (I > N)
        {
            memset(&read_entry, 0, sizeof(read_entry));
        }
        else
        {
            /* Match with the expection. Both the content as well as the index. */
            read_pdb_entry(pdb, (I - 1), read_entry);
        }

        mpz_import(tmp, sizeof(read_entry.element), 1, 1, 1, 0, read_entry.element);
        mpz_class d = mpz_class(tmp);

        /* FHE encrypt both the shares */
        Ciphertext<DCRTPoly> ct_d_alpha = FHE_Enc_SDBElement(d_alpha);
        Ciphertext<DCRTPoly> ct_d_gamma = FHE_Enc_SDBElement(d_gamma);

        /* Homorphically add those shares */
        Ciphertext<DCRTPoly> ct_d = ct_d_alpha + ct_d_gamma;

        /* Decrypt the resulting ciphertext */
        mpz_class dec_block_and_index, dec_block, dec_index;
        FHE_Dec_SDBElement(ct_d, dec_block_and_index);

        dec_index = (dec_block_and_index & ((1U << log_N) - 1U));
        dec_block = (dec_block_and_index >> log_N);

        if (I > N)
        {
            I = 0; // Because in that case the index will be just 0.
            // And in the case of dummy, the content match does not make any sense
        }

        if (dec_index != I)
        {
            PrintLog(LOG_LEVEL_TRACE, __FILE__, __LINE__, "Index does not match with the expected result. Expected: " + mpz_class(I).get_str(16) + " but got: " + mpz_class(dec_index).get_str(16));
            PrintLog(LOG_LEVEL_TRACE, __FILE__, __LINE__, "Randomly selected fetch index is: " + std::to_string(I) + ", and its location in shuffled database is: " + std::to_string(loc));
            PrintLog(LOG_LEVEL_TRACE, __FILE__, __LINE__, "d_alpha is:             " + d_alpha.get_str(16) + "\nd_gamma is:             " + d_gamma.get_str(16)+ "\ndec_block_and_index is: " + dec_block_and_index.get_str(16));

            test_failure++;
        } else {
            PrintLog(LOG_LEVEL_TRACE, __FILE__, __LINE__, "For I=" + std::to_string(I) + "\nd_alpha is:             " + d_alpha.get_str(16) + "\nd_gamma is:             " + d_gamma.get_str(16)+ "\ndec_block_and_index is: " + dec_block_and_index.get_str(16));
        }

        if ((dec_block != d) && (I != 0))
        {
            PrintLog(LOG_LEVEL_TRACE, __FILE__, __LINE__, "Data does not match with the expected result. Expected: " + d.get_str(16) + " but got: " + dec_block.get_str(16));
            PrintLog(LOG_LEVEL_TRACE, __FILE__, __LINE__, "Randomly selected fetch index is: " + std::to_string(I) + ", and its location in shuffled database is: " + std::to_string(loc));
            PrintLog(LOG_LEVEL_TRACE, __FILE__, __LINE__, "d_alpha is:             " + d_alpha.get_str(16) + "\nd_gamma is:             " + d_gamma.get_str(16)+ "\ndec_block_and_index is: " + dec_block_and_index.get_str(16));
        }
    }

    PrintLog(LOG_LEVEL_INFO, __FILE__, __LINE__, "Fetching test complete. Out of: " + std::to_string(NUM_FETCH_TEST_CNT) + " trails, fail count is: " + std::to_string(test_failure));


    pdb.close();
    D_alpha.close();
    D_gamma.close();

    return;
}
#endif

int main(int argc, char *argv[]){
    int ret = -1;

    /* Perform the basic initialization */
    InitSrv_beta();

    /* Process as per the command line arguments */
    if (argc >= 2) {
        if (std::string("gen_db").compare(std::string(argv[1]))==0) {
            ret = CreateRandomDatabase();
        } else if (std::string("one_time_init").compare(std::string(argv[1]))==0) {
            // Perform one-time initialization for server beta
            ret = OneTimeInit_beta();
        } else if (std::string("per_epoch_operations").compare(std::string(argv[1]))==0) {
            // Perform per-epoch initialization for server beta
            ret = PerEpochOperations_beta();
        } else if (std::string("clear_epoch_state").compare(std::string(argv[1]))==0) {
            // Clear the existing state of current epoch, start as if this is the first request of the epoch
            // TODO: Clear the existing state of current epoch, start as if this is the first request of the epoch
        } else if (std::string("process_request").compare(std::string(argv[1]))==0) {
            // Start from last saved state
            ret = ProcessClientRequest_beta();
        } else if (std::string("test").compare(std::string(argv[1]))==0) {
            TestSrv_beta();
        } else if (std::string("perf").compare(std::string(argv[1]))==0) {
            if (argc < 3){
                PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Performance measurement option requires at least three command line parameters. Usage: server_beta perf [srv_avg_online_time]");
            }else{
                if (std::string("srv_avg_online_time").compare(std::string(argv[2]))==0){
                    (void)Perf_avg_online_server_time_beta();
                }
            }
        } else {
            PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Unknown command line argument:"+ std::string(argv[1]));
            PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Improper command line arguments. Usage: server_beta <gen_db|one_time_init|per_epoch_operations|clear_epoch_state|process_request>");
        }
    } else {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Improper command line arguments. Usage: server_beta <gen_db|one_time_init|per_epoch_operations|clear_epoch_state|continue>");
    }

    FinSrv_beta();

    return 0;
}

static int TestShelterDPFSearch_beta() {
    Serial::DeserializeFromFile(ONE_TIME_MATERIALS_LOCATION_BETA + "sk_F.bin", sk_F, SerType::BINARY);
    // Sending the decryption key to the server alpha, so that it can verify the operation of the DPF
    (void)sendAll(sock_beta_alpha_con, Serial::SerializeToString(sk_F).c_str(), Serial::SerializeToString(sk_F).size());

    return 0;
}

static int TestClientProcessing_beta(){
    int ret = -1;
    size_t received_sz = 0;

    /* 1.b.1 Select random h_{\beta 0} */
    mpz_class h_beta0 = ElGamal_randomGroupElement();
    /* Also select its inverse */
    mpz_class h_beta0_1;
    mpz_invert(h_beta0_1.get_mpz_t(), h_beta0.get_mpz_t(), p.get_mpz_t());

    // 2.1.b Receive first component of E(T_I.h_{\alpha 1})
    ret = recvAll(sock_beta_alpha_con, net_buf, sizeof(net_buf), &received_sz);
    if (ret != 0) {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Server Beta: Failed to receive first component of E(T_I.h_{\\alpha 1} from Server Alpha");
        close(sock_beta_alpha_con);
        return ret;
    }
    mpz_class E_T_I_h_alpha1_1 = mpz_class(std::string(net_buf, received_sz));

    //2.2.b Receive second component of E(T_I.h_{\alpha 1})
    ret = recvAll(sock_beta_alpha_con, net_buf, sizeof(net_buf), &received_sz);
    if (ret != 0) {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Server Beta: Failed to receive second component of E(T_I.h_{\\alpha 1} from Server Alpha");
        close(sock_beta_alpha_con);
        return ret;
    }

    //3. Decrypt to T_I.h_{\\alpha 1}
    mpz_class E_T_I_h_alpha1_2 = mpz_class(std::string(net_buf, received_sz));
    std::pair<mpz_class, mpz_class> E_T_I_h_alpha1_local = std::make_pair(E_T_I_h_alpha1_1, E_T_I_h_alpha1_2);
    mpz_class T_I_h_alpha1 = ElGamal_decrypt(E_T_I_h_alpha1_local, sk_E);

    //PrintLog(LOG_LEVEL_INFO, __FILE__, __LINE__, "Decrypted T_I.h_{\\alpha 1}: " + T_I_h_alpha1.get_str());

    /* 4.b.1 Compute T_I.h_{\alpha 1}h_{\beta 0} */
    mpz_class T_I_h_alpha1_h_beta0 = (T_I_h_alpha1*h_beta0) % p;
    /* 4.b.2 Send T_I.h_{\alpha 1}h_{\beta 0} to server alpha */
    (void)sendAll(sock_beta_alpha_con, T_I_h_alpha1_h_beta0.get_str().c_str(), T_I_h_alpha1_h_beta0.get_str().size());

    /********************************************************************
     * For simulating communication with the client*********************/
    ret = recvAll(sock_beta_alpha_con, net_buf, sizeof(net_buf), &received_sz);
    if (ret != 0) {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Server Beta: Failed to receive FHE ciphertext of T_star.h_{\\alpha 2}.h_{\\beta 0} from Server Alpha");
        close(sock_beta_alpha_con);
        return ret;
    }

    mpz_class random_value = rng.get_z_bits(PLAINTEXT_PIR_BLOCK_DATA_SIZE+P_BITS);
    (void)sendAll(sock_beta_alpha_con, random_value.get_str().c_str(), random_value.get_str().size());

    /* End of simulation of client communication */
    return ret;
}

static void Perf_avg_online_server_time_beta() {

    /* Nothing is required for server beta */
    return;
}