#include <assert.h>
#include <stdio.h>
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

#define MATERIALS_LOCATION_CLIENT std::string("/mnt/sumit/PIR_CLIENT/")
#define TMP_FILE std::string("/dev/shm/tmp_client")

static int sock_client_to_alpha = -1, sock_client_to_beta = -1, sock_client_to_gamma = -1;
static char net_buf[NET_BUF_SZ] = {0};

// Function declarations
static int InitClient();
static int OneTimeInit_client();
static int ShelterTagDetermination_Client(uint64_t I);
static int ObliDecReturn_Client(uint64_t* p_received_index);
static int FinClient();

static void TestClient();

// Function definitions
static int InitClient(){
    int ret = -1;
    // Initialize random number generation
    std::random_device rd;
    unsigned long seed = (static_cast<unsigned long>(rd()) << 1) ^ rd();
    rng.seed(seed); // seed() seeds the gmp_randclass    
    
    // Client opens three connecting sockets with three servers
    InitConnectingSocket(SERVER_ALPHA_IP, ALPHA_LISTENING_TO_CLIENT_PORT, &sock_client_to_alpha);
    PrintLog(LOG_LEVEL_TRACE, __FILE__, __LINE__, "Established connection with Server Alpha");

    InitConnectingSocket(SERVER_BETA_IP, BETA_LISTENING_TO_CLIENT_PORT, &sock_client_to_beta);
    PrintLog(LOG_LEVEL_TRACE, __FILE__, __LINE__, "Established connection with Server Beta");

    InitConnectingSocket(SERVER_GAMMA_IP, GAMMA_LISTENING_TO_CLIENT_PORT, &sock_client_to_gamma);
    PrintLog(LOG_LEVEL_TRACE, __FILE__, __LINE__, "Established connection with Server Gamma");

    // This actually corresponds to receiving published materials from the server beta
    ret = OneTimeInit_client();

    PrintLog(LOG_LEVEL_INFO, __FILE__, __LINE__, "Client initialization complete");

    return ret;
}

static int OneTimeInit_client() {
    size_t received_sz = 0;
    int ret_recv = 0;

    // Receive all the parameters from server beta
    // Receive p
    ret_recv = recvAll(sock_client_to_beta, net_buf, sizeof(net_buf), &received_sz);
    if (ret_recv != 0)
    {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Failed to receive p from Server Beta");
        return -1;
    }
    p = mpz_class(std::string(net_buf, received_sz));

    // Receive q
    ret_recv = recvAll(sock_client_to_beta, net_buf, sizeof(net_buf), &received_sz);
    if (ret_recv != 0)
    {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Failed to receive q from Server Beta");
        return -1;
    }
    q = mpz_class(std::string(net_buf, received_sz));

    // Receive g
    ret_recv = recvAll(sock_client_to_beta, net_buf, sizeof(net_buf), &received_sz);
    if (ret_recv != 0)
    {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Failed to receive g from Server Beta");
        return -1;
    }
    g = mpz_class(std::string(net_buf, received_sz));

    // Receive g_dashed
    ret_recv = recvAll(sock_client_to_beta, net_buf, sizeof(net_buf), &received_sz);
    if (ret_recv != 0)
    {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Failed to receive g_dashed from Server Beta");
        return -1;
    }
    g_dashed = mpz_class(std::string(net_buf, received_sz));

    // Receive r
    ret_recv = recvAll(sock_client_to_beta, net_buf, sizeof(net_buf), &received_sz);
    if (ret_recv != 0)
    {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Failed to receive r from Server Beta");
        return -1;
    }
    r = mpz_class(std::string(net_buf, received_sz));

    // Receive pk_E
    ret_recv = recvAll(sock_client_to_beta, net_buf, sizeof(net_buf), &received_sz);
    if (ret_recv != 0)
    {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Failed to receive pk_E from Server Beta");
        return -1;
    }
    pk_E = mpz_class(std::string(net_buf, received_sz));

    // Receive pk_E_dashed
    ret_recv = recvAll(sock_client_to_beta, net_buf, sizeof(net_buf), &received_sz);
    if (ret_recv != 0)
    {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Failed to receive pk_E_dashed from Server Beta");
        return -1;
    }
    pk_E_dashed = mpz_class(std::string(net_buf, received_sz));

    // Receive FHEcryptoContext
    ret_recv = recvFile(sock_client_to_beta, net_buf, sizeof(net_buf), TMP_FILE);
    if (ret_recv != 0)
    {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Failed to receive FHEcryptoContext from Server Beta");
        return -1;
    }

    Serial::DeserializeFromFile(TMP_FILE, FHEcryptoContext, SerType::BINARY);

    // Receive pk_F
    ret_recv = recvFile(sock_client_to_beta, net_buf, sizeof(net_buf), TMP_FILE);
    if (ret_recv != 0)
    {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Failed to receive pk_F from Server Beta");
        return -1;
    }
    Serial::DeserializeFromFile(TMP_FILE, pk_F, SerType::BINARY);

    // Receive E_dashed_Rho.first
    ret_recv = recvAll(sock_client_to_beta, net_buf, sizeof(net_buf), &received_sz);
    if (ret_recv != 0)
    {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Failed to receive E_dashed_Rho.first from Server Beta");
        return -1;
    }
    E_dashed_Rho.first = mpz_class(std::string(net_buf, received_sz));

    // Receive E_dashed_Rho.second
    ret_recv = recvAll(sock_client_to_beta, net_buf, sizeof(net_buf), &received_sz);
    if (ret_recv != 0)
    {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Failed to receive E_dashed_Rho.second from Server Beta");
        return -1;
    }
    E_dashed_Rho.second = mpz_class(std::string(net_buf, received_sz));

    //Save parameters to local files
    export_to_file_from_mpz_class(MATERIALS_LOCATION_CLIENT + "p.bin", p);
    export_to_file_from_mpz_class(MATERIALS_LOCATION_CLIENT + "q.bin", q);
    export_to_file_from_mpz_class(MATERIALS_LOCATION_CLIENT + "g.bin", g);
    export_to_file_from_mpz_class(MATERIALS_LOCATION_CLIENT + "g_dashed.bin", g_dashed);
    export_to_file_from_mpz_class(MATERIALS_LOCATION_CLIENT + "r.bin", r);
    export_to_file_from_mpz_class(MATERIALS_LOCATION_CLIENT + "pk_E.bin", pk_E);
    export_to_file_from_mpz_class(MATERIALS_LOCATION_CLIENT + "pk_E_dashed.bin", pk_E_dashed);
    export_to_file_from_mpz_class(MATERIALS_LOCATION_CLIENT + "E_dashed_Rho_1.bin", E_dashed_Rho.first);
    export_to_file_from_mpz_class(MATERIALS_LOCATION_CLIENT + "E_dashed_Rho_2.bin", E_dashed_Rho.second);

    Serial::SerializeToFile(MATERIALS_LOCATION_CLIENT + "FHEcryptoContext.bin", FHEcryptoContext, SerType::BINARY);
    Serial::SerializeToFile(MATERIALS_LOCATION_CLIENT + "pk_F.bin", pk_F, SerType::BINARY);

    return 0;
}

static int FinClient(){
    int ret = -1;

    // Close the sockets
    if (sock_client_to_alpha != -1) {
        close(sock_client_to_alpha);
        sock_client_to_alpha = -1;
    }
    if (sock_client_to_beta != -1) {
        close(sock_client_to_beta);
        sock_client_to_beta = -1;
    }
    if (sock_client_to_gamma != -1) {
        close(sock_client_to_gamma);
        sock_client_to_gamma = -1;
    }    

    PrintLog(LOG_LEVEL_INFO, __FILE__, __LINE__, "Finalized the client");

    return ret;
}


static int ShelterTagDetermination_Client(uint64_t I){
    int ret = -1;
    mpz_class mpz_I= mpz_class(I);
    mpz_class h_C0;
    mpz_class h_C0_inv;
    mpz_class h_C1;
    mpz_class h_C1_inv;
    const mpz_class q_minus_one = q - 1;
    mpz_class Rho_pow_h_beta0;
    mpz_class h_C0_pow_h_C1;
    mpz_class h_C0_pow_h_beta0_h_C1;
    mpz_class h_C0_pow_h_beta0;
    mpz_class Rho_pow_I_h_beta0;
    mpz_class Rho_pow_I_h_C0__pow_h_beta0;
    std::pair<mpz_class, mpz_class> E_g_pow_Rho_pow_I__h_C0;
    std::pair<mpz_class, mpz_class> E_g_pow_Rho_pow_I;
    size_t received_sz = 0;
    int ret_recv = 0;

    PrintLog(LOG_LEVEL_SPECIAL, __FILE__, __LINE__, "Request fetching start");

    /* Step 1.a: Choose h_C0 and the inverse of it */
    h_C0 = rng.get_z_range(q-1)+1; //i.e., within ZZ*_q
    mpz_invert(h_C0_inv.get_mpz_t(), h_C0.get_mpz_t(), q.get_mpz_t());//Compute the inverse of h_C0 in mod q

    /* Step 2.a: Choose h_C1 and the inverse of it */
    /* For some reason, inverse is not working. So, keep trying till the inverse work */
    do {
        h_C1 = rng.get_z_range(q-2)+1; //i.e., within ZZ*_{q-1}
        mpz_invert(h_C1_inv.get_mpz_t(), h_C1.get_mpz_t(), q_minus_one.get_mpz_t());//Compute the inverse of h_C1
    } while (((h_C1*h_C1_inv)%q_minus_one) != 1);

    PrintLog(LOG_LEVEL_INFO, __FILE__, __LINE__, "HERE");

    /* Step 3.2.2: Receive Rho_pow_h_beta0 from Server Beta */
    ret_recv = recvAll(sock_client_to_beta, net_buf, sizeof(net_buf), &received_sz);
    if (ret_recv != 0)
    {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Failed to receive Rho_pow_h_beta0 from the Server Beta");
        return -1;
    }
    Rho_pow_h_beta0 = mpz_class(std::string(net_buf, received_sz));
    PrintLog(LOG_LEVEL_INFO, __FILE__, __LINE__, "HERE");

    /* Step 4: Compute and send h_C0_pow_h_C1 to the Server Beta */
    /* Step 4.1: Compute h_C0_pow_h_C1 */
    mpz_powm(h_C0_pow_h_C1.get_mpz_t(), h_C0.get_mpz_t(), h_C1.get_mpz_t(), q.get_mpz_t());
    /* Step 4.2.1: Then send to the Server Beta */
    (void)sendAll(sock_client_to_beta, h_C0_pow_h_C1.get_str().c_str(), h_C0_pow_h_C1.get_str().size());

    /* Step 5.2.2: Receive h_C0_pow_h_beta0_h_C1 from Server Beta */
    ret_recv = recvAll(sock_client_to_beta, net_buf, sizeof(net_buf), &received_sz);
    if (ret_recv != 0)
    {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Failed to receive h_C0_pow_h_beta0_h_C1 from the Server Beta");
        return -1;
    }
    h_C0_pow_h_beta0_h_C1 = mpz_class(std::string(net_buf, received_sz));

    PrintLog(LOG_LEVEL_INFO, __FILE__, __LINE__, "HERE");
    /* Step 6: Compute h_C0_pow_h_beta0 */
    mpz_powm(h_C0_pow_h_beta0.get_mpz_t(), h_C0_pow_h_beta0_h_C1.get_mpz_t(), h_C1_inv.get_mpz_t(), q.get_mpz_t());//Apply the inverse

    /* Step 7: Compute and send Rho_pow_I_h_C0__pow_h_beta0 to Server Beta */
    /* Step 7.1: First compute Rho_pow_I_h_C0__pow_h_beta0 */
    mpz_powm(Rho_pow_I_h_beta0.get_mpz_t(), Rho_pow_h_beta0.get_mpz_t(), mpz_I.get_mpz_t(), q.get_mpz_t());
    Rho_pow_I_h_C0__pow_h_beta0 = (Rho_pow_I_h_beta0 * h_C0_pow_h_beta0) % q;
    /* Step 7.2.1: Then send to the Server Beta */
    (void)sendAll(sock_client_to_beta, Rho_pow_I_h_C0__pow_h_beta0.get_str().c_str(), Rho_pow_I_h_C0__pow_h_beta0.get_str().size());

    /* Step 9.3.2 Receive the first component of E_g_pow_Rho_pow_I__h_C0 */
    ret_recv = recvAll(sock_client_to_beta, net_buf, sizeof(net_buf), &received_sz);
    if (ret_recv != 0)
    {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Failed to receive E_g_pow_Rho_pow_I__h_C0.first from the Server Beta");
        return -1;
    }
    PrintLog(LOG_LEVEL_INFO, __FILE__, __LINE__, "HERE");
    E_g_pow_Rho_pow_I__h_C0.first = mpz_class(std::string(net_buf, received_sz));

    /* Step 9.4.2 Receive the second component of E_g_pow_Rho_pow_I__h_C0 */
    ret_recv = recvAll(sock_client_to_beta, net_buf, sizeof(net_buf), &received_sz);
    if (ret_recv != 0)
    {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Failed to receive E_g_pow_Rho_pow_I__h_C0.second from the Server Beta");
        return -1;
    }
    PrintLog(LOG_LEVEL_INFO, __FILE__, __LINE__, "HERE");
    E_g_pow_Rho_pow_I__h_C0.second = mpz_class(std::string(net_buf, received_sz));

    /* Step 10: Homomorphically remove h_C0 and send E_g_pow_Rho_pow_I to Server Alpha */
    /* Step 10.1: First compute E_g_pow_Rho_pow_I */
    E_g_pow_Rho_pow_I = ElGamal_exp_ct(E_g_pow_Rho_pow_I__h_C0, h_C0_inv, pk_E);
    /* Step 10.2.1: Then send the first component of the ciphertext to the Server Alpha */
    (void)sendAll(sock_client_to_alpha, E_g_pow_Rho_pow_I.first.get_str().c_str(), E_g_pow_Rho_pow_I.first.get_str().size());
    /* Step 10.3.1: Then send the first component of the ciphertext to the Server Alpha */
    (void)sendAll(sock_client_to_alpha, E_g_pow_Rho_pow_I.second.get_str().c_str(), E_g_pow_Rho_pow_I.second.get_str().size());



    return ret;
}

static int ObliDecReturn_Client(uint64_t* p_received_index) {
    int ret = -1;
    Ciphertext<DCRTPoly> m_C_ct, refreshed_ct;
    mpz_class m_C;
    size_t received_sz = 0;
    int ret_recv = 0;
    mpz_class received_element, extracted_element, extracted_element_content, extracted_element_index;
    mpz_class extracted_part, received_part, m_C_part, mask;
    *p_received_index = 0;


    /* Step 1.a: Generate random mask */
    m_C = rng.get_z_bits((NUM_BYTES_PER_SDB_ELEMENT*8));

    /* Step 2.1: Generate ciphertext of the random mask */
    m_C_ct = FHE_bitwise_Enc_SDBElement(m_C);

    /* Step 2.2: Send corresponding ciphertext to server gamma */
    Serial::SerializeToFile(TMP_FILE, m_C_ct, SerType::BINARY);
    (void)sendFile(sock_client_to_gamma, net_buf, sizeof(net_buf), TMP_FILE);    

    /* 9.2 Receive decryption result */
    ret_recv = recvAll(sock_client_to_beta, net_buf, sizeof(net_buf), &received_sz);
    if (ret_recv != 0)
    {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Failed to receive the decryption result from the Server Beta");
        return -1;
    }
    PrintLog(LOG_LEVEL_SPECIAL, __FILE__, __LINE__, "Completed receiving requested item");

    received_element = mpz_class(std::string(net_buf, received_sz));

    /* 10. Remove mask */
    /* Similar but reverse logic of per-epoch operations for server beta */
    mpz_xor(extracted_element.get_mpz_t(), received_element.get_mpz_t(), m_C.get_mpz_t());

    /* Extract result */
    extracted_element_content = (extracted_element >> log_N);
    extracted_element_index = (extracted_element & ((1UL << log_N) - 1UL));
    PrintLog(LOG_LEVEL_TRACE, __FILE__, __LINE__, "Received element is (HEX): " + extracted_element.get_str(16));
    PrintLog(LOG_LEVEL_TRACE, __FILE__, __LINE__, "Extracted block content is (HEX): " + extracted_element_content.get_str(16));
    PrintLog(LOG_LEVEL_TRACE, __FILE__, __LINE__, "Received index is (DEC): " + extracted_element_index.get_str());

    *p_received_index = extracted_element_index.get_ui();

exit:
    return ret;
}

int main(int argc, char *argv[])
{
    int ret = -1;
    uint64_t I, received_index;

    /* Process as per the command line arguments */
    if (argc >= 2) {
        try
        {
            I = std::stoull(argv[1]);

            if ((I > N) || (I < 1))
            {
                PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "The database contains: " + std::to_string(N) + " blocks. Please enter a value in between 1 and " + std::to_string(N));
                return ret;
            }
        }
        catch (const std::out_of_range &oor)
        {
            PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Value out of range..!! The database contains: " + std::to_string(N) + " blocks. Please enter a value in between 1 and " + std::to_string(N));
        }
        catch (const std::invalid_argument &ia)
        {
            PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Invalid argument (not a valid number). Usage: pir_client <Requested index from the database>");
        }

        /* Perform the basic initialization */
        InitClient();

        /* Shelter-tag determination */
        ShelterTagDetermination_Client(I);

        /* Oblivious decryption and return */
        ObliDecReturn_Client(&received_index);

        if (received_index == I){
            PrintLog(LOG_LEVEL_INFO, __FILE__, __LINE__, "Success..index of the received block matches with requested index: " + std::to_string(I));
            ret = 0;
        } else {
            PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Received index is: " + std::to_string(received_index) + " which does not matches with requested index: " + std::to_string(I));
        }
    }
    else
    {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Improper command line arguments. Usage: pir_client <Requested index from the database>");
    }

    if (ret == 0) {
        TestClient();
    }

    FinClient();

    return ret;
}

void TestClient(){
    return;
}