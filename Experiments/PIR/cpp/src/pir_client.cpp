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

    // Receive g_q
    ret_recv = recvAll(sock_client_to_beta, net_buf, sizeof(net_buf), &received_sz);
    if (ret_recv != 0)
    {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Failed to receive g_q from Server Beta");
        return -1;
    }
    g_q = mpz_class(std::string(net_buf, received_sz));

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

    // Receive pk_E_q
    ret_recv = recvAll(sock_client_to_beta, net_buf, sizeof(net_buf), &received_sz);
    if (ret_recv != 0)
    {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Failed to receive pk_E_q from Server Beta");
        return -1;
    }
    pk_E_q = mpz_class(std::string(net_buf, received_sz));

    // Receive FHEcryptoContext
    ret_recv = recvAll(sock_client_to_beta, net_buf, sizeof(net_buf), &received_sz);
    if (ret_recv != 0)
    {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Failed to receive FHEcryptoContext from Server Beta");
        return -1;
    }

    Serial::DeserializeFromString(FHEcryptoContext, std::string(net_buf, received_sz));

    // Receive pk_F
    ret_recv = recvAll(sock_client_to_beta, net_buf, sizeof(net_buf), &received_sz);
    if (ret_recv != 0)
    {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Failed to receive pk_F from Server Beta");
        return -1;
    }
    Serial::DeserializeFromString(pk_F, std::string(net_buf, received_sz));

    // Receive E_q_Rho.first
    ret_recv = recvAll(sock_client_to_beta, net_buf, sizeof(net_buf), &received_sz);
    if (ret_recv != 0)
    {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Failed to receive E_q_Rho.first from Server Beta");
        return -1;
    }
    E_q_Rho.first = mpz_class(std::string(net_buf, received_sz));

    // Receive E_q_Rho.second
    ret_recv = recvAll(sock_client_to_beta, net_buf, sizeof(net_buf), &received_sz);
    if (ret_recv != 0)
    {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Failed to receive E_q_Rho.second from Server Beta");
        return -1;
    }
    E_q_Rho.second = mpz_class(std::string(net_buf, received_sz));

    //Save parameters to local files
    export_to_file_from_mpz_class(MATERIALS_LOCATION_CLIENT + "p.bin", p);
    export_to_file_from_mpz_class(MATERIALS_LOCATION_CLIENT + "q.bin", q);
    export_to_file_from_mpz_class(MATERIALS_LOCATION_CLIENT + "g.bin", g);
    export_to_file_from_mpz_class(MATERIALS_LOCATION_CLIENT + "g_q.bin", g_q);
    export_to_file_from_mpz_class(MATERIALS_LOCATION_CLIENT + "r.bin", r);
    export_to_file_from_mpz_class(MATERIALS_LOCATION_CLIENT + "pk_E.bin", pk_E);
    export_to_file_from_mpz_class(MATERIALS_LOCATION_CLIENT + "pk_E_q.bin", pk_E_q);
    export_to_file_from_mpz_class(MATERIALS_LOCATION_CLIENT + "E_q_Rho_1.bin", E_q_Rho.first);
    export_to_file_from_mpz_class(MATERIALS_LOCATION_CLIENT + "E_q_Rho_2.bin", E_q_Rho.second);

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

    PrintLog(LOG_LEVEL_SPECIAL, __FILE__, __LINE__, "Finalized the client");

    return ret;
}


static int ShelterTagDetermination_Client(uint64_t I){
    int ret = -1;
    std::pair<mpz_class, mpz_class> E_g_pow_Rho_pow_I__mul__h_C_h_alpha0;
    std::pair<mpz_class, mpz_class> E_q_Rho_pow_I;
    mpz_class h_C, h_C_1;
    std::pair<mpz_class, mpz_class> E_q_h_C;
    std::pair<mpz_class, mpz_class> E_q_Rho_pow_I__mul__h_C;
    std::pair<mpz_class, mpz_class> E_g_pow_Rho_pow_I__mul_h_alpha0;

    size_t received_sz = 0;
    int ret_recv = 0;
    
    /* Step 1 */
    E_q_Rho_pow_I = ElGamal_q_exp_ct(E_q_Rho, mpz_class(I), pk_E_q);

    /* Step 2.1 */
    h_C = rng.get_z_range(q-1)+1;//i.e., within ZZ_q*
    /* Step 2.1.1 figure out its inverse */
    mpz_invert(h_C_1.get_mpz_t(), h_C.get_mpz_t(), q.get_mpz_t());

    /* Step 2.2.1 */
    E_q_h_C = ElGamal_q_encrypt(h_C, pk_E_q);

    /* Step 2.2.2 */
    E_q_Rho_pow_I__mul__h_C = ElGamal_q_mult_ct(E_q_Rho_pow_I, E_q_h_C);

    /* Step 2.3.1 */
    (void)sendAll(sock_client_to_beta, E_q_Rho_pow_I__mul__h_C.first.get_str().c_str(), E_q_Rho_pow_I__mul__h_C.first.get_str().size());
    (void)sendAll(sock_client_to_beta, E_q_Rho_pow_I__mul__h_C.second.get_str().c_str(), E_q_Rho_pow_I__mul__h_C.second.get_str().size());

    // Step 6.1 Receive the first component of E_g_pow_Rho_pow_I__mul__h_C_h_alpha0
    ret_recv = recvAll(sock_client_to_alpha, net_buf, sizeof(net_buf), &received_sz);
    if (ret_recv != 0)
    {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Failed to receive E_g_pow_Rho_pow_I__mul__h_C_h_alpha0.first from the Server Alpha");
        return -1;
    }
    E_g_pow_Rho_pow_I__mul__h_C_h_alpha0.first = mpz_class(std::string(net_buf, received_sz));

    // Step 6.2 Receive the second component of E_g_pow_Rho_pow_I__mul__h_C_h_alpha0
    ret_recv = recvAll(sock_client_to_alpha, net_buf, sizeof(net_buf), &received_sz);
    if (ret_recv != 0)
    {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Failed to receive E_g_pow_Rho_pow_I__mul__h_C_h_alpha0.second from the Server Alpha");
        return -1;
    }
    E_g_pow_Rho_pow_I__mul__h_C_h_alpha0.second = mpz_class(std::string(net_buf, received_sz));

    // Step 7.1 Semi-homomorphically raises that to h_C_1 under ElGamal encryption in GG to remove h_C
    E_g_pow_Rho_pow_I__mul_h_alpha0 = ElGamal_exp_ct(E_g_pow_Rho_pow_I__mul__h_C_h_alpha0, h_C_1, pk_E);

    /* Step 7.2 Send Both the components of the resulting ciphtext to server Alpha */
    (void)sendAll(sock_client_to_alpha, E_g_pow_Rho_pow_I__mul_h_alpha0.first.get_str().c_str(), E_g_pow_Rho_pow_I__mul_h_alpha0.first.get_str().size());
    (void)sendAll(sock_client_to_alpha, E_g_pow_Rho_pow_I__mul_h_alpha0.second.get_str().c_str(), E_g_pow_Rho_pow_I__mul_h_alpha0.second.get_str().size());



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

    /* Step 1: Generate random mask */
    m_C = rng.get_z_bits((PLAINTEXT_PIR_BLOCK_DATA_SIZE +  log_N));

    /* Step 2.1: Generate ciphertext of the random mask */
    m_C_ct = FHE_Enc_SDBElement(m_C);

    /* Step 2.2: Send corresponding ciphertext to server gamma */
    (void)sendAll(sock_client_to_gamma, Serial::SerializeToString(m_C_ct).c_str(), Serial::SerializeToString(m_C_ct).size());

    /* 6.2 Receive decryption result */
    ret_recv = recvAll(sock_client_to_beta, net_buf, sizeof(net_buf), &received_sz);
    if (ret_recv != 0)
    {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Failed to receive E_g_pow_Rho_pow_I__mul__h_C_h_alpha0.first from the Server Alpha");
        return -1;
    }
    received_element = mpz_class(std::string(net_buf, received_sz));

    /* 7. Remove mask */
    /* Similar but reverse logic of per-epoch operations for server beta */
    extracted_element = mpz_class(0);
    mask = mpz_class((1 << PLAINTEXT_FHE_BLOCK_SIZE) - 1);

    for (unsigned int i = 0; i < TOTAL_NUM_FHE_BLOCKS_PER_ELEMENT; i++)
    {
        /* Extract least significant PLAINTEXT_FHE_BLOCK_SIZE-bits of d and d_alpha */
        m_C_part = (m_C & mask);
        received_part = (received_element & mask);

        /* Compute the difference between two parts. And take only PLAINTEXT_FHE_BLOCK_SIZE-bits */
        extracted_part = (received_part + m_C_part) & mask;

        /* Append the part at the proper location */
        extracted_element = (extracted_element | extracted_part);

        mask = mask << PLAINTEXT_FHE_BLOCK_SIZE;
    }

    /* Extract result */
    extracted_element_content = (extracted_element >> log_N);
    extracted_element_index = (extracted_element & ((1U << log_N) - 1U));
    PrintLog(LOG_LEVEL_TRACE, __FILE__, __LINE__, "Received element is (HEX): " + extracted_element.get_str(16));
    PrintLog(LOG_LEVEL_TRACE, __FILE__, __LINE__, "Extracted block content is (HEX): " + extracted_element_content.get_str(16));
    PrintLog(LOG_LEVEL_TRACE, __FILE__, __LINE__, "Received index is (DEC): " + extracted_element_index.get_str());

    *p_received_index = extracted_element_index.get_ui();

    /* Additional steps for refreshing ciphertext */
    refreshed_ct = FHE_Enc_SDBElement(extracted_element);
    /* Send refreshed ciphertext to Server Gamma */
    (void)sendAll(sock_client_to_gamma, Serial::SerializeToString(refreshed_ct).c_str(), Serial::SerializeToString(refreshed_ct).size());    

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
                return -1;
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
            PrintLog(LOG_LEVEL_SPECIAL, __FILE__, __LINE__, "Index of the received block matches with requested index: " + std::to_string(I));
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

    return 0;
}

void TestClient(){
    return;
}