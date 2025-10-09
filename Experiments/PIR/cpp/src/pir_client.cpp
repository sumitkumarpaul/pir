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

#define ONE_TIME_MATERIALS_LOCATION_CLIENT std::string("/mnt/sumit/PIR_CLIENT/ONE_TIME_MATERIALS/")
#define PER_EPOCH_MATERIALS_LOCATION_CLIENT std::string("/mnt/sumit/PIR_CLIENT/PER_EPOCH_MATERIALS/")

static int sock_client_to_alpha = -1, sock_client_to_beta = -1, sock_client_to_gamma = -1;
static char net_buf[NET_BUF_SZ] = {0};

// Function declarations
static int InitClient();
static int OneTimeInit_client();
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

    //Save parameters to local files
    export_to_file_from_mpz_class(ONE_TIME_MATERIALS_LOCATION_CLIENT + "p.bin", p);
    export_to_file_from_mpz_class(ONE_TIME_MATERIALS_LOCATION_CLIENT + "q.bin", q);
    export_to_file_from_mpz_class(ONE_TIME_MATERIALS_LOCATION_CLIENT + "g.bin", g);
    export_to_file_from_mpz_class(ONE_TIME_MATERIALS_LOCATION_CLIENT + "g_q.bin", g_q);
    export_to_file_from_mpz_class(ONE_TIME_MATERIALS_LOCATION_CLIENT + "r.bin", r);
    export_to_file_from_mpz_class(ONE_TIME_MATERIALS_LOCATION_CLIENT + "pk_E.bin", pk_E);
    export_to_file_from_mpz_class(ONE_TIME_MATERIALS_LOCATION_CLIENT + "pk_E_q.bin", pk_E_q);

    //TODO Segmentation fault Serial::SerializeToFile(ONE_TIME_MATERIALS_LOCATION_CLIENT + "FHEcryptoContext.bin", FHEcryptoContext, SerType::BINARY);
    Serial::SerializeToFile(ONE_TIME_MATERIALS_LOCATION_CLIENT + "pk_F.bin", pk_F, SerType::BINARY);

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

int main(int argc, char *argv[])
{
    int ret = -1;
    uint64_t I;

    /* Process as per the command line arguments */
    if (argc >= 2) {
        try
        {
            I = std::stoull(argv[1]);

            if ((I > N) || (I < 1))
            {
                PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "The database contains: " + std::to_string(N) + " blocks. Please enter a value in between 1 and " + std::to_string(N));
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