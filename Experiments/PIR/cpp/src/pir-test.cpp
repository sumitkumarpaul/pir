#include <chrono>
#include <assert.h>
#include <stdio.h>
#include "fss/fss-common.h"
#include "fss/fss-server.h"
#include "fss/fss-client.h"
#include <immintrin.h>  // Include AVX header
#include <omp.h>
#include <vector>
#include "pir_common.h"
#include <random>

#define PROFILE
#define NUM_CPU_CORES 16

#include <fstream>
#include <iostream>
#include <iterator>
#include <thread>
#include <cstring>
#include <sys/socket.h>
#include <netinet/in.h>
#include <arpa/inet.h>
#include <unistd.h>


// ...existing code...

#define N 50000 // Size of the database, can be adjusted as needed
// The value of N will determine the bitlength during the client initialization
#define NUM_TAG_BITS 3072 // 16 bits can represent up to 65536, which is more than enough for N=50000
#define B 512 // Block size in bits, can be adjusted as needed
// And number of bits determine the evalution time drastically
static mpz_class DB[N][2]; // Database to store values, each entry is a pair, {Tag, Block-content}.

// Client asks for the element located at location i, where 0 <= i < N 
int PIR_Experiment(int I) {
    //Actually I is only used to verify the correctness of the PIR scheme
    //The tester knows that T_sh is the tag of the element at location I

    // Set up variables
    Fss fClient, fServer;
    ServerKeyEq k0;
    ServerKeyEq k1;

    //TODO: Our scheme has the benefit that the entire FSS lookup can be achieved without any disk access
    //Since the entire stash may reside in the RAM

    // Initialize the database
    gmp_randclass r(gmp_randinit_default);
    r.seed(time(NULL));
    
    for(size_t i=0; i<N; i++) {
        // Create random tag and block content
        DB[i][0] = r.get_z_bits(NUM_TAG_BITS); // Create a random tag
        DB[i][1] = r.get_z_bits(B); // Create a random block of B bits
    }  
    
    //Check the possible time requirement for updating the tags of the stash
    mpz_class p = r.get_z_bits(NUM_TAG_BITS); // Create a random tag
    mpz_nextprime(p.get_mpz_t(), p.get_mpz_t());
    mpz_class del_abc = r.get_z_bits(NUM_TAG_BITS); // Create a random multiplier for the tags

    auto t_exp_start = std::chrono::high_resolution_clock::now();
    
    #if 1// Modular multiplication test

    for (size_t i = 0; i < N; i += NUM_CPU_CORES)//TODO: Check the implementation of parallel execution
    {
#pragma omp parallel for
        for (int j = 0; j < NUM_CPU_CORES; ++j)
        {
            if ((i + j) < N)
            {
                mpz_class a = DB[i+j][0]; // Get the tag of the element at location i
                mpz_mul(a.get_mpz_t(), a.get_mpz_t(), del_abc.get_mpz_t()); // Perform multiplication operation
                mpz_mod(DB[i][0].get_mpz_t(), a.get_mpz_t(), p.get_mpz_t()); // Perform modular operation
            }
        }
    }

    #if 0//Non parallel version
    for(size_t i=0; i<N; i++) {
        mpz_class a = DB[i][0]; // Get the tag of the element at location i
        mpz_mul(a.get_mpz_t(), a.get_mpz_t(), del_abc.get_mpz_t()); // Perform multiplication operation
        mpz_mod(DB[i][0].get_mpz_t(), a.get_mpz_t(), p.get_mpz_t()); // Perform modular operation
    }
    #endif
    #endif
    auto t_exp_end = std::chrono::high_resolution_clock::now();
    std::cout << "Time for shelter tag update operation: " <<
     std::chrono::duration<double, std::milli>(t_exp_end - t_exp_start).count()
     << " ms" << endl;

    //Cheat and set the search tag
    mpz_class T_sh = DB[I][0]; // The tag of the element at location I 
    cout << "Search tag is: " << T_sh << endl;

    auto t_begin = std::chrono::high_resolution_clock::now();

    // Initialize client, use 32 bits in domain as example
    initializeClient(&fClient, NUM_TAG_BITS, 2); // If bit length is not set properly, then incorrect answer will be returned
    // bit length must be able to store ans0 and ans1
    
    auto t_initClient = std::chrono::high_resolution_clock::now();

    // Equality FSS test
    generateTreeEq(&fClient, &k0, &k1, T_sh, 1);//So that the point function will evaluate as 1 at location i, and zero elsewhere

    auto t_keyGen = std::chrono::high_resolution_clock::now();    
    
    // Initialize server
    initializeServer(&fServer, &fClient);

    auto t_initServer = std::chrono::high_resolution_clock::now();

    mpz_class ans0, ans1, fin;
    ans0 = 0;
    ans1 = 0;

    std::vector<mpz_class> thread_sums(NUM_CPU_CORES);
    for (size_t i = 0; i < N; i += NUM_CPU_CORES)
    {
        for (int t = 0; t < NUM_CPU_CORES; ++t)
            thread_sums[t] = 0;

#pragma omp parallel for
        for (int j = 0; j < NUM_CPU_CORES; ++j)
        {
            if ((i + j) < N)
            {
                auto y = evaluateEq(&fServer, &k0, DB[i + j][0]);// Evaluate the FSS for the tag
                thread_sums[j] += y * DB[i + j][1]; // Multiply the result with the block content
            }
        }
        for (int t = 0; t < NUM_CPU_CORES; ++t)
            ans0 += thread_sums[t];
    }

    auto t_evalServer0 = std::chrono::high_resolution_clock::now();

    for (size_t i = 0; i < N; i += NUM_CPU_CORES)
    {
        for (int t = 0; t < NUM_CPU_CORES; ++t)
            thread_sums[t] = 0;

#pragma omp parallel for
        for (int j = 0; j < NUM_CPU_CORES; ++j)
        {
            if ((i + j) < N)
            {
                auto y = evaluateEq(&fServer, &k1, DB[i + j][0]);// Evaluate the FSS for the tag
                thread_sums[j] += y * DB[i + j][1]; // Multiply the result with the block content
            }
        }
        for (int t = 0; t < NUM_CPU_CORES; ++t)
            ans1 += thread_sums[t];
    }

    fin = ans0 - ans1;
    cout << "PIR Test: value of ans0: " << ans0 << ", value of ans1: " << ans1 << endl;
    cout << "PIR Test: (it should be the value stored at DB[I]): " << fin << endl;
    cout << "The value stored at DB[I] is: " << DB[I][1] << endl;

    auto t_evalServer1 = std::chrono::high_resolution_clock::now();    
    
    if (fin != mpz_class(DB[I][1])) {
        cout << "!!!!!!!!!!!!!!!!!!! PIR Test: Error, the value returned does not match the expected value !!!!!!!!!!!!!" << endl;
    } else {
        cout << "PIR Test: Success, the value returned matches the expected value!" << endl;
    }    

    auto t_end = std::chrono::high_resolution_clock::now();

    std::cout << "Time to initialize client: " <<
     std::chrono::duration<double, std::milli>(t_initClient - t_begin).count()
     << " ms" << endl;

    std::cout << "Time to generate keys: " <<
     std::chrono::duration<double, std::milli>(t_keyGen - t_initClient).count()
     << " ms" << endl;

    std::cout << "Time to initialize server: " <<
     std::chrono::duration<double, std::milli>(t_initServer - t_keyGen).count()
     << " ms" << endl;

    std::cout << "Time to evaluate by server 0 " <<
     std::chrono::duration<double, std::milli>(t_evalServer0 - t_initServer).count()
     << " ms" << endl;

    std::cout << "Time to evaluate by server 1 " <<
     std::chrono::duration<double, std::milli>(t_evalServer1 - t_evalServer0).count()
     << " ms" << endl;

    std::cout << "Total time: " <<
     std::chrono::duration<double, std::milli>(t_end - t_begin).count()
     << " ms" << endl;
    return 1;
}




// Declaration for external use
void ElGamalEncryptorSocket(int port);


// Decryptor: acts as client, sends public key, receives ciphertext and message
void ElGamalDecryptorSocket(const std::string& server_ip, int port) {
    using namespace std::chrono;
    std::cout << "=== ElGamal Test (Socket) ===" << std::endl;
    ElGamal elgamal;
    int p_bits = 3072, q_bits = 256;

    auto t_start = high_resolution_clock::now();
    elgamal.setup(p_bits, q_bits);
    auto t_setup = high_resolution_clock::now();
    std::cout << "Setup complete." << std::endl;
    std::cout << "p (prime, " << p_bits << " bits): " << elgamal.p.get_str() << std::endl;
    std::cout << "q (subgroup order, " << q_bits << " bits): " << elgamal.q.get_str() << std::endl;
    std::cout << "g (generator): " << elgamal.g.get_str() << std::endl;
    auto t_params = high_resolution_clock::now();

    // Key generation
    auto [priv, pub] = elgamal.keyGen();
    auto t_keygen = high_resolution_clock::now();
    std::cout << "Private key x: " << priv.get_str() << std::endl;
    std::cout << "Public key y: " << pub.get_str() << std::endl;

    // Serialize parameters
    std::string p_str = elgamal.p.get_str();
    std::string q_str = elgamal.q.get_str();
    std::string g_str = elgamal.g.get_str();
    std::string pub_str = pub.get_str();

    int sock;
    char buffer[65536] = {0};
    InitConnectingSocket(server_ip, port, &sock);
    // Send public key and params
    std::string msg = p_str + "\n" + q_str + "\n" + g_str + "\n" + pub_str + "\n";
    send(sock, msg.c_str(), msg.size(), 0);
    // Receive ciphertext and message
    int valread = read(sock, buffer, sizeof(buffer));
    std::string received(buffer, valread);
    // Expect: c1\nc2\nmessage\n
    size_t pos = 0;
    std::vector<std::string> params;
    while ((pos = received.find("\n")) != std::string::npos) {
        params.push_back(received.substr(0, pos));
        received.erase(0, pos + 1);
    }
    if (params.size() < 3) {
        std::cerr << "Client: Invalid ciphertext message" << std::endl;
        close(sock); return;
    }
    std::string ciphertext1_str = params[0], ciphertext2_str = params[1], message_str = params[2];
    close(sock);

    // Deserialize ciphertext and message
    mpz_class c1(ciphertext1_str), c2(ciphertext2_str), message(message_str);

    // Decrypt
    mpz_class decrypted = elgamal.decrypt({c1, c2}, priv);
    std::cout << "Decrypted message: " << decrypted.get_str() << std::endl;
    std::cout << "Original message: " << message.get_str() << std::endl;

    // Verify
    if (decrypted == message)
        std::cout << "ElGamal Test: Success! Decrypted message matches original message." << std::endl;
    else
        std::cout << "ElGamal Test: Failure! Decrypted message does not match." << std::endl;

    auto t_end = high_resolution_clock::now();
    std::cout << "Total time: " << duration<double, std::milli>(t_end - t_start).count() << " ms" << std::endl;
}


// Top-level function to run both threads
void RunElGamalSocketDemo(const std::string& server_ip, int port) {
    ElGamalDecryptorSocket(server_ip, port);
}

int main()
{
    int I;
    // Lets test for location I = 1234
    // First cheat and look for the tag of the element at location I
    //I = 1234;

    //PIR_Experiment(I); // Test with the first element in the database

    //PIR_Experiment(0); // Test with the first element in the database
    //PIR_Experiment(49999); // Test with the last element in the database
    //PIR_Experiment(2555); // Test with the first element in the database
    //PIR_Experiment(3000); // Test with the first element in the database

    //ElGamalDecryptor();

    std::string server_ip = "127.0.0.1";
    int port = 8080;
    RunElGamalSocketDemo(server_ip, port);
    return 0;
}