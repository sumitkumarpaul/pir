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

#define PROFILE
#define NUM_CPU_CORES 16

#include <iterator>
#include <cstring>
#include <unistd.h>
#include "pir_common.h"


// 2. Deserializes and validates input, returns true if valid and fills params
bool ElGamalDeserializeAndValidate(const std::string& data, std::vector<std::string>& params) {
    params.clear();
    std::string temp = data;
    size_t pos = 0;
    while ((pos = temp.find("\n")) != std::string::npos) {
        params.push_back(temp.substr(0, pos));
        temp.erase(0, pos + 1);
    }
    if (params.size() < 4) {
        std::cerr << "Server: Invalid public key message (deserialization)" << std::endl;
        return false;
    }
    // Optionally, check if params are valid numbers
    for (int i = 0; i < 4; ++i) {
        if (params[i].empty()) {
            std::cerr << "Server: Empty parameter at index " << i << std::endl;
            return false;
        }
    }
    return true;
}

// 3. Performs ElGamal operations and sends result on the same socket
void ElGamalPerformAndSend(const std::vector<std::string>& params, int client_socket) {
    std::string p_str = params[0], q_str = params[1], g_str = params[2], pub_str = params[3];
    ElGamal elgamal;
    elgamal.p = mpz_class(p_str);
    elgamal.q = mpz_class(q_str);
    elgamal.g = mpz_class(g_str);
    mpz_class pub(pub_str);
    mpz_class message = elgamal.randomGroupElement();
    auto [c1, c2] = elgamal.encrypt(message, pub);
    std::string msg = c1.get_str() + "\n" + c2.get_str() + "\n" + message.get_str() + "\n";
    send(client_socket, msg.c_str(), msg.size(), 0);
}

void ElGamalEncryptor(int port) {
    int ret = -1;
    int server_fd, new_socket;
    char buffer[65536] = {0};

    ret = InitAcceptingSocket(port, &server_fd, &new_socket);

    if (ret != 0) {
        std::cerr << "Server: Failed to initialize accepting socket" << std::endl;
        return;
    }

    int valread = recv(new_socket, buffer, sizeof(buffer), 0);
    if (valread <= 0) {
        std::cerr << "Server: Read failed" << std::endl;
        close(new_socket); close(server_fd);
        return;
    }

    std::string data(buffer, valread);
    std::vector<std::string> params;
    if (!ElGamalDeserializeAndValidate(data, params)) {
        std::cerr << "Server: Deserialization error, aborting." << std::endl;
        close(new_socket); close(server_fd);
        return;
    }
    ElGamalPerformAndSend(params, new_socket);

    FinishAcceptingSocket(server_fd, new_socket);
}


#define N 50000 // Size of the database, can be adjusted as needed
// The value of N will determine the bitlength during the client initialization
#define NUM_TAG_BITS 3072 // 16 bits can represent up to 65536, which is more than enough for N=50000
#define B 512 // Block size in bits, can be adjusted as needed
// And number of bits determine the evalution time drastically
static mpz_class sh[N][2]; // Database to store values, each entry is a pair, {Tag, Block-content}.

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
        sh[i][0] = r.get_z_bits(NUM_TAG_BITS); // Create a random tag
        sh[i][1] = r.get_z_bits(B); // Create a random block of B bits
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
                mpz_class a = sh[i+j][0]; // Get the tag of the element at location i
                mpz_mul(a.get_mpz_t(), a.get_mpz_t(), del_abc.get_mpz_t()); // Perform multiplication operation
                mpz_mod(sh[i][0].get_mpz_t(), a.get_mpz_t(), p.get_mpz_t()); // Perform modular operation
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
    mpz_class T_sh = sh[I][0]; // The tag of the element at location I 
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
                auto y = evaluateEq(&fServer, &k0, sh[i + j][0]);// Evaluate the FSS for the tag
                thread_sums[j] += y * sh[i + j][1]; // Multiply the result with the block content
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
                auto y = evaluateEq(&fServer, &k1, sh[i + j][0]);// Evaluate the FSS for the tag
                thread_sums[j] += y * sh[i + j][1]; // Multiply the result with the block content
            }
        }
        for (int t = 0; t < NUM_CPU_CORES; ++t)
            ans1 += thread_sums[t];
    }

    fin = ans0 - ans1;
    cout << "PIR Test: value of ans0: " << ans0 << ", value of ans1: " << ans1 << endl;
    cout << "PIR Test: (it should be the value stored at DB[I]): " << fin << endl;
    cout << "The value stored at DB[I] is: " << sh[I][1] << endl;

    auto t_evalServer1 = std::chrono::high_resolution_clock::now();    
    
    if (fin != mpz_class(sh[I][1])) {
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

int main()
{
    //int I;
    // Lets test for location I = 1234
    // First cheat and look for the tag of the element at location I
    //I = 1234;

    //PIR_Experiment(I); // Test with the first element in the database

    //PIR_Experiment(0); // Test with the first element in the database
    //PIR_Experiment(49999); // Test with the last element in the database
    //PIR_Experiment(2555); // Test with the first element in the database
    //PIR_Experiment(3000); // Test with the first element in the database

    //TestElGamal();

    ElGamalEncryptor(8080);

    return 0;
}