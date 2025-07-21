#include <chrono>
#include <assert.h>
#include <stdio.h>
#include "fss-common.h"
#include "fss-server.h"
#include "fss-client.h"
#include <immintrin.h>  // Include AVX header
#include <omp.h>
#include <vector>
#include <gmp.h>

#define PROFILE
#define NUM_CPU_CORES 16

#include <fstream>
#include <iostream>
#include <iterator>

#define N 50000 // Size of the database, can be adjusted as needed
// The value of N will determine the bitlength during the client initialization
#define num_bits_in_N 16 // 16 bits can represent up to 65536, which is more than enough for N=50000
// And number of bits determine the evalution time drastically
static mpz_class DB[N];

// Client asks for the element located at location i, where 0 <= i < N 
int PIR_Experiment(int I)
{
    // Set up variables
    Fss fClient, fServer;
    ServerKeyEq k0;
    ServerKeyEq k1;

    // Initialize the database
    for(size_t i=0; i<N; i++) {
        mpz_class entry_val;
        mpz_ui_pow_ui(entry_val.get_mpz_t(), 2, 496);
        DB[i] = mpz_class(i)*(entry_val); // Each location stores a multiple of 10 of its index
    }  

    auto t_begin = std::chrono::high_resolution_clock::now();

    // Initialize client, use 32 bits in domain as example
    initializeClient(&fClient, num_bits_in_N, 2); // If bit length is not set properly, then incorrect answer will be returned
    // bit length must be able to store ans0 and ans1
    
    auto t_initClient = std::chrono::high_resolution_clock::now();

    // Equality FSS test
    generateTreeEq(&fClient, &k0, &k1, I, 1);//So that the point function will evaluate as 1 at location i, and zero elsewhere

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
                auto x = evaluateEq(&fServer, &k0, i + j);
                thread_sums[j] += x * mpz_class(DB[i + j]);
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
                auto x = evaluateEq(&fServer, &k1, i + j);
                thread_sums[j] += x * mpz_class(DB[i + j]);
            }
        }
        for (int t = 0; t < NUM_CPU_CORES; ++t)
            ans1 += thread_sums[t];
    }

    fin = ans0 - ans1;
    cout << "PIR Test: value of ans0: " << ans0 << ", value of ans1: " << ans1 << endl;
    cout << "PIR Test: (it should be the value stored at DB[I]): " << fin << endl;
    cout << "The value stored at DB[I] is: " << DB[I] << endl;

    auto t_evalServer1 = std::chrono::high_resolution_clock::now();    
    
    if (fin != mpz_class(DB[I])) {
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

    PIR_Experiment(20);

    return 0;
}