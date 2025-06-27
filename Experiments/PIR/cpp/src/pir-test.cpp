#include <chrono>
#include <assert.h>
#include <stdio.h>
#include "fss-common.h"
#include "fss-server.h"
#include "fss-client.h"
#include "tfhe.h"

#define N 10000 // Size of the database, can be adjusted as needed
// The value of N wiil determine the bitlength during the client initialization
#define BITS 14 
// And number of bits determine the evalution time drastically
static u_int64_t DB[N];

// Client asks for the element located at location i, where 0 <= i < N 
int PIR_Experiment(int I)
{
    // Set up variables
    Fss fClient, fServer;
    ServerKeyEq k0;
    ServerKeyEq k1;

    // Initialize the database
    for(size_t i=0; i<N; i++) {
        DB[i] = i*10000; // Each location stores a multiple of 10 of its index
    }  

    auto t_begin = std::chrono::high_resolution_clock::now();

    // Initialize client, use 32 bits in domain as example
    initializeClient(&fClient, BITS, 2); // If bit length is not set properly, then incorrect answer will be returned
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

    // Server 0, evalute this for all i in [0, N)
    for(size_t i=0; i<N; i++) {
        auto x = evaluateEq(&fServer, &k0, i);
        ans0 += x * mpz_class(DB[i]); // Accumulate the value at each index
    }

    auto t_evalServer0 = std::chrono::high_resolution_clock::now();

    // Server 1, evalute this for all i in [0, N)
    for(size_t i=0; i<N; i++) {
        auto x = evaluateEq(&fServer, &k1, i);
        ans1 += x * mpz_class(DB[i]); // Accumulate the value at each index
    }

    fin = ans0 - ans1;
    cout << "PIR Test: value of ans0: " << ans0 << ", value of ans1: " << ans1 << endl;
    cout << "PIR Test: (it should be the value stored at DB[I]): " << fin << endl;
    cout << "The value stored at DB[I] is: " << DB[I] << endl;

    auto t_evalServer1 = std::chrono::high_resolution_clock::now();    
    
    if (fin != mpz_class(DB[I])) {
        cout << "!!!!!!!!!!!!!!!!!!! PIR Test: Error, the value returned does not match the expected value !!!!!!!!!!!!!" << endl;
        return 0; // Error
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

int TFHE_Example(void)
{
    int ok = 0;
    // Prepare the config builder for the high level API and choose which types to enable
    ConfigBuilder *builder;
    Config *config;

    // Put the builder in a default state without any types enabled
    config_builder_default(&builder);
    // Populate the config
    config_builder_build(builder, &config);

    ClientKey *client_key = NULL;
    ServerKey *server_key = NULL;

    // Generate the keys using the config
    generate_keys(config, &client_key, &server_key);
    // Set the server key for the current thread
    set_server_key(server_key);

    FheUint128 *lhs = NULL;
    FheUint128 *rhs = NULL;
    FheUint128 *result = NULL;

    auto t_FHEEncStart = std::chrono::high_resolution_clock::now(); 
    // A 128-bit unsigned integer containing value: 20 << 64 | 10
    U128 clear_lhs = { .w0 = 10, .w1 = 20 };
    // A 128-bit unsigned integer containing value: 2 << 64 | 1
    U128 clear_rhs = { .w0 = 1, .w1 = 2 };

    ok = fhe_uint128_try_encrypt_with_client_key_u128(clear_lhs, client_key, &lhs);
    assert(ok == 0);
    auto t_FHEEncEnd = std::chrono::high_resolution_clock::now();
    std::cout << "Time to encrypt FHE lhs: " <<
     std::chrono::duration<double, std::milli>(t_FHEEncEnd - t_FHEEncStart).count()
     << " ms" << endl;

    ok = fhe_uint128_try_encrypt_with_client_key_u128(clear_rhs, client_key, &rhs);
    assert(ok == 0);

    auto t_FHESubStart = std::chrono::high_resolution_clock::now();
    // Compute the subtraction
    ok = fhe_uint128_sub(lhs, rhs, &result);
    assert(ok == 0);
    auto t_FHESubEnd = std::chrono::high_resolution_clock::now();
    std::cout << "Time to compute FHE subtraction: " <<
     std::chrono::duration<double, std::milli>(t_FHESubEnd - t_FHESubStart).count()
     << " ms" << endl;

    U128 clear_result;

    auto t_FHEDecStart = std::chrono::high_resolution_clock::now();    
    // Decrypt
    ok = fhe_uint128_decrypt(result, client_key, &clear_result);
    assert(ok == 0);
    auto t_FHEDecEnd = std::chrono::high_resolution_clock::now();
    std::cout << "Time to compute FHE decryption: " <<
     std::chrono::duration<double, std::milli>(t_FHEDecEnd - t_FHEDecStart).count()
     << " ms" << endl;

    // Here the subtraction allows us to compare each word
    assert(clear_result.w0 == 9);
    assert(clear_result.w1 == 18);

    // Destroy the ciphertexts
    fhe_uint128_destroy(lhs);
    fhe_uint128_destroy(rhs);
    fhe_uint128_destroy(result);

    // Destroy the keys
    client_key_destroy(client_key);
    server_key_destroy(server_key);

    printf("FHE computation successful!\n");
    return EXIT_SUCCESS;
}

int main()
{
    // Set up variables
    uint64_t a = 3;
    uint64_t b = 2;
    Fss fClient, fServer;
    ServerKeyEq k0;
    ServerKeyEq k1;

    PIR_Experiment(N-1);

    TFHE_Example();

    return 0;//Skip the rest of the tests for now

    // Initialize client, use 10 bits in domain as example
    initializeClient(&fClient, 10, 2); 
    
    // Equality FSS test
    generateTreeEq(&fClient, &k0, &k1, a, b);
    
    // Initialize server
    initializeServer(&fServer, &fClient);
    mpz_class ans0, ans1, fin;
    
    ans0 = evaluateEq(&fServer, &k0, a);
    ans1 = evaluateEq(&fServer, &k1, a);
    fin = ans0 - ans1;
    cout << "FSS Eq Match (should be non-zero): " << fin << endl;
    
    ans0 = evaluateEq(&fServer, &k0, (a-1));
    ans1 = evaluateEq(&fServer, &k1, (a-1));
    fin = ans0 - ans1;
    cout << "FSS Eq No Match (should be 0): " << fin << endl;

    // Less than FSS test
    ServerKeyLt lt_k0;
    ServerKeyLt lt_k1;
    
    initializeClient(&fClient, 10, 2);
    generateTreeLt(&fClient, &lt_k0, &lt_k1, a, b);

    initializeServer(&fServer, &fClient);
    uint64_t lt_ans0, lt_ans1, lt_fin;

    lt_ans0 = evaluateLt(&fServer, &lt_k0, (a-1));
    lt_ans1 = evaluateLt(&fServer, &lt_k1, (a-1));
    lt_fin = lt_ans0 - lt_ans1;
    cout << "FSS Lt Match (should be non-zero): " << lt_fin << endl;

    lt_ans0 = evaluateLt(&fServer, &lt_k0, a);
    lt_ans1 = evaluateLt(&fServer, &lt_k1, a);
    lt_fin = lt_ans0 - lt_ans1;
    cout << "FSS Lt No Match (should be 0): " << lt_fin << endl;

    // Equality FSS test for multi-parties
    MPKey mp_keys[3];
    initializeClient(&fClient, 10, 3);
    generateTreeEqMParty(&fClient, a, b, mp_keys);

    initializeServer(&fServer, &fClient);
    uint32_t mp_ans0 = evaluateEqMParty(&fServer, &mp_keys[0], a);
    uint32_t mp_ans1 = evaluateEqMParty(&fServer, &mp_keys[1], a);
    uint32_t mp_ans2 = evaluateEqMParty(&fServer, &mp_keys[2], a);
    uint32_t xor_mp = mp_ans0 ^ mp_ans1 ^ mp_ans2;
    cout << "FSS Eq Multi-Party Match (should be non-zero): " << xor_mp << endl;

    mp_ans0 = evaluateEqMParty(&fServer, &mp_keys[0], (a+1));
    mp_ans1 = evaluateEqMParty(&fServer, &mp_keys[1], (a+1));
    mp_ans2 = evaluateEqMParty(&fServer, &mp_keys[2], (a+1));
    xor_mp = mp_ans0 ^ mp_ans1 ^ mp_ans2;
    cout << "FSS Eq Multi-Party No Match (should be 0): " << xor_mp << endl;

    size_t rounds = 100000;
    auto t_begin = std::chrono::high_resolution_clock::now();
    for(size_t i=0; i<rounds; i++) {
        volatile auto x = evaluateEq(&fServer, &k0, i);
    }
    for(size_t i=0; i<rounds; i++) {
        volatile auto x = evaluateLt(&fServer, &lt_k0, i);
    }
    for(size_t i=0; i<rounds; i++) {
        volatile auto x = evaluateEqMParty(&fServer, &mp_keys[1], a);
    }
    auto t_end = std::chrono::high_resolution_clock::now();
    std::cout << "Benchmark result: " <<
     std::chrono::duration<double, std::milli>(t_end - t_begin).count()
     << " ms" << endl;
    return 1;
}
