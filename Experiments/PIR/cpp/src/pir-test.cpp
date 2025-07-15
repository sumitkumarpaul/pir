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

int TFHE_128_bit_Example(void)
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
    FheBool *selBit = NULL;

    auto t_FHEBoolEncStart = std::chrono::high_resolution_clock::now(); 
    // A 128-bit unsigned integer containing value: 20 << 64 | 10
    U128 clear_lhs = { .w0 = 10, .w1 = 20 };
    // A 128-bit unsigned integer containing value: 2 << 64 | 1
    U128 clear_rhs = { .w0 = 1, .w1 = 2 };
    
    // A boolean value to select the bit
    bool clear_selBit = true; // This will select the first bit of the subtraction result

    ok = fhe_bool_try_encrypt_with_client_key_bool(clear_selBit, client_key, &selBit);
    assert(ok == 0);
    auto t_FHEBoolEncEnd = std::chrono::high_resolution_clock::now();
    std::cout << "Time to encrypt selectBit: " <<
     std::chrono::duration<double, std::milli>(t_FHEBoolEncEnd - t_FHEBoolEncStart).count()
     << " ms" << endl;

    auto t_FHEEncStart = std::chrono::high_resolution_clock::now(); 

    ok = fhe_uint128_try_encrypt_with_client_key_u128(clear_lhs, client_key, &lhs);
    assert(ok == 0);
    auto t_FHEEncEnd = std::chrono::high_resolution_clock::now();
    std::cout << "Time to encrypt FHE lhs: " <<
     std::chrono::duration<double, std::milli>(t_FHEEncEnd - t_FHEEncStart).count()
     << " ms" << endl;

    ok = fhe_uint128_try_encrypt_with_client_key_u128(clear_rhs, client_key, &rhs);
    assert(ok == 0);

    ok = fhe_uint128_if_then_else(selBit, lhs, rhs, &result);
    assert(ok == 0);
    auto t_FHEIfThenElseEnd = std::chrono::high_resolution_clock::now();
    std::cout << "Time to compute FHE if-then-else: " <<
     std::chrono::duration<double, std::milli>(t_FHEIfThenElseEnd - t_FHEEncEnd).count()
     << " ms" << endl;

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

uint8_t shortint_homomorphic_and(uint8_t a, uint8_t b) {
    // Use a predefined parameter set
    //const ShortintPBSParameters params = SHORTINT_V1_2_PARAM_MESSAGE_6_CARRY_1_KS_PBS_GAUSSIAN_2M128;
    const ShortintPBSParameters params = {
        .lwe_dimension = 2,
        .glwe_dimension = 1,
        .polynomial_size = 1,
        .lwe_noise_distribution = {0},
        .glwe_noise_distribution = {0},
        .pbs_base_log = 1,
        .pbs_level = 1,
        .ks_base_log = 1,
        .ks_level = 1,
        .message_modulus = 1,
        .carry_modulus = 1,
        .max_noise_level = 1,
        .log2_p_fail = 1.0,
        .modulus_power_of_2_exponent = 1,
        .encryption_key_choice = ShortintEncryptionKeyChoiceSmall,
        .modulus_switch_noise_reduction_params = {0},
        },
    };
    ShortintClientKey* client_key = nullptr;
    ShortintPublicKey* public_key = nullptr;
    ShortintServerKey* server_key = nullptr;
    ShortintCiphertext *ct_a = nullptr, *ct_b = nullptr, *ct_res = nullptr;
    uint8_t result = 0;
    int err = 0;
    uint64_t tmp = 0;
    auto t0 = std::chrono::high_resolution_clock::now();
    auto t1 = t0, t2 = t0, t3 = t0, t4 = t0, t5 = t0;

    // Key generation
    err = shortint_gen_client_key(params, &client_key);
    if (err) { std::cerr << "Client key gen failed\n"; goto cleanup; }
    err = shortint_gen_public_key(client_key, &public_key);
    if (err) { std::cerr << "Public key gen failed\n"; goto cleanup; }
    err = shortint_gen_server_key(client_key, &server_key);
    if (err) { std::cerr << "Server key gen failed\n"; goto cleanup; }

    // Encrypt
    t0 = std::chrono::high_resolution_clock::now();
    err = shortint_public_key_encrypt(public_key, a, &ct_a);
    err |= shortint_public_key_encrypt(public_key, b, &ct_b);
    t1 = std::chrono::high_resolution_clock::now();
    if (err) { std::cerr << "Encryption failed\n"; goto cleanup; }

    // Homomorphic AND (unchecked)
    t2 = std::chrono::high_resolution_clock::now();
    err = shortint_server_key_unchecked_bitand(server_key, ct_a, ct_b, &ct_res);
    t3 = std::chrono::high_resolution_clock::now();
    if (err) { std::cerr << "Homomorphic AND failed\n"; goto cleanup; }

    // Decrypt
    t4 = std::chrono::high_resolution_clock::now();
    err = shortint_client_key_decrypt(client_key, ct_res, &tmp);
    t5 = std::chrono::high_resolution_clock::now();
    if (err) { std::cerr << "Decryption failed\n"; goto cleanup; }

    std::cout << "Encryption time: " << std::chrono::duration<double, std::milli>(t1-t0).count() << " ms\n";
    std::cout << "Homomorphic AND time: " << std::chrono::duration<double, std::milli>(t3-t2).count() << " ms\n";
    std::cout << "Decryption time: " << std::chrono::duration<double, std::milli>(t5-t4).count() << " ms\n";
    std::cout << "Result: " << tmp << std::endl;

cleanup:
    shortint_destroy_ciphertext(ct_a);
    shortint_destroy_ciphertext(ct_b);
    shortint_destroy_ciphertext(ct_res);
    shortint_destroy_public_key(public_key);
    shortint_destroy_client_key(client_key);
    shortint_destroy_server_key(server_key);
    return result;
}

#if 1

int shortint_mul_demo() {
    int err = 0;

    // 1. Generate parameters and keys
    struct ShortintClientKey *client_key = NULL;
    struct ShortintServerKey *server_key = NULL;

    // Use a parameter set (adjust as needed)
    err = shortint_gen_keys_with_parameters(
        //SHORTINT_V1_1_PARAM_MESSAGE_4_CARRY_4_KS_PBS_GAUSSIAN_2M128,
        SHORTINT_V1_2_PARAM_MESSAGE_4_CARRY_1_KS_PBS_GAUSSIAN_2M128,
        &client_key,
        &server_key
    );
    if (err) { printf("Key generation failed\n"); return err; }

    // 2. Encrypt two shortints
    uint64_t a = 2, b = 7;
    
    struct ShortintCiphertext *ct_a = NULL, *ct_b = NULL;
    err = shortint_client_key_encrypt(client_key, a, &ct_a);
    if (err) { printf("Encrypt a with client key failed\n"); }
    err = shortint_client_key_encrypt(client_key, b, &ct_b);
    if (err) { printf("Encrypt b with client key failed\n");  }

    /*
    struct ShortintCiphertext *ct_a_srv = NULL, *ct_b_srv = NULL;
    err = shortint_server_key_encrypt(server_key, a, &ct_a_srv);
    if (err) { printf("Encrypt a with server key failed\n"); }
    err = shortint_server_key_encrypt(server_key, b, &ct_b);
    if (err) { printf("Encrypt b with server key failed\n");  }
    */
    
    // 3. Multiply (unchecked)
    struct ShortintCiphertext *ct_unchecked = NULL;
    auto t_UncheckedMulStart = std::chrono::high_resolution_clock::now();
    err = shortint_server_key_unchecked_mul(server_key, ct_a, ct_b, &ct_unchecked);
    auto t_UncheckedMulEnd = std::chrono::high_resolution_clock::now();
    if (err) { printf("Unchecked mul failed\n"); }
    std::cout << "Time to multiply with unchecked: "
        << std::chrono::duration<double, std::milli>(t_UncheckedMulEnd - t_UncheckedMulStart).count()
        << " ms" << std::endl;

    auto t_UncheckedAddStart = std::chrono::high_resolution_clock::now();
    err = shortint_server_key_unchecked_add(server_key, ct_a, ct_b, &ct_unchecked);
    auto t_UncheckedAddEnd = std::chrono::high_resolution_clock::now();
    if (err) { printf("Unchecked add failed\n"); }
    std::cout << "Time to add with unchecked: "
        << std::chrono::duration<double, std::milli>(t_UncheckedAddEnd - t_UncheckedAddStart).count()
        << " ms" << std::endl;


    // 4. Multiply (smart/checked)
    struct ShortintCiphertext *ct_smart = NULL;
    auto t_SmartMulStart = std::chrono::high_resolution_clock::now();
    err = shortint_server_key_smart_mul(server_key, ct_a, ct_b, &ct_smart);
    auto t_SmartMulEnd = std::chrono::high_resolution_clock::now();
    if (err) { printf("Smart mul failed\n");  }
    std::cout << "Time to multiply with smart: "
        << std::chrono::duration<double, std::milli>(t_SmartMulEnd - t_SmartMulStart).count()
        << " ms" << std::endl;

    // 5. Decrypt and print results
    uint64_t res_unchecked = 0, res_smart = 0;
    shortint_client_key_decrypt(client_key, ct_unchecked, &res_unchecked);
    shortint_client_key_decrypt(client_key, ct_smart, &res_smart);
    printf("Unchecked mul result: %lu\n", res_unchecked);
    printf("Smart mul result: %lu\n", res_smart);

    // 6. Cleanup
    shortint_destroy_ciphertext(ct_smart);
cleanup_unchecked:
    shortint_destroy_ciphertext(ct_unchecked);
cleanup_ct_b:
    shortint_destroy_ciphertext(ct_b);
cleanup_ct_a:
    shortint_destroy_ciphertext(ct_a);
cleanup_keys:
    shortint_destroy_client_key(client_key);
    shortint_destroy_server_key(server_key);

    return err;
}

#endif

int TFHE_64_bit_Example(void)
{
    int ok = 0;
    // Prepare the config builder for the high level API and choose which types to enable
    ConfigBuilder *builder;
    Config *config;

    // Put the builder in a default state without any types enabled
    config_builder_default(&builder);
    
    // Use the custom parameters 
    config_builder_use_custom_parameters(&builder, SHORTINT_V0_11_PARAM_MESSAGE_2_CARRY_2_KS_PBS_GAUSSIAN_2M64);

    //use_dedicated_compact_public_key_parameters(&builder, SHORTINT_PARAM_PKE_MESSAGE_2_CARRY_2_KS_PBS_TUNIFORM_2M128);

    // Populate the config
    config_builder_build(builder, &config);

    ClientKey *client_key = NULL;
    ServerKey *server_key = NULL;

    // Generate the keys using the config
    generate_keys(config, &client_key, &server_key);
    // Set the server key for the current thread
    set_server_key(server_key);

    FheUint64 *lhs = NULL;
    FheUint64 *rhs = NULL;
    FheUint64 *result = NULL;
    FheBool *selBit = NULL;

    auto t_FHEBoolEncStart = std::chrono::high_resolution_clock::now(); 

    uint64_t clear_lhs = 10;

    uint64_t clear_rhs = 20;
    
    // A boolean value to select the bit
    bool clear_selBit = true; // This will select the first bit of the subtraction result

    ok = fhe_bool_try_encrypt_with_client_key_bool(clear_selBit, client_key, &selBit);
    assert(ok == 0);
    auto t_FHEBoolEncEnd = std::chrono::high_resolution_clock::now();
    std::cout << "Time to encrypt selectBit: " <<
     std::chrono::duration<double, std::milli>(t_FHEBoolEncEnd - t_FHEBoolEncStart).count()
     << " ms" << endl;

    auto t_FHEEncStart = std::chrono::high_resolution_clock::now(); 

    ok = fhe_uint64_try_encrypt_with_client_key_u64(clear_lhs, client_key, &lhs);
    assert(ok == 0);
    auto t_FHEEncEnd = std::chrono::high_resolution_clock::now();
    std::cout << "Time to encrypt FHE lhs: " <<
     std::chrono::duration<double, std::milli>(t_FHEEncEnd - t_FHEEncStart).count()
     << " ms" << endl;

    ok = fhe_uint64_try_encrypt_with_client_key_u64(clear_rhs, client_key, &rhs);
    assert(ok == 0);

    ok = fhe_uint64_if_then_else(selBit, lhs, rhs, &result);
    assert(ok == 0);
    auto t_FHEIfThenElseEnd = std::chrono::high_resolution_clock::now();
    std::cout << "Time to compute FHE if-then-else: " <<
     std::chrono::duration<double, std::milli>(t_FHEIfThenElseEnd - t_FHEEncEnd).count()
     << " ms" << endl;

    uint64_t clear_result;

    auto t_FHEDecStart = std::chrono::high_resolution_clock::now();    
    // Decrypt
    ok = fhe_uint64_decrypt(result, client_key, &clear_result);
    assert(ok == 0);
    auto t_FHEDecEnd = std::chrono::high_resolution_clock::now();
    std::cout << "The decryption result is:" << clear_result << ", Time to compute FHE decryption: " <<
     std::chrono::duration<double, std::milli>(t_FHEDecEnd - t_FHEDecStart).count()
     << " ms" << endl;

    // Here the subtraction allows us to compare each word
    assert(clear_result == clear_lhs);

    // Destroy the ciphertexts
    fhe_uint64_destroy(lhs);
    fhe_uint64_destroy(rhs);
    fhe_uint64_destroy(result);
    fhe_bool_destroy(selBit);

    // Destroy the keys
    client_key_destroy(client_key);
    server_key_destroy(server_key);

    printf("64-bit FHE computation successful!\n");
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

    //PIR_Experiment(N-1);

    //TFHE_64_bit_Example();

    shortint_homomorphic_and(63,63);

    //shortint_mul_demo();

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
