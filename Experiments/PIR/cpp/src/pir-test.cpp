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
#include <gmpxx.h>
#include <random>

#define PROFILE
#define NUM_CPU_CORES 16

#include <fstream>
#include <iostream>
#include <iterator>

class ElGamal {
private:
    
public:
    mpz_class p, q, g;  // p: large prime, q: subgroup order, g: generator
    // Initialize with safe primes p, q where q|(p-1)
    void setup(int p_bits = 3072, int q_bits = 256) {
        // Generate q (256-bit prime)
        gmp_randclass rng(gmp_randinit_default);
        rng.seed(time(NULL));
        
        // First generate q
        do {
            // Generate a random prime q of q_bits size
            q = rng.get_z_bits(q_bits);
            // Get the next prime number
            mpz_nextprime(q.get_mpz_t(), q.get_mpz_t());
        } while (mpz_sizeinbase(q.get_mpz_t(), 2) != q_bits); // Ensure q is exactly q_bits long
        
        // Then generate p such that q|(p-1)
        mpz_class temp;
        do {
            temp = rng.get_z_bits(p_bits - q_bits);
            p = temp * q + 1;
        } while (!mpz_probab_prime_p(p.get_mpz_t(), 25));// Repetition count of 25 for the primality test
        
        // Find generator g of order q
        mpz_class h, exp;
        exp = (p - 1) / q;
        
        // TODO: Verify the generator generation logic
        do {
            h = rng.get_z_range(p - 1) + 1;
            mpz_powm(g.get_mpz_t(), h.get_mpz_t(), exp.get_mpz_t(), p.get_mpz_t());
        } while (g == 1);
    }
    
    // Generate random element in subgroup of order q
    mpz_class randomGroupElement() {
        gmp_randclass rng(gmp_randinit_default);
        mpz_class r = rng.get_z_range(q);
        mpz_class result;
        mpz_powm(result.get_mpz_t(), g.get_mpz_t(), r.get_mpz_t(), p.get_mpz_t());
        return result;
    }
    
    // ElGamal key generation
    pair<mpz_class, mpz_class> keyGen() {
        gmp_randclass rng(gmp_randinit_default);
        mpz_class x = rng.get_z_range(q);  // private key
        mpz_class y;
        mpz_powm(y.get_mpz_t(), g.get_mpz_t(), x.get_mpz_t(), p.get_mpz_t());  // public key
        return make_pair(x, y);
    }
    
    // ElGamal encryption
    pair<mpz_class, mpz_class> encrypt(const mpz_class& message, const mpz_class& publicKey) {
        gmp_randclass rng(gmp_randinit_default);
        mpz_class k = rng.get_z_range(q);  // random exponent
        
        mpz_class c1, c2;
        mpz_powm(c1.get_mpz_t(), g.get_mpz_t(), k.get_mpz_t(), p.get_mpz_t());
        
        mpz_class temp;
        mpz_powm(temp.get_mpz_t(), publicKey.get_mpz_t(), k.get_mpz_t(), p.get_mpz_t());
        c2 = (message * temp) % p;
        
        return make_pair(c1, c2);
    }
    
    // ElGamal decryption
    mpz_class decrypt(const pair<mpz_class, mpz_class>& ciphertext, const mpz_class& privateKey) {
        mpz_class c1 = ciphertext.first;
        mpz_class c2 = ciphertext.second;
        
        mpz_class temp, inv_temp;
        mpz_powm(temp.get_mpz_t(), c1.get_mpz_t(), privateKey.get_mpz_t(), p.get_mpz_t());
        mpz_invert(inv_temp.get_mpz_t(), temp.get_mpz_t(), p.get_mpz_t());
        
        return (c2 * inv_temp) % p;
    }

    // ElGamal multiplication of ciphertexts
    pair<mpz_class, mpz_class> mult_ct(const pair<mpz_class, mpz_class>& ciphertext1, const pair<mpz_class, mpz_class>& ciphertext2) {
        mpz_class cm1 = (ciphertext1.first*ciphertext2.first) % p;
        mpz_class cm2 = (ciphertext1.second*ciphertext2.second) % p;
        
        return make_pair(cm1, cm2);
    }
};

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

    //TODO: Our scheme has the benefit that the entire FSS lookup can be achieved without any disk access
    //Since the entire stash may reside in the RAM
    // Initialize the database
    for(size_t i=0; i<N; i++) {
        mpz_class entry_val;
        mpz_ui_pow_ui(entry_val.get_mpz_t(), 2, 496);//A 496-bit multiplier 
        DB[i] = mpz_class(i)*(entry_val - 1); // And the value stored at DB[i] is i*(2^496) is a 512 bit sized value
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

void TestElGamal() {
    using namespace std::chrono;
    std::cout << "=== ElGamal Test ===" << std::endl;
    ElGamal elgamal;
    int p_bits = 3072, q_bits = 256;

    auto t_start = high_resolution_clock::now();
    elgamal.setup(p_bits, q_bits);
    auto t_setup = high_resolution_clock::now();
    std::cout << "Setup complete." << std::endl;

    // Print parameters
    std::cout << "p (prime, " << p_bits << " bits): " << elgamal.p.get_str() << std::endl;
    std::cout << "q (subgroup order, " << q_bits << " bits): " << elgamal.q.get_str() << std::endl;
    std::cout << "g (generator): " << elgamal.g.get_str() << std::endl;

    auto t_params = high_resolution_clock::now();

    // Key generation
    auto [priv, pub] = elgamal.keyGen();
    auto t_keygen = high_resolution_clock::now();
    std::cout << "Private key x: " << priv.get_str() << std::endl;
    std::cout << "Public key y: " << pub.get_str() << std::endl;

    // Choose a valid random message which belongs to the subgroup G
    gmp_randclass rng(gmp_randinit_default);
    rng.seed(time(NULL));
    mpz_class exp = rng.get_z_range(elgamal.q);
    mpz_class m1;
    mpz_powm(m1.get_mpz_t(), elgamal.g.get_mpz_t(), exp.get_mpz_t(), elgamal.p.get_mpz_t());

    // Generate another random message
    exp = rng.get_z_range(elgamal.q);
    mpz_class m2;
    mpz_powm(m2.get_mpz_t(), elgamal.g.get_mpz_t(), exp.get_mpz_t(), elgamal.p.get_mpz_t());

    auto t_message = high_resolution_clock::now();
    std::cout << "Message to encrypt: " << m1.get_str() << std::endl;

    // Encrypt m1
    auto [c11, c12] = elgamal.encrypt(m1, pub);
    std::cout << "Ciphertext c1: " << c11.get_str() << std::endl;
    std::cout << "Ciphertext c2: " << c12.get_str() << std::endl;

    // Encrypt m1
    auto [c21, c22] = elgamal.encrypt(m2, pub);
    auto t_encrypt = high_resolution_clock::now();
    std::cout << "Ciphertext c1: " << c21.get_str() << std::endl;
    std::cout << "Ciphertext c2: " << c22.get_str() << std::endl;

    //Multiply the ciphertexts
    auto [cm1, cm2] = elgamal.mult_ct({c11,c12}, {c21, c22});
    auto t_mul = high_resolution_clock::now();

    // Decrypt m1
    mpz_class decrypted = elgamal.decrypt({c11, c12}, priv);
    std::cout << "Decrypted message 1: " << decrypted.get_str() << std::endl;

    // Check
    if (decrypted == m1)
        std::cout << "ElGamal Test: Success! Decrypted message matches original message 1." << std::endl;
    else
        std::cout << "ElGamal Test: Failure! Decrypted message does not match." << std::endl;

    // Decrypt m2
    decrypted = elgamal.decrypt({c21, c22}, priv);
    std::cout << "Decrypted message 2: " << decrypted.get_str() << std::endl;

    // Check
    if (decrypted == m2)
        std::cout << "ElGamal Test: Success! Decrypted message matches original message 2." << std::endl;
    else
        std::cout << "ElGamal Test: Failure! Decrypted message does not match." << std::endl;

    // Decrypt {cm1, cm2}
    decrypted = elgamal.decrypt({cm1, cm2}, priv);
    auto t_decrypt = high_resolution_clock::now();
    std::cout << "Decrypted message from multiplied ciphertexts is: " << decrypted.get_str() << std::endl;

    // Check
    if (decrypted == ((m1*m2) % elgamal.p))
        std::cout << "ElGamal Test: Success! Decrypted message matches m1*m2." << std::endl;
    else
        std::cout << "ElGamal Test: Failure! Decrypted message does not match." << std::endl;

    // Print timing info
    std::cout << "Time for setup: " << duration<double, std::milli>(t_setup - t_start).count() << " ms" << std::endl;
    std::cout << "Time for printing parameters: " << duration<double, std::milli>(t_params - t_setup).count() << " ms" << std::endl;
    std::cout << "Time for key generation: " << duration<double, std::milli>(t_keygen - t_params).count() << " ms" << std::endl;
    std::cout << "Time for message generation: " << duration<double, std::milli>(t_message - t_keygen).count() << " ms" << std::endl;
    std::cout << "Time for encryption: " << duration<double, std::milli>(t_encrypt - t_message).count()/2 << " ms" << std::endl;
    std::cout << "Time for decryption: " << duration<double, std::milli>(t_decrypt - t_mul).count()/3 << " ms" << std::endl;
    std::cout << "Time for ciphertext multiplication: " << duration<double, std::milli>(t_mul - t_encrypt).count()/2 << " ms" << std::endl;
    std::cout << "Total time: " << duration<double, std::milli>(t_decrypt - t_start).count() << " ms" << std::endl;

    return;
}

int main()
{
    //PIR_Experiment(0); // Test with the first element in the database
    //PIR_Experiment(49999); // Test with the last element in the database

    TestElGamal();

    return 0;
}