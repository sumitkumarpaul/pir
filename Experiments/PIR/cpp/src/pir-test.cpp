#include <chrono>
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

#include <fstream>
#include <iostream>
#include <iterator>
#include <thread>
#include <cstring>
#include <sys/socket.h>
#include <netinet/in.h>
#include <arpa/inet.h>
#include <unistd.h>

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

    // ElGamal exponentiation of ciphertexts
    pair<mpz_class, mpz_class> exp_ct(const pair<mpz_class, mpz_class>& ciphertext, const mpz_class& exp, const mpz_class& publicKey) {
        mpz_class c1, c2;

        // Exponentiate both the components of the ciphertexts
        mpz_powm(c1.get_mpz_t(), ciphertext.first.get_mpz_t(), exp.get_mpz_t(), p.get_mpz_t());       
        mpz_powm(c2.get_mpz_t(), ciphertext.second.get_mpz_t(), exp.get_mpz_t(), p.get_mpz_t());

        //TODO: Only for security purpose
        // Generate a ciphertext of 1 with new randomness, which will make the scheme IND-CPA secure
        auto [cI1, cI2] = encrypt(mpz_class(1), publicKey); // Encrypt 1 with public key

        // Return the exponentiated ciphertext and the multiplcation of the ciphertext of 1
        return mult_ct({c1, c2}, {cI1, cI2}); // Multiply the ciphertexts
    }
};

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



// Encryptor: acts as server, receives public key, sends ciphertext and message
void ElGamalEncryptorSocket(int port) {
    int server_fd, new_socket;
    struct sockaddr_in address;
    int opt = 1;
    int addrlen = sizeof(address);
    char buffer[65536] = {0};

    // Create socket file descriptor
    server_fd = socket(AF_INET, SOCK_STREAM, 0);
    setsockopt(server_fd, SOL_SOCKET, SO_REUSEADDR | SO_REUSEPORT, &opt, sizeof(opt));
    address.sin_family = AF_INET;
    address.sin_addr.s_addr = INADDR_ANY;
    address.sin_port = htons(port);
    bind(server_fd, (struct sockaddr *)&address, sizeof(address));
    listen(server_fd, 1);
    new_socket = accept(server_fd, (struct sockaddr *)&address, (socklen_t*)&addrlen);

    // Receive public key and params
    int valread = read(new_socket, buffer, sizeof(buffer));
    std::string received(buffer, valread);
    // Expect: p\nq\ng\npub\n
    size_t pos = 0;
    std::vector<std::string> params;
    while ((pos = received.find("\n")) != std::string::npos) {
        params.push_back(received.substr(0, pos));
        received.erase(0, pos + 1);
    }
    if (params.size() < 4) {
        std::cerr << "Server: Invalid public key message" << std::endl;
        close(new_socket); close(server_fd); return;
    }
    std::string p_str = params[0], q_str = params[1], g_str = params[2], pub_str = params[3];
    ElGamal elgamal;
    elgamal.p = mpz_class(p_str);
    elgamal.q = mpz_class(q_str);
    elgamal.g = mpz_class(g_str);
    mpz_class pub(pub_str);
    // Generate random message
    mpz_class message = elgamal.randomGroupElement();
    auto [c1, c2] = elgamal.encrypt(message, pub);
    std::string msg = c1.get_str() + "\n" + c2.get_str() + "\n" + message.get_str() + "\n";
    send(new_socket, msg.c_str(), msg.size(), 0);
    close(new_socket); close(server_fd);
}


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

    // Connect to server
    int sock = 0;
    struct sockaddr_in serv_addr;
    char buffer[65536] = {0};
    sock = socket(AF_INET, SOCK_STREAM, 0);
    serv_addr.sin_family = AF_INET;
    serv_addr.sin_port = htons(port);
    inet_pton(AF_INET, server_ip.c_str(), &serv_addr.sin_addr);
    while (connect(sock, (struct sockaddr *)&serv_addr, sizeof(serv_addr)) < 0) {
        std::this_thread::sleep_for(std::chrono::milliseconds(100)); // Wait for server
    }
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
    std::thread server_thread([port]() {
        ElGamalEncryptorSocket(port);
    });
    std::this_thread::sleep_for(std::chrono::milliseconds(200)); // Ensure server starts first
    std::thread client_thread([server_ip, port]() {
        ElGamalDecryptorSocket(server_ip, port);
    });
    client_thread.join();
    server_thread.join();
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