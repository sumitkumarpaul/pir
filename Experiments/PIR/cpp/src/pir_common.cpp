#include "pir_common.h"


#define TIC1(t)    t = std::chrono::high_resolution_clock::now()
#define TOC1(t)    std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - t).count()

// Randomization
gmp_randclass rng(gmp_randinit_default);

// Global ElGamal parameters
mpz_class p, q, r, g, g_q;

// El-Gamal encryption keys
mpz_class pk_E, pk_E_q;
mpz_class sk_E, sk_E_q;

// FHE related
PublicKey<DCRTPoly> pk_F;
PrivateKey<DCRTPoly> sk_F;
CryptoContext<DCRTPoly> FHEcryptoContext;
Ciphertext<DCRTPoly> vectorOnesforElement_ct;
Ciphertext<DCRTPoly> vectorOnesforTag_ct;
Ciphertext<DCRTPoly> fnd_ct;

std::string ready_for_epoch_message = "READY_FOR_EPOCH";

void PrintLog(int log_level, const char* file, int line, const std::string& message) {

    if (log_level <= SET_LOG_LEVEL) {
        const char* LogTypeStr;

        if (log_level == LOG_LEVEL_ERROR) {
            LogTypeStr = "ERROR: ";
        } else if (log_level == LOG_LEVEL_INFO) {
            LogTypeStr = "INFO: ";
        } else if (log_level == LOG_LEVEL_DEBUG) {
            LogTypeStr= "DEBUG: ";
        } else if (log_level == LOG_LEVEL_TRACE) {
            LogTypeStr= "TRACE: ";
        } else {
            LogTypeStr= "***** ";
        }
        // Get current time
        auto now = std::chrono::system_clock::now();
        auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(now.time_since_epoch()) % 1000;
        std::time_t t_now = std::chrono::system_clock::to_time_t(now);
        std::tm tm_now;
#ifdef _WIN32
        localtime_s(&tm_now, &t_now);
#else
        localtime_r(&t_now, &tm_now);
#endif
        char time_buf[32];
        std::strftime(time_buf, sizeof(time_buf), "%d-%m-%Y %H:%M:%S", &tm_now);
        const char *const_file = (strrchr(file, '/') ? strrchr(file, '/') + 1 : file);
        std::cout << "[" << time_buf << ":" << std::setfill('0') << std::setw(3) << ms.count() << "] [" << const_file << ":" << line << "] " << LogTypeStr << message << std::endl;
    }
}

// El-Gamal related functions
mpz_class ElGamal_randomGroupElement() {
    mpz_class localr = rng.get_z_range(q);
    mpz_class result;
    mpz_powm(result.get_mpz_t(), g.get_mpz_t(), localr.get_mpz_t(), p.get_mpz_t());

    return result;
}

std::pair<mpz_class, mpz_class> ElGamal_keyGen() {
    //Randomness is already initialized during the initialization of the servers
    mpz_class x = rng.get_z_range(q);
    mpz_class y;
    mpz_powm(y.get_mpz_t(), g.get_mpz_t(), x.get_mpz_t(), p.get_mpz_t());

    return std::make_pair(y, x);
}

std::pair<mpz_class, mpz_class> ElGamal_q_keyGen() {//q and g_q are global parameters and set previously
    //Randomness is already initialized during the initialization of the servers
    mpz_class x = rng.get_z_range(q-1)+1;//i.e., within ZZ_q*
    mpz_class y;
    mpz_powm(y.get_mpz_t(), g_q.get_mpz_t(), x.get_mpz_t(), q.get_mpz_t());

    return std::make_pair(y, x);
}

std::pair<mpz_class, mpz_class> ElGamal_encrypt(const mpz_class& message, const mpz_class& publicKey) {
    mpz_class k = rng.get_z_range(q);
    mpz_class c1, c2;
    mpz_powm(c1.get_mpz_t(), g.get_mpz_t(), k.get_mpz_t(), p.get_mpz_t());
    mpz_class temp;
    mpz_powm(temp.get_mpz_t(), publicKey.get_mpz_t(), k.get_mpz_t(), p.get_mpz_t());
    c2 = (message * temp) % p;
    return std::make_pair(c1, c2);
}

std::pair<mpz_class, mpz_class> ElGamal_q_encrypt(const mpz_class& message, const mpz_class& publicKey) {
    mpz_class k = rng.get_z_range(q-1)+1;//Deliberately choosing it in ZZ_q*, instead of ZZ_q
    mpz_class c1, c2;
    mpz_powm(c1.get_mpz_t(), g_q.get_mpz_t(), k.get_mpz_t(), q.get_mpz_t());
    mpz_class temp;
    mpz_powm(temp.get_mpz_t(), publicKey.get_mpz_t(), k.get_mpz_t(), q.get_mpz_t());
    c2 = (message * temp) % q;
    return std::make_pair(c1, c2);
}

mpz_class ElGamal_decrypt(const std::pair<mpz_class, mpz_class>& ciphertext, const mpz_class& privateKey) {
    mpz_class c1 = ciphertext.first;
    mpz_class c2 = ciphertext.second;
    mpz_class temp, inv_temp;
    mpz_powm(temp.get_mpz_t(), c1.get_mpz_t(), privateKey.get_mpz_t(), p.get_mpz_t());
    mpz_invert(inv_temp.get_mpz_t(), temp.get_mpz_t(), p.get_mpz_t());
    return (c2 * inv_temp) % p;
}

mpz_class ElGamal_q_decrypt(const std::pair<mpz_class, mpz_class>& ciphertext, const mpz_class& privateKey) {
    mpz_class c1 = ciphertext.first;
    mpz_class c2 = ciphertext.second;
    mpz_class temp, inv_temp;
    mpz_powm(temp.get_mpz_t(), c1.get_mpz_t(), privateKey.get_mpz_t(), q.get_mpz_t());
    mpz_invert(inv_temp.get_mpz_t(), temp.get_mpz_t(), q.get_mpz_t());
    return (c2 * inv_temp) % q;
}

std::pair<mpz_class, mpz_class> ElGamal_mult_ct(const std::pair<mpz_class, mpz_class>& ciphertext1, const std::pair<mpz_class, mpz_class>& ciphertext2) {
    mpz_class cm1 = (ciphertext1.first * ciphertext2.first) % p;
    mpz_class cm2 = (ciphertext1.second * ciphertext2.second) % p;
    return std::make_pair(cm1, cm2);
}

std::pair<mpz_class, mpz_class> ElGamal_q_mult_ct(const std::pair<mpz_class, mpz_class>& ciphertext1, const std::pair<mpz_class, mpz_class>& ciphertext2) {
    mpz_class cm1 = (ciphertext1.first * ciphertext2.first) % q;
    mpz_class cm2 = (ciphertext1.second * ciphertext2.second) % q;
    return std::make_pair(cm1, cm2);
}

std::pair<mpz_class, mpz_class> ElGamal_exp_ct(const std::pair<mpz_class, mpz_class>& ciphertext, const mpz_class& exp, const mpz_class& publicKey) {
    mpz_class c1, c2;
    mpz_powm(c1.get_mpz_t(), ciphertext.first.get_mpz_t(), exp.get_mpz_t(), p.get_mpz_t());
    mpz_powm(c2.get_mpz_t(), ciphertext.second.get_mpz_t(), exp.get_mpz_t(), p.get_mpz_t());

#if 0 /* For the time being donot multiply with ciphertext of 1, since anyway our protocol will multiply with E(h_C) */
    auto [cI1, cI2] = ElGamal_encrypt(mpz_class(1), publicKey);
    return ElGamal_mult_ct({c1, c2}, {cI1, cI2});
#else
    return std::make_pair(c1, c2);
#endif
}

std::pair<mpz_class, mpz_class> ElGamal_q_exp_ct(const std::pair<mpz_class, mpz_class>& ciphertext, const mpz_class& exp, const mpz_class& publicKey) {
    mpz_class c1, c2;
    mpz_powm(c1.get_mpz_t(), ciphertext.first.get_mpz_t(), exp.get_mpz_t(), q.get_mpz_t());
    mpz_powm(c2.get_mpz_t(), ciphertext.second.get_mpz_t(), exp.get_mpz_t(), q.get_mpz_t());

#if 0 /* For the time being donot multiply with ciphertext of 1, since anyway our protocol will multiply with E(h_C) */
    auto [cI1, cI2] = ElGamal_encrypt(mpz_class(1), publicKey);
    return ElGamal_mult_ct({c1, c2}, {cI1, cI2});
#else
    return std::make_pair(c1, c2);
#endif
}
// Networking related functions
int InitAcceptingSocket(int port, int* p_server_fd, int* p_new_socket) {
    struct sockaddr_in address;
    int opt = 1;
    int addrlen = sizeof(address);
    *p_new_socket = -1;
    *p_server_fd = -1;
    int ret = -1;

    *p_server_fd = socket(AF_INET, SOCK_STREAM, 0);
    if (*p_server_fd < 0) {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "socket() returned: " + std::to_string(*p_server_fd));
        return ret;
    }
    if (setsockopt(*p_server_fd, SOL_SOCKET, SO_REUSEADDR | SO_REUSEPORT, &opt, sizeof(opt)) < 0) {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "setsockopt() failed!!");
        close(*p_server_fd);
        *p_server_fd = -1;
        return ret;
    }
    address.sin_family = AF_INET;
    address.sin_addr.s_addr = INADDR_ANY;
    address.sin_port = htons(port);
    if (bind(*p_server_fd, (struct sockaddr *)&address, sizeof(address)) < 0) {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "bind() failed!!");
        close(*p_server_fd);
        *p_server_fd = -1;
        return ret;
    }
    if (listen(*p_server_fd, 1) < 0) {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "listen() returned: " + std::to_string(*p_new_socket));
        close(*p_server_fd);
        *p_server_fd = -1;
        return ret;
    }

    PrintLog(LOG_LEVEL_SPECIAL, __FILE__, __LINE__, "Listening on port: " + std::to_string(port));

    *p_new_socket = accept(*p_server_fd, (struct sockaddr *)&address, (socklen_t*)&addrlen);
    if (*p_new_socket < 0) {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "accept() returned: " + std::to_string(*p_new_socket));
        close(*p_server_fd);
        *p_server_fd = -1;
        return ret;
    } else {
        char client_ip[INET_ADDRSTRLEN];
        inet_ntop(AF_INET, &(address.sin_addr), client_ip, INET_ADDRSTRLEN);
        int client_port = ntohs(address.sin_port);
        PrintLog(LOG_LEVEL_INFO, __FILE__, __LINE__, "Accepted connection from IP: " + std::string(client_ip) + ", Port: " + std::to_string(client_port));
    }

    return 0;
}

void FinishAcceptingSocket(int server_fd, int new_socket) {
    close(new_socket);
    close(server_fd);
}

void InitConnectingSocket(const std::string& server_ip, int port, int* p_sock) {
    // Connect to server, blocking call
    struct sockaddr_in serv_addr;
    *p_sock = socket(AF_INET, SOCK_STREAM, 0);
    serv_addr.sin_family = AF_INET;
    serv_addr.sin_port = htons(port);
    inet_pton(AF_INET, server_ip.c_str(), &serv_addr.sin_addr);

    while (connect(*p_sock, (struct sockaddr *)&serv_addr, sizeof(serv_addr)) < 0) {
        std::this_thread::sleep_for(std::chrono::seconds(1)); // Wait before retrying
    }

    return;
}

// Send all the data
int sendAll(int sock, const char* data, size_t sz) {
    size_t totalSent = 0;
    char buf[20];// The size must not be larger than that

    // First send the size to the receiving end
    if (send(sock, std::to_string(sz).c_str(), std::to_string(sz).size(), 0) < std::to_string(sz).size()) {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Failed to sending data size");
        return -1;
    }

    // The receiving end should echo back with the same message
    ssize_t received = recv(sock, buf, sizeof(buf), 0);
    if ((received <= 0) || (std::string(buf, received) != std::to_string(sz))) {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Failed to receive echo from the other end");
        return -1;
    }

    // Now start sending the entire message
    while (totalSent < sz) {
        ssize_t sent = send(sock, data + totalSent, sz - totalSent, 0);
        if (sent <= 0) {
            PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "After sending " + std::to_string(totalSent) + " bytes, send() returned: " + std::to_string(sent));
            return -1; // Error occurred
        }
        totalSent += sent;
    }

    return 0;
}

// Receive all the data
int recvAll(int sock, char* data, size_t max_sz, size_t* received_sz) {
    char buf[20];// The size must not be larger than that
    size_t totalReceived = 0;

    // Step 1: Receive the size
    ssize_t sz_len = recv(sock, buf, sizeof(buf), 0);
    if (sz_len <= 0) {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Failed to receive data size");
        return -1;
    }

    std::string sz_str(buf, sz_len);
    size_t sz = std::stoull(sz_str);

    if (sz > max_sz) {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Received size " + std::to_string(sz) + " exceeds maximum allowed size " + std::to_string(max_sz));
        return -1; // Size exceeds maximum allowed
    }

    // Step 2: Echo the size back
    if (send(sock, sz_str.c_str(), sz_str.size(), 0) < sz_str.size()) {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "Failed to echo data size");
        return -1;
    }

    // Step 3: Receive the full message
    while (totalReceived < sz) {
        ssize_t recvd = recv(sock, data + totalReceived, (sz - totalReceived), 0);
        if (recvd <= 0) {
            PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "After receiving " + std::to_string(totalReceived) + " bytes, recv() returned: " + std::to_string(recvd));
            return -1;
        }
        totalReceived += recvd;
    }

    if (received_sz != NULL) {
        *received_sz = totalReceived;
    }

    return 0;
}

// FHE related functions
int FHE_keyGen(){
    ////////////////////////////////////////////////////////////
    // Set-up of parameters
    ////////////////////////////////////////////////////////////
    KeyPair<DCRTPoly> keyPair;

    // Crypto Parameters
    // # of evalMults = 3 (first 3) is used to support the multiplication of 7
    // ciphertexts, i.e., ceiling{log2{7}} Max depth is set to 3 (second 3) to
    // generate homomorphic evaluation multiplication keys for s^2 and s^3
    CCParams<CryptoContextBGVRNS> parameters;
    parameters.SetMultiplicativeDepth(1);
    parameters.SetPlaintextModulus(65537);//TODO, 65537, 536903681 these values must have special properties.

    /*****************************************************************
     * Since, the plaintext modulus is 65537, hence upto 16-bit number
     * can be represented in a single ciphertext. However, we may add
     * two ciphertexts as well and the result must be within 16-bit.
     * Hence, each individual plaintext must remain within 15-bit.
     * ***************************************************************/

    //At this moment, using the value mentioned in the original example.
    parameters.SetMaxRelinSkDeg(1);// Initially 1 What does this value mean?
    parameters.SetScalingTechnique(FIXEDAUTO);//Only this is not giving any exception and giving good result
    
    //parameters.SetSecurityLevel(HEStd_128_classic);
    //parameters.SetRingDim(8192);

    FHEcryptoContext = GenCryptoContext(parameters);
    // enable features that you wish to use
    FHEcryptoContext->Enable(PKE);
    //FHEcryptoContext->Enable(KEYSWITCH);
    FHEcryptoContext->Enable(LEVELEDSHE);
    //FHEcryptoContext->Enable(ADVANCEDSHE);

    // Initialize Public Key Containers
    keyPair = FHEcryptoContext->KeyGen();

    if (!keyPair.good()) {
        PrintLog(LOG_LEVEL_ERROR, __FILE__, __LINE__, "FHE Key generation failed!");
        return -1;
    } else {
        pk_F = keyPair.publicKey;
        sk_F = keyPair.secretKey;
    }

    return 0;
}

Ciphertext<DCRTPoly> FHE_Enc_DBElement(const mpz_class block_content, const mpz_class block_index) {
    std::vector<int64_t> DBElementVector;/* FHE can encrypt 15-bits. But we must use int64_t vector, since this is what the existing function takes */
    
    /* Add the block content */
    mpz_class tmp = block_content;
    for (unsigned i = 0; i < NUM_FHE_BLOCKS_PER_PIR_BLOCK; ++i) {
        mpz_class rem;
        mpz_fdiv_r_2exp(rem.get_mpz_t(), tmp.get_mpz_t(), PLAINTEXT_FHE_BLOCK_SIZE); // rem = tmp % 2^15
        DBElementVector.push_back(static_cast<int64_t>(rem.get_ui()));
        mpz_fdiv_q_2exp(tmp.get_mpz_t(), tmp.get_mpz_t(), PLAINTEXT_FHE_BLOCK_SIZE); // tmp >>= 15
    }

    /* Append the index */
    tmp = block_index;
    for (unsigned i = 0; i < NUM_FHE_BLOCKS_PER_PIR_INDEX; ++i) {
        mpz_class rem;
        mpz_fdiv_r_2exp(rem.get_mpz_t(), tmp.get_mpz_t(), PLAINTEXT_FHE_BLOCK_SIZE); // rem = tmp % 2^15
        DBElementVector.push_back(static_cast<int64_t>(rem.get_ui()));
        mpz_fdiv_q_2exp(tmp.get_mpz_t(), tmp.get_mpz_t(), PLAINTEXT_FHE_BLOCK_SIZE); // tmp >>= 15
    }

    Plaintext FHEPackedPlaintext = FHEcryptoContext->MakePackedPlaintext(DBElementVector);

    /* TODO: If problem occurs, check whether it is due to this compression or not */
    //return FHEcryptoContext->Compress(FHEcryptoContext->Encrypt(pk_F, FHEPackedPlaintext), 1);

    /* Let's not compress at this moment, since it might take extra time, as well as can break functionality */
    return FHEcryptoContext->Encrypt(pk_F, FHEPackedPlaintext);
}

// Decrypts a ciphertext, unpacks the packed vector, and reconstructs block_content and block_index as mpz_class
void FHE_Dec_DBElement(const Ciphertext<DCRTPoly>& ct, mpz_class& block_content, mpz_class& block_index) {
    Plaintext pt;
    FHEcryptoContext->Decrypt(sk_F, ct, &pt);
    pt->SetLength(TOTAL_NUM_FHE_BLOCKS_PER_ELEMENT);
    const std::vector<int64_t>& packed = pt->GetPackedValue();

    // Reconstruct block_content from first NUM_FHE_BLOCKS_PER_PIR_BLOCK elements
    block_content = 0;
    for (int i = NUM_FHE_BLOCKS_PER_PIR_BLOCK - 1; i >= 0; --i) {
        block_content <<= PLAINTEXT_FHE_BLOCK_SIZE;
        block_content += packed[i] & ((1 << PLAINTEXT_FHE_BLOCK_SIZE) - 1);
    }

    // Reconstruct block_index from next NUM_FHE_BLOCKS_PER_PIR_INDEX elements
    block_index = 0;
    for (int i = TOTAL_NUM_FHE_BLOCKS_PER_ELEMENT - 1; i >= NUM_FHE_BLOCKS_PER_PIR_BLOCK; --i) {
        block_index <<= PLAINTEXT_FHE_BLOCK_SIZE;
        block_index += packed[i] & ((1 << PLAINTEXT_FHE_BLOCK_SIZE) - 1);
    }
}

// Encrypt the tag
Ciphertext<DCRTPoly> FHE_Enc_Tag(const mpz_class tag) {
    std::vector<int64_t> TagVector;/* FHE can encrypt 15-bits. But we must use int64_t vector, since this is what the existing function takes */
    
    mpz_class tmp = tag;
    for (unsigned i = 0; i < NUM_FHE_BLOCKS_PER_TAG; ++i) {
        mpz_class rem;
        mpz_fdiv_r_2exp(rem.get_mpz_t(), tmp.get_mpz_t(), PLAINTEXT_FHE_BLOCK_SIZE); // rem = tmp % 2^15
        TagVector.push_back(static_cast<int64_t>(rem.get_ui()));
        mpz_fdiv_q_2exp(tmp.get_mpz_t(), tmp.get_mpz_t(), PLAINTEXT_FHE_BLOCK_SIZE); // tmp >>= 15
    }

    Plaintext FHEPackedPlaintext = FHEcryptoContext->MakePackedPlaintext(TagVector);

    return FHEcryptoContext->Encrypt(pk_F, FHEPackedPlaintext);
}

// Decrypt the tag
void FHE_Dec_Tag(const Ciphertext<DCRTPoly>& ct, mpz_class& tag) {
    Plaintext pt;
    FHEcryptoContext->Decrypt(sk_F, ct, &pt);
    pt->SetLength(NUM_FHE_BLOCKS_PER_TAG);
    const std::vector<int64_t>& packed = pt->GetPackedValue();

    // Reconstruct block_content from first NUM_FHE_BLOCKS_PER_TAG elements
    tag = 0;
    for (int i = NUM_FHE_BLOCKS_PER_TAG - 1; i >= 0; --i) {
        tag <<= PLAINTEXT_FHE_BLOCK_SIZE;
        tag += packed[i] & ((1 << PLAINTEXT_FHE_BLOCK_SIZE) - 1);
    }
}

// selElementBits_ct must be encryption of select bit but extended over TOTAL_NUM_FHE_BLOCKS_PER_ELEMENT
Ciphertext<DCRTPoly> FHE_SelectElement(const Ciphertext<DCRTPoly>& selElementBits_ct, const Ciphertext<DCRTPoly>& A_ct, const Ciphertext<DCRTPoly>& B_ct){
    ////////////////////////////////////////////////////////////
    // Homomorphic selection between two ciphertexts w/o any relinearization
    // Select(a,b, select_bit) = select_bit ? a : b
    //                         = ((1 - select_bit) * a) + (select_bit * b)
    ////////////////////////////////////////////////////////////
    
    auto fnd_not_ct = FHEcryptoContext->EvalSub(vectorOnesforElement_ct, selElementBits_ct);
    FHEcryptoContext->ModReduceInPlace(fnd_not_ct);

    auto A_not_ct = FHEcryptoContext->EvalMultNoRelin(fnd_not_ct, A_ct);
    FHEcryptoContext->ModReduceInPlace(A_not_ct);

    auto B_mul_fnd_ct = FHEcryptoContext->EvalMultNoRelin(selElementBits_ct, B_ct);
    FHEcryptoContext->ModReduceInPlace(B_mul_fnd_ct);

    auto C_ct = FHEcryptoContext->EvalAdd(A_not_ct, B_mul_fnd_ct);
    FHEcryptoContext->ModReduceInPlace(C_ct);

    return C_ct;
}

// selectTagBits_ct must be encryption of select bit but extended over NUM_FHE_BLOCKS_PER_TAG
Ciphertext<DCRTPoly> FHE_SelectTag(const Ciphertext<DCRTPoly>& selectTagBits_ct, const Ciphertext<DCRTPoly>& A_ct, const Ciphertext<DCRTPoly>& B_ct){
    ////////////////////////////////////////////////////////////
    // Homomorphic selection between two ciphertexts w/o any relinearization
    // Select(a,b, select_bit) = select_bit ? a : b
    //                         = ((1 - select_bit) * a) + (select_bit * b)
    ////////////////////////////////////////////////////////////
    
    auto fnd_not_ct = FHEcryptoContext->EvalSub(vectorOnesforTag_ct, selectTagBits_ct);
    FHEcryptoContext->ModReduceInPlace(fnd_not_ct);

    auto A_not_ct = FHEcryptoContext->EvalMultNoRelin(fnd_not_ct, A_ct);
    FHEcryptoContext->ModReduceInPlace(A_not_ct);

    auto B_mul_fnd_ct = FHEcryptoContext->EvalMultNoRelin(selectTagBits_ct, B_ct);
    FHEcryptoContext->ModReduceInPlace(B_mul_fnd_ct);

    auto C_ct = FHEcryptoContext->EvalAdd(A_not_ct, B_mul_fnd_ct);
    FHEcryptoContext->ModReduceInPlace(C_ct);

    return C_ct;
}

void FHE_EncOfOnes(Ciphertext<DCRTPoly>& OnesforElement_ct, Ciphertext<DCRTPoly>& OnesforTag_ct){
    std::vector<int64_t> vectorOfOnes;
    Plaintext plaintextOnes;

    // Use a for loop to add elements to the vector
    for (int i = 0; i < TOTAL_NUM_FHE_BLOCKS_PER_ELEMENT; ++i) {
        vectorOfOnes.push_back(1);
    }
    plaintextOnes = FHEcryptoContext->MakePackedPlaintext(vectorOfOnes);
    OnesforElement_ct = FHEcryptoContext->Encrypt(pk_F, plaintextOnes);

    vectorOfOnes.clear();
    // Use a for loop to add elements to the vector
    for (int i = 0; i < NUM_FHE_BLOCKS_PER_TAG; ++i) {
        vectorOfOnes.push_back(1);
    }
    
    plaintextOnes = FHEcryptoContext->MakePackedPlaintext(vectorOfOnes);
    OnesforTag_ct = FHEcryptoContext->Encrypt(pk_F, plaintextOnes);

    return;
}

// bytes: most-significant byte first (big-endian)
mpz_class import_from_bytes(const std::string &bytes) {
    PrintLog(LOG_LEVEL_DEBUG, __FILE__, __LINE__, "Size of the serialized ciphertext: " + std::to_string(bytes.size()) + " bytes");
    mpz_t tmp;
    mpz_init(tmp);
    // count = bytes.size(), size=1 (bytes), order=1 (most-significant word first),
    // endian = 1 (most-significant byte first), nails = 0
    mpz_import(tmp, bytes.size(), 1, 1, 1, 0, bytes.data());
    mpz_class out(tmp);
    mpz_clear(tmp);
    return out;
}

mpz_class serialized_ct_to_mpz_class(const std::string& filename) {
//void serialized_ct_to_mpz_class(const std::string& filename) {
    std::ifstream in(filename, std::ios::binary);
    if (!in) throw std::runtime_error("Cannot open file");

    std::vector<unsigned char> buf(
        (std::istreambuf_iterator<char>(in)),
        std::istreambuf_iterator<char>()
    );

    mpz_t tmp;
    mpz_init(tmp);
    // Interpret the buffer as big-endian bytes:
    mpz_import(tmp, buf.size(), 1, 1, 1, 0, buf.data());

    mpz_class result(tmp);
    mpz_clear(tmp);
    in.close();
    return result;
}

// Serializes an mpz_class to a file as big-endian bytes
void mpz_class_to_serialized_ct(const mpz_class& value, const std::string& filename) {
    std::ofstream out(filename, std::ios::binary);
    if (!out) throw std::runtime_error("Cannot open file for writing");

    // Export mpz_class to big-endian bytes
    size_t count = 0;
    void* data = mpz_export(nullptr, &count, 1, 1, 1, 0, value.get_mpz_t());
    if (data && count > 0) {
        out.write(reinterpret_cast<const char*>(data), count);
        free(data);
    }
    out.close();
}

void insert_pdb_entry(std::fstream& pdb, uint64_t id, const plain_db_entry& entry) {
    //TODO: Error check
    pdb.seekp(static_cast<std::streampos>(id) * sizeof(plain_db_entry));
    pdb.write(reinterpret_cast<const char*>(&entry), sizeof(plain_db_entry));

    return;
}

void read_pdb_entry(std::fstream& pdb, uint64_t id, plain_db_entry& out_entry) {
    //TODO: Error check
    pdb.seekg(static_cast<std::streampos>(id) * sizeof(plain_db_entry));
    pdb.read(reinterpret_cast<char*>(&out_entry), sizeof(plain_db_entry));
    return;
}

void insert_sdb_entry(std::fstream& sdb, uint64_t id, const shuffled_db_entry& entry) {
    //TODO: Error check
    sdb.seekp(static_cast<std::streampos>(id) * sizeof(shuffled_db_entry));
    sdb.write(reinterpret_cast<const char*>(&entry), sizeof(shuffled_db_entry));

    return;
}

void read_sdb_entry(std::fstream& sdb, uint64_t id, shuffled_db_entry& out_entry) {
    //TODO: Error check
    sdb.seekg(static_cast<std::streampos>(id) * sizeof(shuffled_db_entry));
    sdb.read(reinterpret_cast<char*>(&out_entry), sizeof(shuffled_db_entry));
    return;
}
