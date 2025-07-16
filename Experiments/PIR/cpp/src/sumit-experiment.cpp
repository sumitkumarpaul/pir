//==================================================================================
// BSD 2-Clause License
//
// Copyright (c) 2014-2022, NJIT, Duality Technologies Inc. and other contributors
//
// All rights reserved.
//
// Author TPOC: contact@openfhe.org
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
// 1. Redistributions of source code must retain the above copyright notice, this
//    list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright notice,
//    this list of conditions and the following disclaimer in the documentation
//    and/or other materials provided with the distribution.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//==================================================================================

/*
  Example of a computation circuit of depth 3
  BGVrns demo for a homomorphic multiplication of depth 6 and three different approaches for depth-3 multiplications
 */

#define PROFILE

#include <chrono>
#include <fstream>
#include <iostream>
#include <iterator>

#include "openfhe.h"

using namespace lbcrypto;

int OpenFHEBGVrns_example(){
    ////////////////////////////////////////////////////////////
    // Set-up of parameters
    ////////////////////////////////////////////////////////////

    // benchmarking variables
    TimeVar t;
    double processingTime(0.0);

    // Crypto Parameters
    // # of evalMults = 3 (first 3) is used to support the multiplication of 7
    // ciphertexts, i.e., ceiling{log2{7}} Max depth is set to 3 (second 3) to
    // generate homomorphic evaluation multiplication keys for s^2 and s^3
    CCParams<CryptoContextBGVRNS> parameters;
    parameters.SetMultiplicativeDepth(1);
    parameters.SetPlaintextModulus(536903681);//TODO, this value must have special properties.
    //At this moment, using the value mentioned in the original example.
    parameters.SetMaxRelinSkDeg(2);//What does this value mean?

    CryptoContext<DCRTPoly> cryptoContext = GenCryptoContext(parameters);
    // enable features that you wish to use
    cryptoContext->Enable(PKE);
    //cryptoContext->Enable(KEYSWITCH);
    cryptoContext->Enable(LEVELEDSHE);
    //cryptoContext->Enable(ADVANCEDSHE);

    std::cout << "\np = " << cryptoContext->GetCryptoParameters()->GetPlaintextModulus() << std::endl;
    std::cout << "n = " << cryptoContext->GetCryptoParameters()->GetElementParams()->GetCyclotomicOrder() / 2
              << std::endl;
    std::cout << "log2 q = "
              << log2(cryptoContext->GetCryptoParameters()->GetElementParams()->GetModulus().ConvertToDouble())
              << std::endl;

    // Initialize Public Key Containers
    KeyPair<DCRTPoly> keyPair;

    // Perform Key Generation Operation
    TIC(t);

    keyPair = cryptoContext->KeyGen();

    processingTime = TOC(t);
    std::cout << "Key generation time: " << processingTime << "ms" << std::endl;

    if (!keyPair.good()) {
        std::cout << "Key generation failed!" << std::endl;
        exit(1);
    }

    ////////////////////////////////////////////////////////////
    // Encode source data
    ////////////////////////////////////////////////////////////
    //Since the plaintext modulus (536903681) is a 30-bit number.
    // So the multiplication result can go upto 30-bits
    // For safety(TODO) we are keeping each number in 14-bit range.
    int min = 0x2FFF;
    int max = 0x3FFF;
    int range = max - min + 1;
    int vector_sz = 16384; // Size of the vector to be generated. TODO: Maximum usable is: 16384

    std::vector<int64_t> vectorOfPt1;
    // Use a for loop to add elements to the vector
    for (int i = 0; i < vector_sz; ++i) {
        int num = rand() % range + min;
        vectorOfPt1.push_back(num);
    }

    Plaintext plaintext1 = cryptoContext->MakePackedPlaintext(vectorOfPt1);

    std::vector<int64_t> vectorOfPt2;
    // Use a for loop to add elements to the vector
    for (int i = 0; i < vector_sz; ++i) {
        int num = rand() % range + min;
        vectorOfPt2.push_back(num);
    }
    Plaintext plaintext2 = cryptoContext->MakePackedPlaintext(vectorOfPt2);

    #if 0
    std::cout << "\nOriginal Plaintext #1: \n";
    std::cout << plaintext1 << std::endl;

    std::cout << "\nOriginal Plaintext #2: \n";
    std::cout << plaintext2 << std::endl;
    #endif

    ////////////////////////////////////////////////////////////
    // Encryption
    ////////////////////////////////////////////////////////////

    //std::cout << "\nRunning encryption of all plaintexts... ";

    std::vector<Ciphertext<DCRTPoly>> ciphertexts;

    TIC(t);

    ciphertexts.push_back(cryptoContext->Encrypt(keyPair.publicKey, plaintext1));
    ciphertexts.push_back(cryptoContext->Encrypt(keyPair.publicKey, plaintext2));

    processingTime = TOC(t);

    //std::cout << "Completed\n";

    std::cout << "\nAverage encryption time: " << processingTime / 2 << "ms" << std::endl;


    ////////////////////////////////////////////////////////////
    // Homomorphic multiplication of two ciphertexts w/o any relinearization
    ////////////////////////////////////////////////////////////

    //std::cout << "\nRunning a multiplication of two ciphertexts w/o relinearization...";
    
    TIC(t);

    auto ciphertextMult12 = cryptoContext->EvalMultNoRelin(ciphertexts[0], ciphertexts[1]);
    cryptoContext->ModReduceInPlace(ciphertextMult12);

    processingTime = TOC(t);
    std::cout << "Time for multiplying two ciphertexts w/o relinearization: " << processingTime << "ms" << std::endl;

    //std::cout << "Completed\n";  

    Plaintext plaintextDecMult12;

    TOC(t);
    cryptoContext->Decrypt(keyPair.secretKey, ciphertextMult12, &plaintextDecMult12);
    processingTime = TOC(t);
    std::cout << "Decryption time: " << processingTime << "ms" << std::endl;

    plaintextDecMult12->SetLength(plaintext1->GetLength());

    //std::cout << "\nResult of homomorphic multiplication of ciphertexts #1 and #2: \n";
    //TODO: I found evan multiplying 16384-14bit numbers, are done in 5ms..!!

    //std::cout << plaintextDecMult12 << std::endl;

    for (int i = 0; i < vector_sz; ++i) {
        if (plaintextDecMult12->GetPackedValue()[i] !=
            (vectorOfPt1[i] * vectorOfPt2[i])) {
            std::cout << "Error in multiplication of ciphertexts #1 and #2 at index " << i << std::endl;
            std::cout << "Expected: " << (vectorOfPt1[i] * vectorOfPt2[i]) << ", got: "
                      << plaintextDecMult12->GetPackedValue()[i] << std::endl;
            return 1;
        }
    }    

    return 0;
}

int OpenFHEBGVrns_select(int select_bit){
    ////////////////////////////////////////////////////////////
    // Set-up of parameters
    ////////////////////////////////////////////////////////////

    // benchmarking variables
    TimeVar t;
    double processingTime(0.0);

    // Crypto Parameters
    // # of evalMults = 3 (first 3) is used to support the multiplication of 7
    // ciphertexts, i.e., ceiling{log2{7}} Max depth is set to 3 (second 3) to
    // generate homomorphic evaluation multiplication keys for s^2 and s^3
    CCParams<CryptoContextBGVRNS> parameters;
    parameters.SetMultiplicativeDepth(1);
    parameters.SetPlaintextModulus(65537);//TODO, 536903681 this value must have special properties.
    //At this moment, using the value mentioned in the original example.
    parameters.SetMaxRelinSkDeg(2);//What does this value mean?

    CryptoContext<DCRTPoly> cryptoContext = GenCryptoContext(parameters);
    // enable features that you wish to use
    cryptoContext->Enable(PKE);
    //cryptoContext->Enable(KEYSWITCH);
    cryptoContext->Enable(LEVELEDSHE);
    //cryptoContext->Enable(ADVANCEDSHE);

    std::cout << "\np = " << cryptoContext->GetCryptoParameters()->GetPlaintextModulus() << std::endl;
    std::cout << "n = " << cryptoContext->GetCryptoParameters()->GetElementParams()->GetCyclotomicOrder() / 2
              << std::endl;
    std::cout << "log2 q = "
              << log2(cryptoContext->GetCryptoParameters()->GetElementParams()->GetModulus().ConvertToDouble())
              << std::endl;

    // Initialize Public Key Containers
    KeyPair<DCRTPoly> keyPair;

    // Perform Key Generation Operation
    TIC(t);

    keyPair = cryptoContext->KeyGen();

    processingTime = TOC(t);
    std::cout << "Key generation time: " << processingTime << "ms" << std::endl;

    if (!keyPair.good()) {
        std::cout << "Key generation failed!" << std::endl;
        exit(1);
    }

    ////////////////////////////////////////////////////////////
    // Encode source data
    ////////////////////////////////////////////////////////////
    //Since the plaintext modulus (536903681) is a 30-bit number.
    // So the multiplication result can go upto 30-bits
    // For safety(TODO) we are keeping each number in 14-bit range.
    //int min = 0x2FFF;
    //int max = 0x3FFF;
    //int range = max - min + 1;
    int vector_sz = 2; // Size of the vector to be generated. TODO: Maximum usable is: 16384

    std::vector<int64_t> vectorOfPt1;
    // Use a for loop to add elements to the vector
    for (int i = 0; i < vector_sz; ++i) {
        //int num = rand() % range + min;
        int num = 0;
        vectorOfPt1.push_back(num);
    }

    Plaintext plaintext1 = cryptoContext->MakePackedPlaintext(vectorOfPt1);

    std::vector<int64_t> vectorOfPt2;
    // Use a for loop to add elements to the vector
    for (int i = 0; i < vector_sz; ++i) {
        //int num = rand() % range + min;
        int num = 32767;
        vectorOfPt2.push_back(num);
    }
    Plaintext plaintext2 = cryptoContext->MakePackedPlaintext(vectorOfPt2);

    std::vector<int64_t> vectorOfSelect;
    // Use a for loop to add elements to the vector
    for (int i = 0; i < vector_sz; ++i) {
        vectorOfSelect.push_back(select_bit);
    }
    Plaintext plaintextZeros = cryptoContext->MakePackedPlaintext(vectorOfSelect);

    #if 1
    std::cout << "\nOriginal Plaintext #1: \n";
    std::cout << plaintext1 << std::endl;

    std::cout << "\nOriginal Plaintext #2: \n";
    std::cout << plaintext2 << std::endl;
    #endif

    ////////////////////////////////////////////////////////////
    // Encryption
    ////////////////////////////////////////////////////////////

    //std::cout << "\nRunning encryption of all plaintexts... ";

    std::vector<Ciphertext<DCRTPoly>> ciphertexts;

    TIC(t);

    ciphertexts.push_back(cryptoContext->Encrypt(keyPair.publicKey, plaintext1));
    ciphertexts.push_back(cryptoContext->Encrypt(keyPair.publicKey, plaintext2));
    ciphertexts.push_back(cryptoContext->Encrypt(keyPair.publicKey, plaintextZeros));

    processingTime = TOC(t);

    //std::cout << "Completed\n";

    std::cout << "\nAverage encryption time: " << processingTime / 3 << "ms" << std::endl;


    ////////////////////////////////////////////////////////////
    // Homomorphic selection between two ciphertexts w/o any relinearization
    ////////////////////////////////////////////////////////////

    //std::cout << "\nRunning a multiplication of two ciphertexts w/o relinearization...";
    
    TIC(t);

    auto ciphertextMult12 = cryptoContext->EvalMultNoRelin(ciphertexts[0], ciphertexts[1]);
    cryptoContext->ModReduceInPlace(ciphertextMult12);

    processingTime = TOC(t);
    std::cout << "Time for multiplying two ciphertexts w/o relinearization: " << processingTime << "ms" << std::endl;

    TIC(t);

    auto ciphertextAdd12 = cryptoContext->EvalAdd(ciphertexts[0], ciphertexts[1]);
    cryptoContext->ModReduceInPlace(ciphertextAdd12);

    processingTime = TOC(t);
    std::cout << "Time for adding two ciphertexts: " << processingTime << "ms" << std::endl;

    TIC(t);

    auto ciphertextSub12 = cryptoContext->EvalSub(ciphertexts[0], ciphertexts[1]);
    cryptoContext->ModReduceInPlace(ciphertextSub12);

    processingTime = TOC(t);
    std::cout << "Time for subtracting two ciphertexts: " << processingTime << "ms" << std::endl;

    //std::cout << "Completed\n";  

    Plaintext plaintextDecMult12;

    TOC(t);
    cryptoContext->Decrypt(keyPair.secretKey, ciphertextMult12, &plaintextDecMult12);
    processingTime = TOC(t);
    std::cout << "Decryption time: " << processingTime << "ms" << std::endl;

    plaintextDecMult12->SetLength(plaintext1->GetLength());

    std::cout << "\nResult of homomorphic multiplication of ciphertexts #1 and #2: \n";
    //TODO: I found evan multiplying 16384-14bit numbers, are done in 5ms..!!

    std::cout << plaintextDecMult12 << std::endl;

    for (int i = 0; i < vector_sz; ++i) {
        if (plaintextDecMult12->GetPackedValue()[i] !=
            (vectorOfPt1[i] * vectorOfPt2[i])) {
            std::cout << "Error in multiplication of ciphertexts #1 and #2 at index " << i << std::endl;
            std::cout << "Expected: " << (vectorOfPt1[i] * vectorOfPt2[i]) << ", got: "
                      << plaintextDecMult12->GetPackedValue()[i] << std::endl;
            return 1;
        }
    }    

    //Decrypt addition result and check correctness

    Plaintext plaintextDecAdd12;

    TOC(t);
    cryptoContext->Decrypt(keyPair.secretKey, ciphertextAdd12, &plaintextDecAdd12);
    processingTime = TOC(t);
    std::cout << "Decryption time: " << processingTime << "ms" << std::endl;

    plaintextDecAdd12->SetLength(plaintext1->GetLength());

    std::cout << "\nResult of homomorphic addition of ciphertexts #1 and #2: \n";
    //TODO: I found evan multiplying 16384-14bit numbers, are done in 5ms..!!

    std::cout << plaintextDecAdd12 << std::endl;

    for (int i = 0; i < vector_sz; ++i) {
        if (plaintextDecAdd12->GetPackedValue()[i] !=
            (vectorOfPt1[i] + vectorOfPt2[i])) {
            std::cout << "Error in addition of ciphertexts #1 and #2 at index " << i << std::endl;
            std::cout << "Expected: " << (vectorOfPt1[i] + vectorOfPt2[i]) << ", got: "
                      << plaintextDecAdd12->GetPackedValue()[i] << std::endl;
            return 1;
        }
    }    

    //Decrypt subtraction result and check correctness

    Plaintext plaintextDecSub12;

    TOC(t);
    cryptoContext->Decrypt(keyPair.secretKey, ciphertextSub12, &plaintextDecSub12);
    processingTime = TOC(t);
    std::cout << "Decryption time: " << processingTime << "ms" << std::endl;

    plaintextDecSub12->SetLength(plaintext1->GetLength());//TODO: What to set?

    std::cout << "\nResult of homomorphic subtraction of ciphertexts #1 and #2: \n";
    //TODO: I found evan multiplying 16384-14bit numbers, are done in 5ms..!!

    std::cout << plaintextDecSub12 << std::endl;

    for (int i = 0; i < vector_sz; ++i) {
        if (plaintextDecSub12->GetPackedValue()[i] !=
            (vectorOfPt1[i] - vectorOfPt2[i])) {
            std::cout << "Error in subtraction of ciphertexts #1 and #2 at index " << i << std::endl;
            std::cout << "Expected: " << (vectorOfPt1[i] - vectorOfPt2[i]) << ", got: "
                      << plaintextDecSub12->GetPackedValue()[i] << std::endl;
            return 1;
        }
    }    

    return 0;
}


int main(int argc, char* argv[]) {
    //OpenFHEBGVrns_example();
    OpenFHEBGVrns_select(true);
}