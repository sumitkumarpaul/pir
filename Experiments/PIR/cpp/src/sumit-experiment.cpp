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
// 1. Redistributions of source code must retain the above copyright noTIC1e, this
//    list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright noTIC1e,
//    this list of conditions and the following disclaimer in the documentation
//    and/or other materials provided with the distribution.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTIC1ULAR PURPOSE ARE
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

#define TIC1(t)    t = std::chrono::high_resolution_clock::now()
#define TOC1(t)    std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - t).count()

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
    TIC1(t);

    keyPair = cryptoContext->KeyGen();

    processingTime = TOC1(t);
    std::cout << "Key generation time: " << processingTime << "us" << std::endl;

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

    TIC1(t);

    ciphertexts.push_back(cryptoContext->Encrypt(keyPair.publicKey, plaintext1));
    ciphertexts.push_back(cryptoContext->Encrypt(keyPair.publicKey, plaintext2));

    processingTime = TOC1(t);

    //std::cout << "Completed\n";

    std::cout << "\nAverage encryption time: " << processingTime / 2 << "us" << std::endl;


    ////////////////////////////////////////////////////////////
    // Homomorphic multiplication of two ciphertexts w/o any relinearization
    ////////////////////////////////////////////////////////////

    //std::cout << "\nRunning a multiplication of two ciphertexts w/o relinearization...";
    
    TIC1(t);

    auto ciphertextMult12 = cryptoContext->EvalMultNoRelin(ciphertexts[0], ciphertexts[1]);
    cryptoContext->ModReduceInPlace(ciphertextMult12);

    processingTime = TOC1(t);
    std::cout << "Time for multiplying two ciphertexts w/o relinearization: " << processingTime << "us" << std::endl;

    //std::cout << "Completed\n";  

    Plaintext plaintextDecMult12;

    TOC1(t);
    cryptoContext->Decrypt(keyPair.secretKey, ciphertextMult12, &plaintextDecMult12);
    processingTime = TOC1(t);
    std::cout << "Decryption time: " << processingTime << "us" << std::endl;

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
    parameters.SetPlaintextModulus(65537);//TODO, 65537, 536903681 these values must have special properties.
    //At this moment, using the value mentioned in the original example.
    parameters.SetMaxRelinSkDeg(1);//What does this value mean?
    parameters.SetScalingTechnique(FIXEDAUTO);//Only this is not giving any exception and giving good result

    CryptoContext<DCRTPoly> cryptoContext = GenCryptoContext(parameters);
    // enable features that you wish to use
    cryptoContext->Enable(PKE);
    //cryptoContext->Enable(KEYSWITCH);
    cryptoContext->Enable(LEVELEDSHE);
    //cryptoContext->Enable(ADVANCEDSHE);

    // Initialize Public Key Containers
    KeyPair<DCRTPoly> keyPair;

    // Perform Key Generation Operation
    TIC1(t);

    keyPair = cryptoContext->KeyGen();

    processingTime = TOC1(t);
    std::cout << "Key generation time: " << processingTime << " us" << std::endl;

    if (!keyPair.good()) {
        std::cout << "Key generation failed!" << std::endl;
        exit(1);
    }

    ////////////////////////////////////////////////////////////
    // Encode source data
    ////////////////////////////////////////////////////////////
    int vector_sz = 5000;//314; // Size of the vector to be generated. TODO: Maximum usable is: 16384. Each FHE block is 14-bits. So 51200-bit size data block will take 3657 vector elements. Increased that to 5000.

    std::vector<int64_t> vectorOfPt1;
    // Use a for loop to add elements to the vector
    for (int i = 0; i < vector_sz; ++i) {
        //int num = rand() % range + min;
        int num = 1;
        vectorOfPt1.push_back(num);
    }

    Plaintext plaintext1 = cryptoContext->MakePackedPlaintext(vectorOfPt1);

    std::vector<int64_t> vectorOfPt2;
    // Use a for loop to add elements to the vector
    for (int i = 0; i < vector_sz; ++i) {
        //int num = rand() % range + min;
        int num = 255;//Previous was 32767
        vectorOfPt2.push_back(num);
    }
    Plaintext plaintext2 = cryptoContext->MakePackedPlaintext(vectorOfPt2);

    std::vector<int64_t> vectorOfSelect;
    // Use a for loop to add elements to the vector
    for (int i = 0; i < vector_sz; ++i) {
        vectorOfSelect.push_back(select_bit);
    }
    Plaintext plaintextSelect = cryptoContext->MakePackedPlaintext(vectorOfSelect);
    std::vector<int64_t> vectorOfOnes;
    // Use a for loop to add elements to the vector
    for (int i = 0; i < vector_sz; ++i) {
        vectorOfOnes.push_back(1);
    }
    Plaintext plaintextOnes = cryptoContext->MakePackedPlaintext(vectorOfOnes);

    ////////////////////////////////////////////////////////////
    // Encryption
    ////////////////////////////////////////////////////////////

    std::vector<Ciphertext<DCRTPoly>> ciphertexts;

    TIC1(t);

    auto ct1 = cryptoContext->Encrypt(keyPair.publicKey, plaintext1);
    auto ct2 = cryptoContext->Encrypt(keyPair.publicKey, plaintext2);
    auto ctSel = cryptoContext->Encrypt(keyPair.publicKey, plaintextSelect);
    auto ctOnes = cryptoContext->Encrypt(keyPair.publicKey, plaintextOnes);

    processingTime = TOC1(t);

    std::cout << "Encryption time: " << (processingTime/4) << " us" << std::endl;

    ////////////////////////////////////////////////////////////
    // Homomorphic selection between two ciphertexts w/o any relinearization
    // Select(a,b, select_bit) = select_bit ? a : b
    //                         = ((1 - select_bit) * a) + (select_bit * b)
    ////////////////////////////////////////////////////////////
    
    TIC1(t);
    
    auto ciphertextSelNot = cryptoContext->EvalSub(ctOnes, ctSel);
    //cryptoContext->ModReduceInPlace(ciphertextSelNot);

    //processingTime = TOC1(t);
    //std::cout << "Time for first sub between two ciphertexts: " << processingTime << "us" << std::endl;

    //TIC1(t);

    auto ciphertextSelNota = cryptoContext->EvalMultNoRelin(ciphertextSelNot, ct1);
    //cryptoContext->ModReduceInPlace(ciphertextSelNota);

    //processingTime = TOC1(t);
    //std::cout << "Time for first mult. between two ciphertexts w/o relin: " << processingTime << "us" << std::endl;

    //TIC1(t);

    auto ciphertextSelb = cryptoContext->EvalMultNoRelin(ctSel, ct2);
    //cryptoContext->ModReduceInPlace(ciphertextSelb);

    //processingTime = TOC1(t);
    //std::cout << "Time for sencond mult. between two ciphertexts w/o relin: " << processingTime << "us" << std::endl;

    //TIC1(t);

    auto ciphertextSel = cryptoContext->EvalAdd(ciphertextSelNota, ciphertextSelb);
    //cryptoContext->ModReduceInPlace(ciphertextSel);

    //processingTime = TOC1(t);
    //std::cout << "Time for addition: " << processingTime << "us" << std::endl;

    processingTime = TOC1(t);
    std::cout << "Time for homomorphic selection: " << processingTime << " us" << std::endl;

    Plaintext plaintextDecSel;

    TIC1(t);
    cryptoContext->Decrypt(keyPair.secretKey, ciphertextSel, &plaintextDecSel);
    processingTime = TOC1(t);
    std::cout << "Decryption time: " << processingTime << "us" << std::endl;

    plaintextDecSel->SetLength(plaintext1->GetLength());//TODO: What to set?

    for (int i = 0; i < vector_sz; ++i) {
        if (plaintextDecSel->GetPackedValue()[i] !=
            ((select_bit == 1) ? vectorOfPt2[i]: vectorOfPt1[i])) {
            std::cout << "Error in selection of ciphertexts #1 and #2 at index " << i << std::endl;
            std::cout << "Expected: " << ((select_bit == 1) ? vectorOfPt2[i]: vectorOfPt1[i]) << ", got: "
                      << plaintextDecSel->GetPackedValue()[i] << std::endl;
            return 1;
        }
    }    

    return 0;
}


int main(int argc, char* argv[]) {
    //OpenFHEBGVrns_example();
    #if 0
    OpenFHEBGVrns_select(true);
    OpenFHEBGVrns_select(false);
    OpenFHEBGVrns_select(true);
    OpenFHEBGVrns_select(false);
    OpenFHEBGVrns_select(false);
    OpenFHEBGVrns_select(true);
    #endif
    OpenFHEBGVrns_select(true);
}