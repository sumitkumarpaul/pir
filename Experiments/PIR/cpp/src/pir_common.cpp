#include <gmp.h>
#include <gmpxx.h>
#include <ctime>
#include <iomanip>
#include "pir_common.h"


// Global ElGamal parameters
mpz_class p, q, r, g;

// El-Gamal encryption keys
mpz_class pk_E;
mpz_class sk_E;

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

mpz_class ElGamal_randomGroupElement() {
    gmp_randclass rng(gmp_randinit_default);
    // good-ish seed source: std::random_device (combine two samples)
    std::random_device rd;
    unsigned long seed = (static_cast<unsigned long>(rd()) << 1) ^ rd();
    rng.seed(seed); // seed() seeds the gmp_randclass

    mpz_class localr = rng.get_z_range(q);
    mpz_class result;
    mpz_powm(result.get_mpz_t(), g.get_mpz_t(), localr.get_mpz_t(), p.get_mpz_t());

    return result;
}

std::pair<mpz_class, mpz_class> ElGamal_keyGen(const mpz_class& p, const mpz_class& q, const mpz_class& g) {
    gmp_randclass rng(gmp_randinit_default);
    std::random_device rd;
    unsigned long seed = (static_cast<unsigned long>(rd()) << 1) ^ rd();
    rng.seed(seed); // seed() seeds the gmp_randclass

    mpz_class x = rng.get_z_range(q);
    mpz_class y;
    mpz_powm(y.get_mpz_t(), g.get_mpz_t(), x.get_mpz_t(), p.get_mpz_t());

    return std::make_pair(y, x);
}

std::pair<mpz_class, mpz_class> ElGamal_encrypt(const mpz_class& message, const mpz_class& publicKey) {
    gmp_randclass rng(gmp_randinit_default);
    std::random_device rd;
    unsigned long seed = (static_cast<unsigned long>(rd()) << 1) ^ rd();
    rng.seed(seed); // seed() seeds the gmp_randclass

    mpz_class k = rng.get_z_range(q);
    mpz_class c1, c2;
    mpz_powm(c1.get_mpz_t(), g.get_mpz_t(), k.get_mpz_t(), p.get_mpz_t());
    mpz_class temp;
    mpz_powm(temp.get_mpz_t(), publicKey.get_mpz_t(), k.get_mpz_t(), p.get_mpz_t());
    c2 = (message * temp) % p;
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

std::pair<mpz_class, mpz_class> ElGamal_mult_ct(const std::pair<mpz_class, mpz_class>& ciphertext1, const std::pair<mpz_class, mpz_class>& ciphertext2) {
    mpz_class cm1 = (ciphertext1.first * ciphertext2.first) % p;
    mpz_class cm2 = (ciphertext1.second * ciphertext2.second) % p;
    return std::make_pair(cm1, cm2);
}

std::pair<mpz_class, mpz_class> ElGamal_exp_ct(const std::pair<mpz_class, mpz_class>& ciphertext, const mpz_class& exp, const mpz_class& publicKey) {
    mpz_class c1, c2;
    mpz_powm(c1.get_mpz_t(), ciphertext.first.get_mpz_t(), exp.get_mpz_t(), p.get_mpz_t());
    mpz_powm(c2.get_mpz_t(), ciphertext.second.get_mpz_t(), exp.get_mpz_t(), p.get_mpz_t());
    auto [cI1, cI2] = ElGamal_encrypt(mpz_class(1), publicKey);
    return ElGamal_mult_ct({c1, c2}, {cI1, cI2});
}

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

