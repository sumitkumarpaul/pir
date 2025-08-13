#include "pir_common.h"

int InitAcceptingSocket(int port, int* p_server_fd, int* p_new_socket) {
    struct sockaddr_in address;
    int opt = 1;
    int addrlen = sizeof(address);
    char buffer[65536] = {0};
    *p_new_socket = -1;
    *p_server_fd = -1;
    int ret = -1;

    *p_server_fd = socket(AF_INET, SOCK_STREAM, 0);
    if (*p_server_fd < 0) {
        std::cerr << "Server: Socket creation failed" << std::endl;
        return ret;
    }
    if (setsockopt(*p_server_fd, SOL_SOCKET, SO_REUSEADDR | SO_REUSEPORT, &opt, sizeof(opt)) < 0) {
        std::cerr << "Server: setsockopt failed" << std::endl;
        close(*p_server_fd);
        *p_server_fd = -1;
        return ret;
    }
    address.sin_family = AF_INET;
    address.sin_addr.s_addr = INADDR_ANY;
    address.sin_port = htons(port);
    if (bind(*p_server_fd, (struct sockaddr *)&address, sizeof(address)) < 0) {
        std::cerr << "Server: Bind failed" << std::endl;
        close(*p_server_fd);
        *p_server_fd = -1;
        return ret;
    }
    if (listen(*p_server_fd, 1) < 0) {
        std::cerr << "Server: Listen failed" << std::endl;
        close(*p_server_fd);
        *p_server_fd = -1;
        return ret;
    }
    *p_new_socket = accept(*p_server_fd, (struct sockaddr *)&address, (socklen_t*)&addrlen);
    if (*p_new_socket < 0) {
        std::cerr << "Server: Accept failed" << std::endl;
        close(*p_server_fd);
        *p_server_fd = -1;
        return ret;
    }

    return 0;
}

void FinishAcceptingSocket(int server_fd, int new_socket) {
    close(new_socket);
    close(server_fd);
}

void InitConnectingSocket(const std::string& server_ip, int port, int* p_sock) {
    // Connect to server
    struct sockaddr_in serv_addr;
    *p_sock = socket(AF_INET, SOCK_STREAM, 0);
    serv_addr.sin_family = AF_INET;
    serv_addr.sin_port = htons(port);
    inet_pton(AF_INET, server_ip.c_str(), &serv_addr.sin_addr);
    while (connect(*p_sock, (struct sockaddr *)&serv_addr, sizeof(serv_addr)) < 0) {
        std::this_thread::sleep_for(std::chrono::milliseconds(100)); // Wait for server
    }

    return;
}