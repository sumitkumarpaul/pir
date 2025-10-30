cmake_minimum_required(VERSION 3.16)
set(CMAKE_CXX_COMPILER "g++")
project(pir_experiments LANGUAGES C CXX ASM)
add_subdirectory(Kuku)

set(CMAKE_CXX_STANDARD 17)
enable_language(ASM OPTIONAL)

# Compiler flags from Makefile
set(COMMON_FLAGS "-fomit-frame-pointer -fno-strict-aliasing -O3 -s -DNDEBUG -fopenmp")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${COMMON_FLAGS}")

include(CheckCXXCompilerFlag)
check_cxx_compiler_flag("-msse4.2" COMPILER_SUPPORTS_MSSE)
check_cxx_compiler_flag("-maes" COMPILER_SUPPORTS_AESNI)
check_cxx_compiler_flag("-mavx" COMPILER_SUPPORTS_AVX)
check_cxx_compiler_flag("-mavx2" COMPILER_SUPPORTS_AVX2)
check_cxx_compiler_flag("-mavx512f" COMPILER_SUPPORTS_AVX512F)

# TODO: Check details about this particular flag
if(COMPILER_SUPPORTS_MSSE)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -msse4.2")
endif()

if(COMPILER_SUPPORTS_AESNI)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -maes")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DAESNI")   
endif()
if(COMPILER_SUPPORTS_AVX)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mavx")
endif()
if(COMPILER_SUPPORTS_AVX2)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mavx2")
endif()
if(COMPILER_SUPPORTS_AVX512F)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mavx512f")
endif()

# Find OpenFHE (require at least version 1.3.0)
find_package(OpenFHE 1.3.0 REQUIRED CONFIG)
if(OpenFHE_FOUND)
    message(STATUS "OpenFHE Version: ${BASE_OPENFHE_VERSION}")
else()
    message(FATAL_ERROR "PACKAGE OpenFHE NOT FOUND")
endif()

# Manually search the GMP version
# Find GMP include and libraries
find_path(GMP_INCLUDE_DIR gmp.h)
find_library(GMP_LIBRARY NAMES gmp)
find_library(GMPXX_LIBRARY NAMES gmpxx)

# Extract GMP version from gmp.h
if(GMP_INCLUDE_DIR)
    file(READ "${GMP_INCLUDE_DIR}/gmp.h" GMP_H_CONTENTS)
    string(REGEX MATCH "#define[ \t]+__GNU_MP_VERSION[ \t]+([0-9]+)" _ ${GMP_H_CONTENTS})
    set(GMP_MAJOR_VERSION ${CMAKE_MATCH_1})
    string(REGEX MATCH "#define[ \t]+__GNU_MP_VERSION_MINOR[ \t]+([0-9]+)" _ ${GMP_H_CONTENTS})
    set(GMP_MINOR_VERSION ${CMAKE_MATCH_1})
    string(REGEX MATCH "#define[ \t]+__GNU_MP_VERSION_PATCHLEVEL[ \t]+([0-9]+)" _ ${GMP_H_CONTENTS})
    set(GMP_PATCH_VERSION ${CMAKE_MATCH_1})
    set(GMP_VERSION "${GMP_MAJOR_VERSION}.${GMP_MINOR_VERSION}.${GMP_PATCH_VERSION}")
endif()

# Print found version
message(STATUS "Found GMP version: ${GMP_VERSION}")

# Check version is greater than 6.3
if(GMP_VERSION VERSION_LESS "6.3")
    message(FATAL_ERROR "GMP version 6.3 or higher is required, but found ${GMP_VERSION}")
endif()

set(GMP_LIBRARIES ${GMP_LIBRARY})
set(GMPXX_LIBRARIES ${GMPXX_LIBRARY})

# Find OpenMP
find_package(OpenMP REQUIRED)

# Find Threads (for pthread)
find_package(Threads REQUIRED)

# Find OpenSSL (for ssl and crypto)
find_package(OpenSSL 1.1 REQUIRED)

# System libraries (dl is always present on Linux)
list(APPEND EXTRA_LIBS dl)

# Include directories
include_directories(
    ${OpenFHE_INCLUDE}
    ${GMP_INCLUDE_DIR}
    /usr/local/include/openfhe/core
    /usr/local/include/openfhe/binfhe
    /usr/local/include/openfhe
    /usr/local/include/openfhe/pke
)

# Link directories
link_directories(
    ${OpenFHE_LIBDIR}
    /usr/local/lib
)

# Add all your executables here
add_executable(sumit-experiment sumit-experiment.cpp)
#add_executable(pir-test fss/aesni-x86_64.S fss/aesni_wrapper.c pir-test.cpp fss/fss-client.cpp fss/fss-server.cpp fss/fss-common.cpp pir_common.cpp)
add_executable(server_alpha fss/aesni-x86_64.S fss/aesni_wrapper.c server_alpha.cpp fss/fss-client.cpp fss/fss-server.cpp fss/fss-common.cpp pir_common.cpp)
add_executable(server_beta fss/aesni-x86_64.S fss/aesni_wrapper.c server_beta.cpp fss/fss-client.cpp fss/fss-server.cpp fss/fss-common.cpp pir_common.cpp)
add_executable(server_gamma fss/aesni-x86_64.S fss/aesni_wrapper.c server_gamma.cpp fss/fss-client.cpp fss/fss-server.cpp fss/fss-common.cpp pir_common.cpp)
add_executable(pir_client pir_client.cpp pir_common.cpp)

# Link libraries for sumit-experiment
target_link_libraries(sumit-experiment
    ${GMPXX_LIBRARIES}
    ${GMP_LIBRARIES}
    OPENFHEpke
    OPENFHEcore
    OpenMP::OpenMP_CXX
    Threads::Threads
    OpenSSL::SSL
    OpenSSL::Crypto
    crypto
    ${EXTRA_LIBS}
)

# Link libraries for pir-test
#target_link_libraries(pir-test
#    ${GMPXX_LIBRARIES}
#    ${GMP_LIBRARIES}
#    OpenMP::OpenMP_CXX
#    Threads::Threads
#    OpenSSL::SSL
#    OpenSSL::Crypto
#    crypto
#    ${EXTRA_LIBS}
#)

# Link libraries for server_alpha
target_link_libraries(server_alpha
    ${GMPXX_LIBRARIES}
    ${GMP_LIBRARIES}
    OPENFHEpke
    OPENFHEcore
    OpenMP::OpenMP_CXX
    Threads::Threads
    OpenSSL::SSL
    OpenSSL::Crypto
    crypto
    Kuku::kuku
    ${EXTRA_LIBS}
)

# Link libraries for server_beta
target_link_libraries(server_beta
    OPENFHEpke
    OPENFHEcore
    ${GMPXX_LIBRARIES}
    ${GMP_LIBRARIES}
    OpenMP::OpenMP_CXX
    Threads::Threads
    OpenSSL::SSL
    OpenSSL::Crypto
    crypto
    Kuku::kuku
    ${EXTRA_LIBS}
)

# Link libraries for server_gamma
target_link_libraries(server_gamma
    ${GMPXX_LIBRARIES}
    ${GMP_LIBRARIES}
    OPENFHEpke
    OPENFHEcore
    OpenMP::OpenMP_CXX
    Threads::Threads
    OpenSSL::SSL
    OpenSSL::Crypto
    crypto
    Kuku::kuku
    ${EXTRA_LIBS}
)

# Link libraries for pir_client
target_link_libraries(pir_client
    ${GMPXX_LIBRARIES}
    ${GMP_LIBRARIES}
    OPENFHEpke
    OPENFHEcore
    OpenSSL::SSL
    OpenSSL::Crypto
    crypto
    Kuku::kuku
    ${EXTRA_LIBS}
)
# Optionally, add more executables as needed:
# add_executable(other-exe other-source.cpp)
# target_link_libraries(other-exe ...)

# Optionally, install targets
# install(TARGETS sumit-experiment pir-test DESTINATION bin)