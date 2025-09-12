# Install OpenSSL library
```
sudo apt-get install libssl-dev
```

# Install openmp (Required for parallelization)
```
sudo apt install libomp-dev
```

# Install m4 (Required for text processing)
```
sudo apt-get install m4
```

# Build and install OpenFHE library (for FHE encryption) in the default path

Details can be found [here] (https://openfhe-development.readthedocs.io/en/latest/sphinx_rsts/intro/installation/linux.html).

Specifically, first install the required dependencies.
```
sudo apt-get install build-essential
sudo apt-get install cmake
```

Download and compile the latest source code, and then install
```
git clone https://github.com/openfheorg/openfhe-development.git
cd openfhe-development/
mkdir build
cd build
cmake ..
make
sudo make install
```

Make sure the newly installed library is accessible

```
echo "/usr/local/lib" | sudo tee /etc/ld.so.conf.d/openfhe.conf
sudo ldconfig
```


# TODO use OpenFHE for best version with INTELHEXL?

TODO
https://github.com/openfheorg/openfhe-development/blob/main/docs/static_docs/Best_Performance.md

# Install GMP (to perform arithmatic operations with large numbers) with AVX512 support

```
wget https://ftp.gnu.org/gnu/gmp/gmp-6.3.0.tar.xz
tar -xvf gmp-6.3.0.tar.xz
cd gmp-6.3.0
./configure CFLAGS="-mavx512f -O3" CXXFLAGS="-mavx512f -O3"  --enable-cxx
make
sudo make install
```

# Modify the Kuku library (Required for Cuckoo hash functions)
File: Kuku-2.1/kuku/common.h
```
    //using table_size_type = location_type;
    using table_size_type = std::uint64_t;//To support larger table_size
    :
    //constexpr table_size_type max_table_size = table_size_type(1) << 30;
    constexpr table_size_type max_table_size = table_size_type(1) << 32;//To support larger table_size to support 100 GB database
```


# Install Kuku library (Required for Cuckoo hash functions)
```
cmake -S . -B build
cmake --build build
sudo cmake --install build
```
Details can be found [here](https://github.com/microsoft/Kuku.git).

# Download and compile the pir source code
```
git clone https://github.com/sumitkumarpaul/pir.git
cd pir/Experiments/PIR/cpp/src/
mkdir build
cmake ..
make
```
