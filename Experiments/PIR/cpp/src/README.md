# Install OpenSSL library
```
sudo apt-get install libssl-dev
```

# Build and install OpenFHE library (for FHE encryption) in the default path

Details can be found [here] (https://openfhe-development.readthedocs.io/en/latest/sphinx_rsts/intro/installation/installation.html).

# TODO use OpenFHE for best version with INTELHEXL?

TODO
https://github.com/openfheorg/openfhe-development/blob/main/docs/static_docs/Best_Performance.md

# Install GMP (to perform arithmatic operations with large numbers) with AVX512 support

Follow this [link](https://gmplib.org/manual/Installing-GMP) but use the following commandline arguments.
```
./configure CFLAGS="-mavx512f -O3" CXXFLAGS="-mavx512f -O3"  --enable-cxx
sudo apt-get install m4
```

# Install openmp (Required for parallelization)
```
sudo apt install libomp-dev
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