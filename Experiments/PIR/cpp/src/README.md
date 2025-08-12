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
