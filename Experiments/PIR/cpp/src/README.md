# Build and install OpenFHE library (for FHE encryption) in the default path

Details can be found [here] (https://openfhe-development.readthedocs.io/en/latest/sphinx_rsts/intro/installation/installation.html).

# TODO use OpenFHE for best version with INTELHEXL?

TODO
https://github.com/openfheorg/openfhe-development/blob/main/docs/static_docs/Best_Performance.md

# Install GMP (to perform arithmatic operations with large numbers) with AVX512 support

Follow this [link](https://gmplib.org/manual/Installing-GMP) but use the following commandline arguments.
```
./configure CFLAGS="-mavx512f -O3" CXXFLAGS="-mavx512f -O3"
sudo apt-get install m4
```

# Install openmp (Required for parallelization)
```
sudo apt install libomp-dev
```

# Install FLINT library (for number theory (modular operation, group, field etc.) operations)
Check their Readme documentation [here](https://github.com/flintlib/flint)
Enable avx512 support
```
./configure --enable-avx512
```

## As a dependency you may require to install latest MPFR library
Check documentation from [here](https://www.mpfr.org/mpfr-current/mpfr.html#Installing-MPFR)