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
File: Kuku-2.1/kuku/kuku.h
```
:
#include <set>
#include <stdexcept>
#include <vector>
// Add from here 
#include <ostream>
#include <istream>
// To here 
:
namespace kuku
        inline double fill_rate() const noexcept
        {
            return static_cast<double>(inserted_items_) /
                   (static_cast<double>(table_size()) + static_cast<double>(stash_size_));
        }
        :
        // Add from here 
        inline void serialize(std::ostream& out) const {
            // Write constructor parameters
            out.write(reinterpret_cast<const char*>(&table_size_), sizeof(table_size_));
            out.write(reinterpret_cast<const char*>(&stash_size_), sizeof(stash_size_));
            std::uint32_t loc_func_count_val = loc_func_count();
            out.write(reinterpret_cast<const char*>(&loc_func_count_val), sizeof(loc_func_count_val));
            out.write(reinterpret_cast<const char*>(loc_func_seed_.data()), sizeof(loc_func_seed_));
            out.write(reinterpret_cast<const char*>(&max_probe_), sizeof(max_probe_));
            out.write(reinterpret_cast<const char*>(empty_item_.data()), sizeof(empty_item_));

            // Write table contents
            std::uint64_t table_vec_size = table_.size();
            out.write(reinterpret_cast<const char*>(&table_vec_size), sizeof(table_vec_size));
            for (const auto& item : table_) {
                out.write(reinterpret_cast<const char*>(item.data()), sizeof(item));
            }

            // Write stash contents
            std::uint64_t stash_vec_size = stash_.size();
            out.write(reinterpret_cast<const char*>(&stash_vec_size), sizeof(stash_vec_size));
            for (const auto& item : stash_) {
                out.write(reinterpret_cast<const char*>(item.data()), sizeof(item));
            }

            // Write leftover_item_
            out.write(reinterpret_cast<const char*>(leftover_item_.data()), sizeof(leftover_item_));

            // Write inserted_items_
            out.write(reinterpret_cast<const char*>(&inserted_items_), sizeof(inserted_items_));

            // Write total_probe_count_
            out.write(reinterpret_cast<const char*>(&total_probe_count_), sizeof(total_probe_count_));
    }        

    static inline std::unique_ptr<KukuTable> deserialize(std::istream& in) {
        table_size_type table_size;
        table_size_type stash_size;
        std::uint32_t loc_func_count_val;
        item_type loc_func_seed;
        std::uint64_t max_probe;
        item_type empty_item;

        // Read constructor parameters
        in.read(reinterpret_cast<char*>(&table_size), sizeof(table_size));
        in.read(reinterpret_cast<char*>(&stash_size), sizeof(stash_size));
        in.read(reinterpret_cast<char*>(&loc_func_count_val), sizeof(loc_func_count_val));
        in.read(reinterpret_cast<char*>(loc_func_seed.data()), sizeof(loc_func_seed));
        in.read(reinterpret_cast<char*>(&max_probe), sizeof(max_probe));
        in.read(reinterpret_cast<char*>(empty_item.data()), sizeof(empty_item));

        // Construct the table
        auto table = std::make_unique<KukuTable>(table_size, stash_size, loc_func_count_val, loc_func_seed, max_probe, empty_item);

        // Read table contents
        std::uint64_t table_vec_size;
        in.read(reinterpret_cast<char*>(&table_vec_size), sizeof(table_vec_size));
        table->table_.resize(table_vec_size);
        for (auto& item : table->table_) {
            in.read(reinterpret_cast<char*>(item.data()), sizeof(item));
        }

        // Read stash contents
        std::uint64_t stash_vec_size;
        in.read(reinterpret_cast<char*>(&stash_vec_size), sizeof(stash_vec_size));
        table->stash_.resize(stash_vec_size);
        for (auto& item : table->stash_) {
            in.read(reinterpret_cast<char*>(item.data()), sizeof(item));
        }

        // Read leftover_item_
        in.read(reinterpret_cast<char*>(table->leftover_item_.data()), sizeof(table->leftover_item_));

        // Read inserted_items_
        in.read(reinterpret_cast<char*>(&table->inserted_items_), sizeof(table->inserted_items_));

        // Read total_probe_count_
        in.read(reinterpret_cast<char*>(&table->total_probe_count_), sizeof(table->total_probe_count_));

        // Regenerate location functions
        table->generate_loc_funcs(loc_func_count_val, loc_func_seed);

        return table;
    }

    // To here
    private:
        KukuTable(const KukuTable &copy) = delete;
    :
```
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
