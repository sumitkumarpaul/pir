# About PPIR:
PPIR is a multi-server PIR system to hide a client's access pattern from a public database.
PPIR uses three non-colluding servers, $\mathsf{S_{\alpha},\ S_{\beta}\ \&\ S_{\gamma}}$.
By combining the concepts from DORAM, it achieves server computation overhead of $O(\mathsf{\sqrt{N}})$, where $\mathsf{N}$ is the total number of items in the public database.

In our implementation, we utilized the [GMP 6.3.0 library](https://gmplib.org/) for large number arithmetic.
For DPF operations, we modified and repurposed the library created by [Wang et al.](https://github.com/frankw2/libfss.git).
Our protocol employs [BGV](https://dl.acm.org/doi/10.1145/2090236.2090262) as the underlying fully homomorphic encryption (FHE) scheme, leveraging the available implementation from the [OpenFHE library](https://openfhe-development.readthedocs.io/en/latest/). Additionally, for cuckoo hashing, we used the [Kuku library](https://github.com/microsoft/Kuku.git).

For our experiment, we utilized servers with 16-core CPUs and 120GB of RAM, while simulating a low-power client using a single-core machine with 4GB of RAM, although actual RAM utilization was just a few megabytes.

PPIR is a computational PIR system, which we configured for $\mathsf{128}$-bit security during our experimentation. Consequently, we have selected the tag size ($\mathsf{\Vert p\Vert}$) to be 3072 bits and the group size of El-Gamal encryption ($\mathsf{\Vert q\Vert}$) to be (256 + 2), which totals 258 bits. Additionally, we have chosen the size of the DPF-search tag ($\mathsf{\Vert r\Vert}$) to be 64 bits, resulting in a collision probability of approximately [$\mathsf{4.55\times 10^{-11}}$](https://www.perplexity.ai/search/consider-a-list-l-t-1-t-2-ldot-CLLEHd.1QGCsRi0JEZ15Kg#0) during the shelter update operation.


# Requirement:
To set up PPIR, three server computers and one client computer, all with internet connectivity, are required.
All four computers must have a fresh installation of Ubuntu 20.04 Operating System.

# How to compile PPIR:

## Prepare the system

On each computers perform the following steps to compile PPIR.
### Install necessary build tools
```
sudo apt-get install build-essential
sudo apt-get install cmake
```

### Install pre-built libraries
For cryptographic operations, the **OpenSSL** library is required. The **m4** library is required as a dependency of other required libraries. The **OpenMP** library is required to exploit CPU-level parallelism wherever possible. These pre-built libraries can be installed with the following commands.
```
sudo apt-get install libssl-dev
sudo apt install libomp-dev
sudo apt-get install m4
```


### Build and install OpenFHE library
PPIR utilizes homomorphic encryption. **OpenFHE** library is used for that purpose. You can download it from the GitHub repository and compile it for installation on your system. For detailed instructions, please refer to the documentation [here](https://openfhe-development.readthedocs.io/en/latest/sphinx_rsts/intro/installation/linux.html). However, the following commands are sufficient for our experimentation.
```
git clone https://github.com/openfheorg/openfhe-development.git
cd openfhe-development/
mkdir build
cd build
cmake ..
make
sudo make install
```

To make sure the newly installed library is accessible by other applications, please issue the following commands.

```
echo "/usr/local/lib" | sudo tee /etc/ld.so.conf.d/openfhe.conf
sudo ldconfig
```


### Install GMP (to perform arithmatic operations with large numbers) with AVX512 support
This software relies on the **GMP** library for large number arithmetic operations. If the hardware supports AVX512F, the performance of this library can be enhanced by compiling it with additional flags. The following commands can be used for that purpose.
```
wget https://ftp.gnu.org/gnu/gmp/gmp-6.3.0.tar.xz
tar -xvf gmp-6.3.0.tar.xz
cd gmp-6.3.0
./configure CFLAGS="-mavx512f -O3" CXXFLAGS="-mavx512f -O3"  --enable-cxx
make
sudo make install
```

### Modify the Kuku library (Required for Cuckoo hash functions)
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

### Install Kuku library (Required for Cuckoo hash functions)
```
cmake -S . -B build
cmake --build build
sudo cmake --install build
```
Details can be found [here](https://github.com/microsoft/Kuku.git).

## Download and PPIR source code and modify configuration
After preparing the systems, download the source code from the GitHub repository onto all four computers by issuing the following command:
```
git clone https://github.com/sumitkumarpaul/pir.git
```
Next, open the file: *pir/Experiments/PIR/cpp/src/pir_common.h* to set the IP addresses of the three configured servers.
```
:
#define SERVER_ALPHA_IP "127.0.0.1" // Update with the IP address of Server_alpha
#define SERVER_BETA_IP  "127.0.0.1" // Update with the IP address of Server_beta
#define SERVER_GAMMA_IP "127.0.0.1" // Update with the IP address of Server_gamma
:
```

## Compile the modified source code
You can use the following commands for that.
```
cd pir/Experiments/PIR/cpp/src/
mkdir build
cmake ..
make
```
Upon success, you should see four executables named: *pir_client*, *server_alpha*, *server_beta*, and *server_gamma* in the *pir/Experiments/PIR/cpp/src/build* directory.

# How to run PPIR:
Go to *pir/Experiments/PIR/cpp/src/build* directory on all four computers.

## Perform One-time Initialization
First, the one-time initialization of the system is required. To do that, please execute the following commands:
### On $\mathsf{S_{\alpha}}$

```
server_alpha one_time_init
```

### On $\mathsf{S_{\beta}}$

```
server_beta one_time_init
```
### On $\mathsf{S_{\gamma}}$

```
server_gamma one_time_init
```

## Perform Per-Epoch Initialization
Then per-epoch operations are required to be performed. To do that, please execute the following commands:
### On $\mathsf{S_{\alpha}}$

```
server_alpha per_epoch_operations
```

### On $\mathsf{S_{\beta}}$

```
server_beta per_epoch_operations
```
### On $\mathsf{S_{\gamma}}$

```
server_gamma per_epoch_operations
```

## Process PIR-request
The servers must now be entered into PIR-request processing state. To do that, please execute the following commands:
### On $\mathsf{S_{\alpha}}$

```
server_alpha process_request
```

### On $\mathsf{S_{\beta}}$

```
server_beta process_request
```
### On $\mathsf{S_{\gamma}}$

```
server_gamma process_request
```
After these steps, the servers will stay active for an entire epoch. In other words, the servers can handle up to $\mathsf{\sqrt{N}}$ PIR requests from *any* PIR client.

### On Client computer
Now the PIR client can issue the request. To do that issue the following command:

```
pir_client <index of the requested item>
```
**Note:** Upto $\mathsf{\sqrt{N}}$ PIR requests can be issued.
After that, [Per-Epoch Initialization](#Perform-Per-Epoch-Initialization) is required to be performed again.