# Configure TFHE library

## Download source code from github

```
git clone https://github.com/zama-ai/tfhe-rs.git
```
## Install rust
```
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
```

## Install latest cmake (must be > 3.24)
```
 wget -O - https://apt.kitware.com/keys/kitware-archive-latest.asc 2>/dev/null | gpg --dearmor - | sudo tee /usr/share/keyrings/kitware-archive-keyring.gpg >/dev/null
  echo 'deb [signed-by=/usr/share/keyrings/kitware-archive-keyring.gpg] https://apt.kitware.com/ubuntu/ focal-rc main' | sudo tee -a /etc/apt/sources.list.d/kitware.list >/dev/null    

sudo apt update
sudo apt install cmake
```


## Configure compilation-tool chain
Install nightly toolchain. Details can be found [here] (https://docs.zama.ai/tfhe-rs/configuration/rust_configuration)

```
rustup toolchain install nightly
rustup override set nightly
```
Verify that nightly toolchain is currently active.

```
rustup show
```

## Compile with proper command line arguments to use C-API and avx512 feature
```
RUSTFLAGS="-C target-cpu=native" cargo +nightly build --release --features=high-level-c-api,nightly-avx512 -p tfhe
```

## Locate the generated library and header file
Those can be found under: ./target/release/.
Copy them to TODO


# Install OpenSSL and GMP library
```
sudo apt-get install libssl-dev
sudo apt-get install libgmp3-dev
```

# Configure FSS library
