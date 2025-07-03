use tfhe::prelude::*;
use tfhe::{generate_keys, set_server_key, ConfigBuilder, FheUint32, FheUint64, FheUint8};
use std::time::Instant;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // Basic configuration to use homomorphic integers
    let config = ConfigBuilder::default().build();

    // Key generation
    let (client_key, server_keys) = generate_keys(config);

    let clear_a = 1344u64;
    let clear_b = 5u64;
    let clear_c = 7u8;

    // Encrypting the input data using the (private) client_key
    // FheUint32: Encrypted equivalent to u32
    let start = Instant::now();
    let mut encrypted_a = FheUint32::try_encrypt(clear_a, &client_key)?;
    println!("Elapsed time for encrypting clear_a: {:?}", start.elapsed());

    let start = Instant::now();
    let mut encrypted = FheUint64::try_encrypt(clear_a, &client_key)?;
    println!("Elapsed time for encrypting clear_a: {:?}", start.elapsed());

    let start = Instant::now();
    let encrypted_b = FheUint32::try_encrypt(clear_b, &client_key)?;
    println!("Elapsed time for encrypting clear_b: {:?}", start.elapsed());

    // FheUint8: Encrypted equivalent to u8
    let start = Instant::now();
    let encrypted_c = FheUint8::try_encrypt(clear_c, &client_key)?;
    println!("Elapsed time for encrypting clear_c: {:?}", start.elapsed());

    // On the server side:
    set_server_key(server_keys);

    // Select the greater value between a and b
    let encrypted_comp = &encrypted_a.gt(&encrypted_b);

    let start = Instant::now();
    let encrypted_res = &encrypted_comp.select(&encrypted_a, &encrypted_b);
    println!("Elapsed time for selecting the result: {:?}", start.elapsed());
    
    //let clear_res: FheUint64 = encrypted_res.decrypt(&client_key);    

    // Decrypting on the client side:
    let start = Instant::now();
    let clear_res: u8 = encrypted_res.decrypt(&client_key);
    println!("Elapsed time for decrypting result: {:?}", start.elapsed());
    assert_eq!(clear_res, 1_u8);

    Ok(())
}