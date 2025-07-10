use tfhe::prelude::*;
use tfhe::*;
use tfhe::shortint::prelude::*;
use std::time::Instant;
use tfhe::integer::gen_keys_radix;
use rand::thread_rng;
use tfhe::integer::{U256, RadixClientKey, ServerKey};
use tfhe::integer::parameters::PARAM_MESSAGE_1_CARRY_1_KS_PBS_32_BITS;
use tfhe::integer::prelude::*;
use tfhe::shortint::parameters::PARAM_MESSAGE_2_CARRY_2_KS_PBS_GAUSSIAN_2M128;
use tfhe::shortint::parameters::PARAM_MESSAGE_2_CARRY_2_KS_PBS_TUNIFORM_2M128;

fn radix_if_then_else_pbs() {
    let size = 32; // 32 * 2 = 64 bits of message space
    let (cks, sks) = gen_keys_radix(PARAM_MESSAGE_1_CARRY_1_KS_PBS_32_BITS, size);

    let a = 128u64;
    let b = 55u64;

    let t_start = Instant::now();
    let ct_a = cks.encrypt(a);
    let ct_b = cks.encrypt(b);
    let t_end = Instant::now();
    println!("Encryption time: {:.3} ms", (t_end - t_start).as_secs_f64() * 1000.0);

    let condition = sks.scalar_ge_parallelized(&ct_a, 66);
    reset_pbs_count();
    let t_start = Instant::now();
    let ct_res = sks.if_then_else_parallelized(&condition, &ct_a, &ct_b);//Currently, there is no difference between `if_then_else_parallelized` and `if_then_else` for the radix PBS.
    let t_end = Instant::now();
    println!("If-then-else evaluation time: {:.3} ms", (t_end - t_start).as_secs_f64() * 1000.0);
    let radix_unchecked_if_then_else_pbs_count = get_pbs_count();

    println!("radix_unchecked_if_then_else_pbs_count: {radix_unchecked_if_then_else_pbs_count}");

    // Decrypt:
    let t_start = Instant::now();
    let dec: u64 = cks.decrypt(&ct_res);
    let t_end = Instant::now();
    println!("Decryption time: {:.3} ms", (t_end - t_start).as_secs_f64() * 1000.0);
    assert_eq!(if a >= 66 { a } else { b }, dec);
    assert_ne!(ct_a, ct_res);
    assert_ne!(ct_b, ct_res);
}

fn if_then_else_1() {
    // Choose a parameter set (adjust as needed)
    let params = tfhe::shortint::parameters::PARAM_MESSAGE_2_CARRY_2_KS_PBS_TUNIFORM_2M128;
    let num_blocks = 32; // Sumit: Maybe this one is not simply bitsize/2

    // Generate keys
    let cks = RadixClientKey::new(params, num_blocks);
    let sks = ServerKey::new_radix_server_key(&cks);

    // Encrypt two random U256 values
    let mut rng = thread_rng();
    let clear_0 = U256::from(rand::random::<u128>());
    let clear_1 = U256::from(rand::random::<u128>());
    let ct_0 = cks.encrypt(clear_0);
    let ct_1 = cks.encrypt(clear_1);

    // Encrypt a random boolean as the condition
    let cond = sks.create_trivial_boolean_block(rand::random::<bool>());

    // Time the if_then_else_parallelized operation
    let start = Instant::now();
    let result = sks.if_then_else_parallelized(&cond, &ct_0, &ct_1);
    let elapsed = start.elapsed();

    // Decrypt and print results
    let decrypted: U256 = cks.decrypt(&result);
    
    println!("True branch: {:?}", clear_0);
    println!("False branch: {:?}", clear_1);
    println!("Result: {:?}", decrypted);
    println!("if_then_else_parallelized took {:.3} ms", elapsed.as_secs_f64() * 1000.0);
}

fn radix_add_pbs() {
    let num_block = 32;//32*2 = 64 bits of message space
    let (client_key, server_key) = gen_keys_radix(PARAM_MESSAGE_2_CARRY_2_KS_PBS, num_block);

    let msg1 = 0x1234123412341234_u64;
    let msg2 = 0x3456345634563456_u64;

    // message_modulus^vec_length
    let modulus = client_key.parameters().message_modulus().0.pow(num_block as u32);

    // We use the client key to encrypt two messages:
    let t_start = Instant::now();
    let ct_1 = client_key.encrypt(msg1);
    let ct_2 = client_key.encrypt(msg2);
    let t_end = Instant::now();
    let unchecked_add_count = get_pbs_count();
    println!("Encryption time: {:.3} ms", (t_end - t_start).as_secs_f64() * 1000.0);

    reset_pbs_count();

    let t_start = Instant::now();    
    let ct_res = server_key.unchecked_add(&ct_1, &ct_2);
    let t_end = Instant::now();
    let radix_unchecked_add_pbs_count = get_pbs_count();

    println!("radix_unchecked_add_pbs_count: {radix_unchecked_add_pbs_count}");
    println!("Addition time: {:.3} ms", (t_end - t_start).as_secs_f64() * 1000.0);

    // We use the client key to decrypt the output of the circuit:
    let t_start = Instant::now();      
    let pt_res: u64 = client_key.decrypt(&ct_res);
    let t_end = Instant::now();
    println!("Decryption time: {:.3} ms", (t_end - t_start).as_secs_f64() * 1000.0);

    // The carry buffer has been overflowed, the result is not correct
    assert_ne!(pt_res, (msg1 + msg2) % modulus);
}

fn radix_and_pbs() {
    let num_block = 32;//32*2 = 64 bits of message space
    let (client_key, server_key) = gen_keys_radix(PARAM_MESSAGE_2_CARRY_2_KS_PBS, num_block);

    let msg1 = 0x1234123412341234_u64;
    let msg2 = 0x3456345634563456_u64;

    // message_modulus^vec_length
    let modulus = client_key.parameters().message_modulus().0.pow(num_block as u32);

    // We use the client key to encrypt two messages:
    let t_start = Instant::now();
    let ct_1 = client_key.encrypt(msg1);
    let ct_2 = client_key.encrypt(msg2);
    let t_end = Instant::now();
    let unchecked_add_count = get_pbs_count();
    println!("Encryption time: {:.3} ms", (t_end - t_start).as_secs_f64() * 1000.0);

    reset_pbs_count();

    let t_start = Instant::now();    
    let ct_res = server_key.unchecked_bitand_parallelized(&ct_1, &ct_2);
    let t_end = Instant::now();
    let radix_unchecked_bitand_pbs_count = get_pbs_count();

    println!("radix_unchecked_bitand_pbs_count: {radix_unchecked_bitand_pbs_count}");
    println!("Bitwise AND time: {:.3} ms", (t_end - t_start).as_secs_f64() * 1000.0);

    // We use the client key to decrypt the output of the circuit:
    let t_start = Instant::now();      
    let pt_res: u64 = client_key.decrypt(&ct_res);
    let t_end = Instant::now();
    println!("Decryption time: {:.3} ms", (t_end - t_start).as_secs_f64() * 1000.0);

    // The carry buffer has been overflowed, the result is not correct
    assert_ne!(pt_res, (msg1 + msg2) % modulus);
}

fn FheUint32_pbs(){
    let config = ConfigBuilder::default().build();

    let (cks, sks) = generate_keys(config);

    let a = FheUint32::encrypt(42u32, &cks);
    let b = FheUint32::encrypt(69u32, &cks);

    set_server_key(sks);

    reset_pbs_count();
    
    let c = &a * &b;
    let mul_32_count = get_pbs_count();

    reset_pbs_count();
    let d = &a & &b;
    let and_32_count = get_pbs_count();

    println!("mul_32_count: {mul_32_count}");
    println!("and_32_count: {and_32_count}");

    let c_dec: u32 = c.decrypt(&cks);
    let d_dec: u32 = d.decrypt(&cks);

    assert_eq!(42 * 69, c_dec);
    assert_eq!(42 & 69, d_dec);
}

pub fn main() {
    //FheUint32_pbs();
    //unchecked_add_pbs();
    //radix_add_pbs();
    radix_if_then_else_pbs();
    //radix_and_pbs();
    //if_then_else_1();
}

