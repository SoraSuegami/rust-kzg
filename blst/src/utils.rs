extern crate alloc;

use alloc::vec::Vec;

use kzg::eip_4844::hash_to_bls_field;
use kzg::{Fr, G1Mul, G2Mul};

use crate::consts::{G1_GENERATOR, G2_GENERATOR, SCALE2_ROOT_OF_UNITY};
use crate::types::fft_settings::expand_root_of_unity;
use crate::types::fr::FsFr;
use crate::types::g1::FsG1;
use crate::types::g2::FsG2;

pub fn generate_trusted_setup(n: usize, secret: [u8; 32usize]) -> (Vec<FsG1>, Vec<FsG2>) {
    let s = hash_to_bls_field(&secret);
    let mut s_pow = Fr::one();

    let mut s1 = Vec::with_capacity(n);
    let mut s2 = Vec::with_capacity(n);

    for _ in 0..n {
        s1.push(G1_GENERATOR.mul(&s_pow));
        s2.push(G2_GENERATOR.mul(&s_pow));

        s_pow = s_pow.mul(&s);
    }

    (s1, s2)
}

// n = 2^scale
pub fn generate_trusted_setup_eval_form(
    scale: usize,
    secret: [u8; 32usize],
) -> (Vec<FsG1>, Vec<FsG2>) {
    let mut s: FsFr = hash_to_bls_field(&secret);
    let n = 1 << scale;
    let root_of_unity = FsFr::from_u64_arr(&SCALE2_ROOT_OF_UNITY[scale]);

    let s_pow_n = s.pow(n);
    let numerator: FsFr = s_pow_n.sub(&FsFr::one());

    let mut s1 = Vec::with_capacity(n);
    let n_fr = FsFr::from_u64(n as u64);
    for i in 0..n {
        let root_powed = root_of_unity.pow(i);
        let denominator = &n_fr.mul(&s.sub(&root_powed));
        let lag = numerator.mul(&root_powed).div(denominator).unwrap();
        s1.push(G1_GENERATOR.mul(&lag));
    }
    let mut s2 = Vec::with_capacity(n);
    let mut s_pow = Fr::one();
    for _ in 0..n {
        s2.push(G2_GENERATOR.mul(&s_pow));
        s_pow = s_pow.mul(&s);
    }

    (s1, s2)
}
