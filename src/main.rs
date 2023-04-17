pub mod constraintsystem;
pub mod srs;
pub mod dory;

use constraintsystem::*;
use dory::{DoryPrecomputation, DoryStatement, DoryProof};
use srs::*;
use ark_std;
use ark_std::ops::Mul;
use ark_ff::{fields::Field, UniformRand};
use ark_ec::{Group, AffineRepr, pairing::Pairing};
use ark_bls12_381::{Bls12_381, G1Affine, G2Affine, G1Projective as G1, G2Projective as G2, Fr};

#[allow(non_snake_case)]
fn main() {
    let G = G1::generator();
    
    let bi = G.to_string();
    let bytes = bi.as_bytes();

    println!("{:?}", bytes);
}