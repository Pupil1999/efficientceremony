pub mod constraintsystem;
pub mod srs;
pub mod dory;

use constraintsystem::*;
use srs::*;
use dory::*;
use ark_std::ops::Mul;
use ark_ff::{fields::Field, UniformRand};
use ark_ec::{Group, AffineRepr, pairing::Pairing};
use ark_bls12_381::{Bls12_381, G1Affine, G2Affine, G1Projective as G1, G2Projective as G2, Fr};

#[allow(non_snake_case)]
fn main() {
    
}