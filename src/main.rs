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

    let G = G1::generator();
    let H = G2::generator();

    let s_13 = Fr::from(13);
    let s_12 = Fr::from(12);
    let s_156 = Fr::from(156);

    let G_13 = G * s_13;
    let H_12 = H * s_12;
    let G_156 = G * s_156;
    
    let gt_156 = Bls12_381::pairing(G_13, H_12);
    let gt_156_2 = gt_156.mul(s_13);
    assert_eq!(gt_156_2, Bls12_381::pairing(G_156, H * s_13));
}