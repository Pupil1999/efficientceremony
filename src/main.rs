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
    //test the dory proof generate
    let qap = ConstraintSystem::create_simple(24, 8, 3);
    let precom = DoryPrecomputation::generate_from_qap(&qap);
    let G = G1::generator();
    let H = G2::generator();
    let w1 = vec![G; 8];
    let w2 = vec![H; 8];
    let state = DoryStatement::create_state(&precom, &w1, &w2);
    DoryProof::generate(&w1, &w2, &precom, &state);
}