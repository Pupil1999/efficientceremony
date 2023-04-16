use std::ops::Mul;

use crate::constraintsystem::ConstraintSystem;
use ark_ec::Group;
use ark_std::UniformRand;
use ark_bls12_381::{G1Affine, G2Affine, G1Projective as G1, G2Projective as G2, Fr};

pub enum PHASE {
    ONE,
    TWO
}

pub struct Trapdoor {
    pub tau: Fr,
    pub alpha: Fr,
    pub beta: Fr,
    pub delta: Fr,
}
impl Trapdoor {
    pub fn get_random_trapdoors() -> Self {
        let mut rng = ark_std::test_rng();

        let tau = Fr::rand(&mut rng);
        let alpha = Fr::rand(&mut rng);
        let beta = Fr::rand(&mut rng);
        let delta = Fr::rand(&mut rng);

        Self {tau, alpha, beta, delta}
    }
}

pub type U = (Vec<(G1Affine, G2Affine)>, Vec<(G1Affine, G2Affine, G1Affine, G2Affine)>);
pub type S = (G1Affine, G2Affine, Vec<G1Affine>, Vec<G1Affine>);
pub struct SRS {
    pub u: U,
    pub s: S,
    pub n: usize,
    pub m: usize,
    pub l: usize,
}

#[allow(non_snake_case)]
impl SRS {
    pub fn create_from_qap(qap: &ConstraintSystem) -> Self {
        let (m, n, l) = qap.get_size();
        let G = G1::generator().into();
        let H = G2::generator().into();
        let mut x_array: Vec<(G1Affine, G2Affine)> = Vec::new();
        let mut ab_array: Vec<(G1Affine, G2Affine, G1Affine, G2Affine)> = Vec::new();
        let mut gd_1: Vec<G1Affine> = Vec::new();
        let mut gd_2: Vec<G1Affine> = Vec::new();

        for _i in 0..=(n - 1) {
            x_array.push((G, H));
            ab_array.push((G, H, G, H));
        }

        for _i in 0..=(n - 2) {
            x_array.push((G, H));
        }

        for _i in 1..=(m - l) {
            gd_1.push(G);
        }

        for _i in 0..=(n - 2) {
            gd_2.push(G);
        }

        Self {
            u: (x_array, ab_array),
            s: (G, H, gd_1, gd_2),
            n,
            m,
            l
        }
    }

    pub fn update(&mut self, phase: PHASE, trapdoors: &Trapdoor){
        let (m, n, l) = (self.m, self.n, self.l);
        let mut temp: Fr = trapdoors.tau;

        match phase {
            PHASE::ONE => {
                for i in 1..=(2 * n - 2) {
                    self.u.0[i].0 = self.u.0[i].0.mul(temp).into();
                    self.u.0[i].1 = self.u.0[i].1.mul(temp).into();
                    temp *= trapdoors.tau;
                }

                temp = trapdoors.tau;

                //Update the first elements of G_alpha, H_alpha, G_beta and H_beta arrays.
                self.u.1[0].0 = self.u.1[0].0.mul(&trapdoors.alpha).into();
                self.u.1[0].1 = self.u.1[0].1.mul(&trapdoors.alpha).into();
                self.u.1[0].2 = self.u.1[0].2.mul(&trapdoors.beta).into();
                self.u.1[0].3 = self.u.1[0].3.mul(&trapdoors.beta).into();

                for i in 1..=n - 1 {
                    self.u.1[i].0 = self.u.1[i].0.mul(&trapdoors.alpha).into();
                    self.u.1[i].1 = self.u.1[i].1.mul(&trapdoors.alpha).into();
                    self.u.1[i].2 = self.u.1[i].2.mul(&trapdoors.beta).into();
                    self.u.1[i].3 = self.u.1[i].3.mul(&trapdoors.beta).into();

                    self.u.1[i].0 = self.u.1[i].0.mul(&temp).into();
                    self.u.1[i].1 = self.u.1[i].1.mul(&temp).into();
                    self.u.1[i].2 = self.u.1[i].2.mul(&temp).into();
                    self.u.1[i].3 = self.u.1[i].3.mul(&temp).into();
                    temp *= trapdoors.tau;
                }
            }

            PHASE::TWO => {
                for _i in 0..=(m - l) {
                    
                }
            }
        }
    }

}