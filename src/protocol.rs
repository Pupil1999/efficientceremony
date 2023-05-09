use crate::srs::*;
use crate::dory::*;
use ark_std::{Zero, One, UniformRand};
use ark_ff::Field;
use ark_bls12_381::Bls12_381;
use ark_bls12_381::{G1Affine, G2Affine, G1Projective as G1, G2Projective as G2, Fr, Config};
use ark_ec::{pairing::{Pairing, PairingOutput}, bls12::Bls12, CurveGroup, AffineRepr, VariableBaseMSM};
type Gt = PairingOutput<Bls12<Config>>;

pub struct Groth16Ceremony {
    _phase: PHASE,
    pub qap_size: (usize, usize, usize),
    pub x_array_info: (Gt, Gt, G1Affine, G2Affine),
    pub delta_com: Gt,
    pub delta_array_info: (Gt, Gt, G1Affine, G2Affine),
}
#[allow(non_snake_case)]
impl Groth16Ceremony {
    //Initialize the ceremony with the basic parameters of QAP system
    pub fn create_from_origin(crs: &SRS, precom: &DoryPrecomputation) -> Self {
        let (m, n, l) = (crs.m, crs.n, crs.l);
        let mut G_xs = Vec::new();
        let mut H_xs = Vec::new();
        for elm in &crs.u.0 {
            G_xs.push(elm.0);
            H_xs.push(elm.1);
        }
        //Four elements used to check the struct of x array.
        let Cx_L = Bls12_381::multi_pairing(&G_xs[0..2*n], &precom.combases.1[0..2*n]);
        let Cx_R = Bls12_381::multi_pairing(&precom.combases.0[0..2*n], &H_xs[0..2*n]);
        
        let Px_1 = crs.u.0[0].0;
        let Qx_1 = crs.u.0[0].1;

        Self { 
            _phase: PHASE::ONE,
            qap_size: (m, n, l),
            x_array_info: (Cx_L, Cx_R, Px_1, Qx_1),
            delta_com: Gt::zero(),
            delta_array_info: (Gt::zero(), Gt::zero(), G1Affine::zero(), G2Affine::zero())
        }
    }

    pub fn update_1(&self, crs: &mut SRS, trapdoors: &Trapdoor, 
                precom: &DoryPrecomputation) -> Contribution1 {
        SRS::update(crs, &PHASE::ONE, trapdoors);
        let n = crs.n;
        let mut G_xs = Vec::new();
        let mut H_xs = Vec::new();
        for elm in &crs.u.0 {
            G_xs.push(elm.0);
            H_xs.push(elm.1);
        }

        //Four elements used to check the struct of x array.
        let Cx_L = Bls12_381::multi_pairing(&G_xs, &precom.combases.1);
        let Cx_R = Bls12_381::multi_pairing(&precom.combases.0, &H_xs);
        let Px_1 = crs.u.0[0].0;
        let Qx_1 = crs.u.0[0].1;
        let Px_last = crs.u.0[2*n-1].0;
        let Qx_last = crs.u.0[2*n-1].1;

        let mut rng = ark_std::test_rng();
        let mut randvec_1 = Vec::new();
        let mut randvec_2 = Vec::new();
        //For simplicity, we just sampled a random number instead of querying the oracle.
        let rand1 = Fr::rand(&mut rng);
        let rand2 = Fr::rand(&mut rng);
        let mut temp1 = Fr::one();
        let mut temp2 = Fr::one();
        for _i in 0..2*n {
            randvec_1.push(temp1);
            temp1 *= rand1;
            randvec_2.push(temp2);
            temp2 *= rand2;
        }
        let Ex_1 = G1::msm(&G_xs[0..2*n], &randvec_1).unwrap().into_affine();
        let Ex_2 = G2::msm(&H_xs[0..2*n], &randvec_2).unwrap().into_affine();

        let statement = DoryStatement::create_state(&precom, &G_xs, &H_xs, &randvec_2 as &Vec<Fr>, &randvec_1 as &Vec<Fr>);
        let proof_1 = DoryProof::generate(&G_xs, &H_xs, &randvec_2 as &Vec<Fr>, &randvec_1 as &Vec<Fr>, &precom, &statement);

        Contribution1 {
            x_array_info: (Cx_L, Cx_R, Px_1, Qx_1),
            x_array_aux_info: (Px_last, Qx_last),
            x_rand: (rand1, rand2),
            x_dory_extend_statement: (Ex_1, Ex_2),
            x_dory_states: statement,
            x_dory_proofs: proof_1
        }
    }

    pub fn update_2(&self, crs: &mut SRS, trapdoors: &Trapdoor, 
        precom: &DoryPrecomputation) -> Contribution2 {
        SRS::update(crs, &PHASE::TWO, &trapdoors);
        
        let G = G1Affine::generator();
        let mut Gs = Vec::new();
        Gs.push(G);
        for elm in &crs.s.2 {
            Gs.push(*elm);
        }
        for elm in &crs.s.3 {
            Gs.push(*elm);
        }
        let mut Hs = Vec::new();
        Hs.push(crs.s.1);
        for _i in 1..Gs.len() {
            Hs.push(G2Affine::zero());
        }
        let presize = 2_usize.pow(precom.m as u32);
        while Gs.len() < presize {
            Gs.push(G1Affine::zero());
            Hs.push(G2Affine::zero());
        }
        let mut rng = ark_std::test_rng();
        let rand = Fr::rand(&mut rng);
        let mut temp = rand;
        let mut s1 = Vec::new();
        for _i in 0..Gs.len() {
            s1.push(temp);
            temp *= rand;
        }
        let s2 = s1.clone();
        let state = DoryStatement::create_state(&precom, &Gs, &Hs, &s1, &s1);
        let delta_dory_proof = DoryProof::generate(&Gs, &Hs, &s2, &s1, precom, &state);

        Contribution2 {
            delta_array_info: (state.x.0, state.x.1, crs.s.0, crs.s.1),
            delta_array_aux_info: (state.x.2, state.E.1),
            delta_dory_extend_statement: (state.E.0, state.E.0),
            delta_rand: rand,
            delta_dory_state: state,
            delta_dory_proof
        }
    }

    pub fn simple_change_phase(&mut self, crs: &mut SRS, precom: &DoryPrecomputation) {
        self._phase = PHASE::TWO;
        SRS::simple_change_phase(crs);
        let mut origin = Vec::new();
        origin.push(G1Affine::generator());
        for elm in &crs.s.2 {
            origin.push(*elm);
        }
        for elm in &crs.s.3 {
            origin.push(*elm);
        }
        self.delta_com = Bls12_381::multi_pairing(&origin, &precom.combases.1[0..origin.len()]);
    }

    pub fn verify_1(&mut self, contri: &Contribution1, precom: &DoryPrecomputation) -> bool {
        //Check the structure of x array. It has two parts.
        let n = self.qap_size.1;
        let G = G1Affine::generator();
        let H = G2Affine::generator();
        let mut left = Bls12_381::pairing(&contri.x_dory_extend_statement.0, &(contri.x_array_info.3*contri.x_rand.0 - H));
        let rand1_2n = contri.x_rand.0.pow([2*n as u64, 0, 0, 0]);
        let mut right = Bls12_381::pairing(&(contri.x_array_aux_info.0*rand1_2n - G), &contri.x_array_info.3);
        let out1 = left==right;

        let rand2_2n = contri.x_rand.1.pow([2*n as u64, 0, 0, 0]);
        left = Bls12_381::pairing(&(contri.x_array_info.2*contri.x_rand.1 - G), &contri.x_dory_extend_statement.1);
        right = Bls12_381::pairing(&contri.x_array_info.2, &(contri.x_array_aux_info.1*rand2_2n - H));
        let out2 = left==right;
        let out3 = DoryProof::verify(&contri.x_dory_states, &precom, &contri.x_dory_proofs);
        
        out1 && out2 && out3
    }

    pub fn verify_2(&mut self, contri: &Contribution2, precom: &DoryPrecomputation) -> bool {
        let left = Bls12_381::pairing(&contri.delta_dory_extend_statement.0, &contri.delta_array_info.3);
        let G = G1Affine::generator();
        let H = G2Affine::generator();
        let right = Bls12_381::multi_pairing(&[contri.delta_dory_extend_statement.1, (G * Fr::from(-1)).into()], &[contri.delta_array_info.3, H]) + Bls12_381::pairing(&G, &H);
        let out1 = left==right;

        let out2 = DoryProof::verify(&contri.delta_dory_state, precom, &contri.delta_dory_proof);

        out1 && out2
    }
}

pub struct Contribution1 {
    pub x_array_info: (Gt, Gt, G1Affine, G2Affine),
    pub x_array_aux_info: (G1Affine, G2Affine),
    pub x_rand: (Fr, Fr),
    pub x_dory_extend_statement: (G1Affine, G2Affine),
    pub x_dory_states: DoryStatement,
    pub x_dory_proofs: DoryProof,
}

pub struct Contribution2 {
    pub delta_array_info: (Gt, Gt, G1Affine, G2Affine),
    pub delta_array_aux_info: (Gt, G2Affine),
    pub delta_dory_extend_statement: (G1Affine, G1Affine),
    pub delta_rand: Fr,
    pub delta_dory_state: DoryStatement,
    pub delta_dory_proof: DoryProof
}