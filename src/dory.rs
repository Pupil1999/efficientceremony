use crate::ConstraintSystem;
use ark_ff::{UniformRand, Zero};
use ark_std::{self, One};
use ark_ff::Field;
use ark_ec::{pairing::{Pairing, PairingOutput}, bls12::Bls12, CurveGroup};
use ark_bls12_381::{Bls12_381, Config, G1Affine, G2Affine, Fr, G1Projective as G1, G2Projective as G2};
use sha2::{Sha256, Digest};
use ark_serialize::{CanonicalSerialize, CanonicalDeserialize};
use ark_ec::VariableBaseMSM;

type Gt = PairingOutput<Bls12<Config>>;

pub struct DoryStatement{
    pub x: (Gt, Gt, Gt)
}
#[allow(non_snake_case)]
impl DoryStatement {
    pub fn create_state(precom: &DoryPrecomputation, w1: &Vec<G1>, w2: &Vec<G2>) -> Self{
        let C = Bls12_381::multi_pairing(w1, w2);
        let D_1 = Bls12_381::multi_pairing(w1, &precom.combases.1);
        let D_2 = Bls12_381::multi_pairing(&precom.combases.0, w2);
        Self {x:(C, D_1, D_2)}
    }
}

#[allow(non_snake_case)]
pub struct DoryProof{
    pub D: Vec<(Gt, Gt, Gt, Gt)>,
    pub C: Vec<(Gt, Gt)>,
    pub m: usize,
    pub E: (G1Affine, G2Affine)
}
#[allow(non_snake_case)]
impl DoryProof {
    pub fn generate(w1: &Vec<G1>, w2: &Vec<G2>,
            setup: &DoryPrecomputation, statement: &DoryStatement) -> Self{
        assert_eq!(w1.len(), w2.len());
        let m = setup.m;
        let mut D = Vec::<(Gt, Gt, Gt, Gt)>::new();
        let mut C = Vec::<(Gt, Gt)>::new();

        let mut v1 = w1.clone();
        let mut v2 = w2.clone();
        let mut C_ = statement.x.0;
        let mut D_1 = statement.x.1;
        let mut D_2 = statement.x.2;
        for j in 0..m {
            let half = 2_usize.pow((m - j - 1) as u32);
            let d_1_left = Bls12_381::multi_pairing(&v1[0..half], &setup.combases.1[0..half]);
            let d_1_right = Bls12_381::multi_pairing(&v1[half..2*half], &setup.combases.1[0..half]);
            let d_2_left = Bls12_381::multi_pairing(&setup.combases.0[0..half], &v2[0..half]);
            let d_2_right = Bls12_381::multi_pairing(&setup.combases.0[0..half], &v2[half..2*half]);
            D.push((d_1_left, d_1_right, d_2_left, d_2_right));

            let beta: Fr = DoryOracle::fromGts(&[d_1_left, d_1_right, d_2_left, d_2_right]);
            let beta_inv = beta.inverse().unwrap();
            
            for k in 0..half * 2 {
                v1[k] = (v1[k] + setup.combases.0[k] * beta).into();
                v2[k] = (v2[k] + setup.combases.1[k] * beta_inv).into();
            }
            let c_positive = Bls12_381::multi_pairing(&v1[0..half], &v2[half..2*half]);
            let c_negative = Bls12_381::multi_pairing(&v1[half..2*half], &v2[0..half]);

            C.push((c_positive, c_negative));

            let alpha: Fr = DoryOracle::fromGts(&[c_positive, c_negative]);
            let alpha_inv = alpha.inverse().unwrap();
            for k in 0..half {
                v1[k] = (v1[k] * alpha + v1[k + half]).into();
                v2[k] = (v2[k] * alpha_inv + v2[k + half]).into();
            }

            C_ += Gt::msm(&[setup.xs[j], D_2, D_1, c_positive, c_negative], &[Fr::one(), beta, beta_inv, alpha, alpha_inv]).unwrap();
            D_1 = Gt::msm(&[d_1_left, d_1_right, setup.precoms[j].0, setup.precoms[j].1], &[alpha, Fr::one(), alpha*beta, beta]).unwrap();
            D_2 = Gt::msm(&[d_2_left, d_2_right, setup.precoms[j].2, setup.precoms[j].3], &[alpha_inv, Fr::one(), alpha_inv*beta_inv, beta_inv]).unwrap();
        }
    
        Self {D, C, m, E: (v1[0].into_affine(), v2[0].into_affine())}
    }

    pub fn verify(state: &DoryStatement, setup: &DoryPrecomputation, proof: &DoryProof) -> bool{
        let (C, D_1, D_2) = state.x;
        let m = setup.m;
        let mut betas = Vec::<Fr>::with_capacity(m);
        let mut betas_inv = Vec::<Fr>::with_capacity(m);
        let mut alphas = Vec::<Fr>::with_capacity(m);
        let mut alphas_inv = Vec::<Fr>::with_capacity(m);
        for elm in proof.D.iter() {
            let beta = DoryOracle::fromGts(&[elm.0, elm.1, elm.2, elm.3]);
            betas.push(beta);
            betas_inv.push(beta.inverse().unwrap());
        }
        for elm in proof.C.iter() {
            let alpha = DoryOracle::fromGts(&[elm.0, elm.1]);
            alphas.push(alpha);
            alphas_inv.push(alpha.inverse().unwrap());
        }
        let mut C_fin = C + setup.sigmax;
        let mut D_1_fin = D_1;
        let mut D_2_fin = D_2;

        for i in 0..m {
            C_fin += Gt::msm(&[D_2_fin, D_1_fin, proof.C[i].0, proof.C[i].1], &[betas[i], betas_inv[i], alphas[i], alphas_inv[i]]).unwrap();
            D_1_fin = Gt::msm(&[proof.D[i].0, proof.D[i].1, setup.precoms[i].0, setup.precoms[i].1], &[alphas[i], Fr::one(), alphas[i]*betas[i], betas[i]]).unwrap();
            D_2_fin = Gt::msm(&[proof.D[i].2, proof.D[i].3, setup.precoms[i].2, setup.precoms[i].3], &[alphas_inv[i], Fr::one(), alphas_inv[i]*betas_inv[i], betas_inv[i]]).unwrap();
        }

        let mut rng = ark_std::test_rng();
        let d = Fr::rand(&mut rng);
        let left = Bls12_381::pairing(proof.E.0 + setup.combases.0[0] * d, proof.E.1 + setup.combases.1[0] * d.inverse().unwrap());
        let right = C_fin + D_2_fin*d + D_1_fin*d.inverse().unwrap();

        left == right
    }


}

#[allow(non_snake_case)]
pub struct DoryPrecomputation {
    pub m: usize,
    pub combases: (Vec<G1>, Vec<G2>),
    pub precoms: Vec<(Gt, Gt, Gt, Gt)>,
    pub xs: Vec<Gt>,
    pub sigmax: Gt
}
#[allow(non_snake_case)]
impl DoryPrecomputation {
    pub fn generate_from_qap(qap: &ConstraintSystem) -> Self{
        //get round times
        let (_, n, _) = qap.get_size();

        //m = log(n)
        let m = ark_std::log2(n) as usize;

        let mut precoms: Vec<(Gt, Gt, Gt, Gt)> = Vec::new();

        let mut rng = ark_std::test_rng();
        let mut n_iter;
        let mut Gama_left = Vec::<G1>::with_capacity(n);
        let mut Gama_right = Vec::<G2>::with_capacity(n);
        let mut xs = Vec::<Gt>::with_capacity(m + 1);
        for _i in 0..n {
            Gama_left.push(G1::rand(&mut rng));
            Gama_right.push(G2::rand(&mut rng));
        }
        let combases = (Gama_left.clone(), Gama_right.clone());
        let mut sigmax: Gt = Gt::zero();
        for i in (0..=m).rev() {
            n_iter = 2_usize.pow(i as u32);
            let bmo = Bls12_381::multi_pairing(&Gama_left[0..n_iter], &Gama_right[0..n_iter]);
            xs.push(bmo);
            sigmax += bmo;
        }
        for i in (0..m).rev() {
            let half = 2_usize.pow(i as u32);
            let delta_1_left = Bls12_381::multi_pairing(&Gama_left[0..half], &Gama_right[0..half]);
            let delta_1_right = Bls12_381::multi_pairing(&Gama_left[half..2*half], &Gama_right[0..half]);
            let delta_2_left = Bls12_381::multi_pairing(&Gama_left[0..half], &Gama_right[0..half]);
            let delta_2_right = Bls12_381::multi_pairing(&Gama_left[0..half], &Gama_right[half..2*half]);
            precoms.push((delta_1_left, delta_1_right, delta_2_left, delta_2_right));            
        }

        Self {m, combases, precoms, xs, sigmax}
    }
}

pub struct DoryOracle;
#[allow(non_snake_case)]
impl DoryOracle {
    pub fn fromGts(gt: &[Gt]) -> Fr {
        let mut buffers = Vec::new();
        for elm in gt {
            let mut buffer = Vec::new();
            elm.serialize_compressed(&mut buffer).unwrap();
            buffers.push(buffer);
        }
        let content = buffers.concat();

        let mut hasher = Sha256::new();
        hasher.update(&content);
        
        let mut u8array: [u8; 32] = hasher.finalize().as_slice().try_into().unwrap();
        u8array[31] = 0;
        let rn: Fr = Fr::deserialize_compressed(&u8array[..]).unwrap();
        rn
    }
}