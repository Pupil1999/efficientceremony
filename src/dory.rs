use crate::ConstraintSystem;
use std::io::Write;
use ark_ff::UniformRand;
use ark_std;
use ark_ff::{Field, BigInteger256};
use ark_ec::{pairing::{Pairing, PairingOutput}, bls12::Bls12, CurveGroup};
use ark_bls12_381::{Bls12_381, Config, G1Affine, G2Affine, Fr, G1Projective as G1, G2Projective as G2};
use sha2::{Sha256, Digest};
use ark_serialize::{CanonicalSerialize, CanonicalDeserialize};

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

        let mut rng = ark_std::test_rng();
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

            let beta = Fr::rand(&mut rng);
            let beta_inv = beta.inverse().unwrap();
            
            for k in 0..half * 2 {
                v1[k] = (v1[k] + setup.combases.0[k] * beta).into();
                v2[k] = (v2[k] + setup.combases.1[k] * beta_inv).into();
            }
            let c_positive = Bls12_381::multi_pairing(&v1[0..half], &v2[half..2*half]);
            let c_negative = Bls12_381::multi_pairing(&v1[half..2*half], &v2[0..half]);

            C.push((c_positive, c_negative));

            let alpha = Fr::rand(&mut rng);
            let alpha_inv = alpha.inverse().unwrap();
            for k in 0..half {
                v1[k] = (v1[k] * alpha + v1[k + half]).into();
                v2[k] = (v2[k] * alpha_inv + v2[k + half]).into();
            }

            C_ = C_ + setup.xs[j] + D_2*beta + D_1*beta_inv + c_positive*alpha + c_negative*alpha_inv;
            D_1 = d_1_left*alpha + d_1_right + setup.precoms[j].0*alpha*beta + setup.precoms[j].1*beta;
            D_2 = d_2_left*alpha_inv + d_2_right + setup.precoms[j].2*alpha_inv*beta_inv + setup.precoms[j].3*beta_inv;
            //assert_eq!(Bls12_381::multi_pairing(&v1[0..half], &v2[0..half]), C_);
            //assert_eq!(Bls12_381::multi_pairing(&v1[0..half], &setup.combases.1[0..half]), D_1);
            //assert_eq!(Bls12_381::multi_pairing(&setup.combases.0[0..half], &v2[0..half]), D_2);
            // if j == m - 1{
            //     let d = Fr::rand(&mut rng);
            //     let left_out = Bls12_381::pairing(&v1[0] + setup.combases.0[0] * d, &v2[0] + setup.combases.1[0] * d.inverse().unwrap());
            //     let right_out = setup.xs[setup.m] + C_ + D_2 * d + D_1 * d.inverse().unwrap();
            //     assert_eq!(left_out, right_out);
            // }
        }
    
        Self {D, C, m, E: (v1[0].into_affine(), v2[0].into_affine())}
    }

    pub fn verify(state: &DoryStatement, setup: &DoryPrecomputation, proof: &DoryProof) {

    }
}

#[allow(non_snake_case)]
pub struct DoryPrecomputation {
    pub m: usize,
    pub combases: (Vec<G1>, Vec<G2>),
    pub precoms: Vec<(Gt, Gt, Gt, Gt)>,
    pub xs: Vec<Gt>
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
        for i in (0..=m).rev() {
            n_iter = 2_usize.pow(i as u32);
            xs.push(Bls12_381::multi_pairing(&Gama_left[0..n_iter], &Gama_right[0..n_iter]));
        }
        for i in (0..m).rev() {
            let half = 2_usize.pow(i as u32);
            let delta_1_left = Bls12_381::multi_pairing(&Gama_left[0..half], &Gama_right[0..half]);
            let delta_1_right = Bls12_381::multi_pairing(&Gama_left[half..2*half], &Gama_right[0..half]);
            let delta_2_left = Bls12_381::multi_pairing(&Gama_left[0..half], &Gama_right[0..half]);
            let delta_2_right = Bls12_381::multi_pairing(&Gama_left[0..half], &Gama_right[half..2*half]);
            precoms.push((delta_1_left, delta_1_right, delta_2_left, delta_2_right));            
        }

        Self {m, combases, precoms, xs}
    }
}

pub struct DoryOracle;
#[allow(non_snake_case)]
impl DoryOracle {
    pub fn fromGts(gt: &Vec<Gt>) -> Fr {
        //let mut str: String = Default::default();
        // for i in gt {
        //     str += &i.to_string();
        // }

        // let content = str.as_bytes();
        // let mut hasher = Sha256::new();
        // hasher.update(content);
        
        // let u8array: [u8; 64] = hasher.finalize().as_slice().try_into().unwrap();

        // let buffer = u8array.into();
        // let b = BigInteger256::new(buffer);

        
        // let mut hasher = Sha256::new();
        // let mut buffer = Vec::new();
        // let mut bts: &[u8] = gt[0].serialize_compressed(&mut buffer).unwrap().into();
        // hasher.update(& bts);

        Fr::from(10)
    }
}