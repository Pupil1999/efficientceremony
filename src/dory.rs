use ark_ff::{UniformRand, Zero};
use ark_std::{self, One};
use ark_ff::Field;
use ark_ec::{pairing::{Pairing, PairingOutput}, bls12::Bls12, CurveGroup, AffineRepr, Group};
use ark_bls12_381::{Bls12_381, Config, G1Affine, G2Affine, Fr, G1Projective as G1, G2Projective as G2};
use sha2::{Sha256, Digest};
use ark_serialize::{CanonicalSerialize, CanonicalDeserialize};
use ark_ec::VariableBaseMSM;

type Gt = PairingOutput<Bls12<Config>>;

#[allow(non_snake_case)]
pub struct DoryStatement{
    pub x: (Gt, Gt, Gt),
    pub E: (G1Affine, G2Affine)
}
#[allow(non_snake_case)]
impl DoryStatement {
    pub fn create_state(precom: &DoryPrecomputation, w1: &Vec<G1Affine>, w2: &Vec<G2Affine>, s1: &Vec<Fr>, s2: &Vec<Fr>) -> Self{
        let C = Bls12_381::multi_pairing(w1, w2);
        let D_1 = Bls12_381::multi_pairing(w1, &precom.combases.1[0..w1.len()]);
        let D_2 = Bls12_381::multi_pairing(&precom.combases.0[0..w2.len()], w2);
        let m = precom.m;
        let mut E_1 = G1::zero();
        let mut E_2 = G2::zero();

        let s = 2_usize.pow(m as u32);
        if s > s1.len() {
            for i in 0..s1.len() {
                E_1 += w1[i] * s2[i];
                E_2 += w2[i] * s1[i];
            }
        } else {
            for i in 0..2_usize.pow(m as u32) {
                E_1 += w1[i] * s2[i];
                E_2 += w2[i] * s1[i];
            }
        }
        Self {x:(C, D_1, D_2), E:(E_1.into_affine(), E_2.into_affine())}
    }
}

#[allow(non_snake_case)]
pub struct DoryProof{
    pub D: Vec<(Gt, Gt, Gt, Gt)>,
    pub C: Vec<(Gt, Gt)>,
    pub E_beta: Vec<(G1Affine, G2Affine)>,
    pub E: Vec<(G1Affine, G1Affine, G2Affine, G2Affine)>,
    pub m: usize,
    pub v_fin: (G1Affine, G2Affine),
    pub s_fin: (Fr, Fr)
}
#[allow(non_snake_case)]
impl DoryProof {
    pub fn generate(witness1: &Vec<G1Affine>, witness2: &Vec<G2Affine>, scalar1: &Vec<Fr>, scalar2: &Vec<Fr>,
            setup: &DoryPrecomputation, statement: &DoryStatement) -> Self {
        assert_eq!(witness1.len(), witness2.len());
        let m = setup.m;
        let mut D = Vec::<(Gt, Gt, Gt, Gt)>::with_capacity(m);
        let mut C = Vec::<(Gt, Gt)>::with_capacity(m);
        let mut E_beta = Vec::<(G1Affine, G2Affine)>::with_capacity(m);
        let mut E = Vec::<(G1Affine, G1Affine, G2Affine, G2Affine)>::with_capacity(m);

        let mut v1 = witness1.clone();
        let mut v2 = witness2.clone();
        let mut s1: Vec<Fr> = scalar1.clone();
        let mut s2: Vec<Fr> = scalar2.clone();
        let mut C_ = statement.x.0;
        let mut D_1 = statement.x.1;
        let mut D_2 = statement.x.2;
        let mut E_1 = statement.E.0.into_group();
        let mut E_2 = statement.E.1.into_group();
        for j in 0..m {
            let half = 2_usize.pow((m - j - 1) as u32);
            let d_1_left = Bls12_381::multi_pairing(&v1[0..half], &setup.combases.1[0..half]);
            let d_1_right = Bls12_381::multi_pairing(&v1[half..2*half], &setup.combases.1[0..half]);
            let d_2_left = Bls12_381::multi_pairing(&setup.combases.0[0..half], &v2[0..half]);
            let d_2_right = Bls12_381::multi_pairing(&setup.combases.0[0..half], &v2[half..2*half]);
            D.push((d_1_left, d_1_right, d_2_left, d_2_right));

            let mut e_1_beta = G1::zero();
            let mut e_2_beta = G2::zero();
            for i in 0..half*2 {
                e_1_beta += setup.combases.0[i].into_group() * s2[i];
                e_2_beta += setup.combases.1[i].into_group() * s1[i];
            }
            E_beta.push((e_1_beta.into_affine(), e_2_beta.into_affine()));

            
            let beta: Fr = DoryOracle::fromGts(&[d_1_left, d_1_right, d_2_left, d_2_right]);
            let beta_inv = beta.inverse().unwrap();
            for k in 0..half * 2 {
                v1[k] = (v1[k] + setup.combases.0[k] * beta).into();
                v2[k] = (v2[k] + setup.combases.1[k] * beta_inv).into();
            }
            let c_positive = Bls12_381::multi_pairing(&v1[0..half], &v2[half..2*half]);
            let c_negative = Bls12_381::multi_pairing(&v1[half..2*half], &v2[0..half]);
            C.push((c_positive, c_negative));

            let mut e_1_positive = G1::zero();
            let mut e_1_negative = G1::zero();
            let mut e_2_positive = G2::zero();
            let mut e_2_negative = G2::zero();
            for i in 0..half{
                e_1_positive += v1[i] * s2[i + half];
                e_1_negative += v1[i + half] * s2[i];
                e_2_positive += v2[i + half] * s1[i];
                e_2_negative += v2[i] * s1[i + half];
            }
            E.push((e_1_positive.into_affine(), e_1_negative.into_affine(), e_2_positive.into_affine(), e_2_negative.into_affine()));

            let alpha: Fr = DoryOracle::fromGts(&[c_positive, c_negative]);
            let alpha_inv = alpha.inverse().unwrap();
            //MSM in this library seems to have some problems, so we use common group mul and add ops.
            for k in 0..half {
                v1[k] = (v1[k] * alpha + v1[k + half]).into();
                v2[k] = (v2[k] * alpha_inv + v2[k + half]).into();
                s1[k] = s1[k] * alpha + s1[k + half];
                s2[k] = s2[k] * alpha_inv + s2[k + half];
            }

            C_ += Gt::msm(&[setup.xs[j], D_2, D_1, c_positive, c_negative], &[Fr::one(), beta, beta_inv, alpha, alpha_inv]).unwrap();
            D_1 = Gt::msm(&[d_1_left, d_1_right, setup.precoms[j].0, setup.precoms[j].1], &[alpha, Fr::one(), alpha*beta, beta]).unwrap();
            D_2 = Gt::msm(&[d_2_left, d_2_right, setup.precoms[j].2, setup.precoms[j].3], &[alpha_inv, Fr::one(), alpha_inv*beta_inv, beta_inv]).unwrap();
        
            E_1 += G1::msm(&[e_1_beta.into_affine(), e_1_positive.into_affine(), e_1_negative.into_affine()], &[beta, alpha, alpha_inv]).unwrap();
            E_2 += G2::msm(&[e_2_beta.into_affine(), e_2_positive.into_affine(), e_2_negative.into_affine()], &[beta_inv, alpha, alpha_inv]).unwrap();
        }

        Self {D, C, m, E_beta, E, v_fin: (v1[0], v2[0]), s_fin: (s1[0], s2[0])}
    }

    pub fn verify(state: &DoryStatement, setup: &DoryPrecomputation, proof: &DoryProof) -> bool{
        let (C, D_1, D_2) = state.x;
        let (E_1, E_2) = state.E;
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
        let mut E_1_fin = E_1.into_group();
        let mut E_2_fin = E_2.into_group();

        for i in 0..m {
            C_fin += Gt::msm(&[D_2_fin, D_1_fin, proof.C[i].0, proof.C[i].1], &[betas[i], betas_inv[i], alphas[i], alphas_inv[i]]).unwrap();
            D_1_fin = Gt::msm(&[proof.D[i].0, proof.D[i].1, setup.precoms[i].0, setup.precoms[i].1], &[alphas[i], Fr::one(), alphas[i]*betas[i], betas[i]]).unwrap();
            D_2_fin = Gt::msm(&[proof.D[i].2, proof.D[i].3, setup.precoms[i].2, setup.precoms[i].3], &[alphas_inv[i], Fr::one(), alphas_inv[i]*betas_inv[i], betas_inv[i]]).unwrap();
            E_1_fin += G1::msm(&[proof.E_beta[i].0, proof.E[i].0, proof.E[i].1], &[betas[i], alphas[i], alphas_inv[i]]).unwrap();
            E_2_fin += G2::msm(&[proof.E_beta[i].1, proof.E[i].2, proof.E[i].3], &[betas_inv[i], alphas[i], alphas_inv[i]]).unwrap();
        }
        let G = G1::generator();
        let H = G2::generator();

        let mut rng = ark_std::test_rng();
        let r = Fr::rand(&mut rng);
        let r_inv = r.inverse().unwrap();
        let v1 = proof.v_fin.0 + G*r*proof.s_fin.0;
        let v2 = proof.v_fin.1 + H*r_inv*proof.s_fin.1;

        C_fin += Bls12_381::multi_pairing(&[G*proof.s_fin.0*proof.s_fin.1, G*r, E_1_fin*r_inv], &[H, E_2_fin, H]);
        D_1_fin += Bls12_381::pairing(G, setup.combases.1[0]*proof.s_fin.0*r);
        D_2_fin += Bls12_381::pairing(setup.combases.0[0]*proof.s_fin.1*r_inv, H);

        let d = Fr::rand(&mut rng);
        let d_inv = d.inverse().unwrap();
        let left = Bls12_381::pairing(v1 + setup.combases.0[0] * d, v2 + setup.combases.1[0] * d_inv);
        let right = C_fin + D_2_fin*d + D_1_fin*d_inv;

        left == right
    }


}

#[allow(non_snake_case)]
pub struct DoryPrecomputation {
    pub m: usize,
    pub combases: (Vec<G1Affine>, Vec<G2Affine>),
    pub precoms: Vec<(Gt, Gt, Gt, Gt)>,
    pub xs: Vec<Gt>,
    pub sigmax: Gt
}
#[allow(non_snake_case)]
impl DoryPrecomputation {
    pub fn generate_with_size(n: usize) -> Self{
        //m = log(n)
        let m = ark_std::log2(n) as usize;

        let mut precoms: Vec<(Gt, Gt, Gt, Gt)> = Vec::new();

        let mut rng = ark_std::test_rng();
        let mut n_iter;
        let mut Gama_left = Vec::<G1Affine>::with_capacity(n);
        let mut Gama_right = Vec::<G2Affine>::with_capacity(n);
        let mut xs = Vec::<Gt>::with_capacity(m + 1);
        for _i in 0..n {
            Gama_left.push(G1::rand(&mut rng).into_affine());
            Gama_right.push(G2::rand(&mut rng).into_affine());
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