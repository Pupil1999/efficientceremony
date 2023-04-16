use crate::ConstraintSystem;
use ark_ec::{pairing::{Pairing, PairingOutput}, bls12::Bls12};
use ark_bls12_381::{Bls12_381, Config};

type Gt = PairingOutput<Bls12<Config>>;

pub struct DoryProof{
    
}

pub struct DoryPrecomputation {
    _m: u64,
    pub precoms: Vec<(Gt, Gt, Gt, Gt)>
}
impl DoryPrecomputation {
    pub fn generate_from_qap(&mut self, QAP: &ConstraintSystem) {
        let (_, n, _) = QAP.get_size();
    }
}