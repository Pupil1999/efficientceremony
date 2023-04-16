use crate::ConstraintSystem;
use ark_ec::{pairing::{Pairing, PairingOutput}, bls12::Bls12};
use ark_bls12_381::{Bls12_381, Config};

type Gt = PairingOutput<Bls12<Config>>;

pub struct DoryProof{
    
}

pub struct DoryPrecomputation {
    _m: usize,
    pub precoms: Vec<(Gt, Gt, Gt, Gt)>
}
impl DoryPrecomputation {
    pub fn generate_from_qap(&mut self, qap: &ConstraintSystem) {
        //get round times
        let (_, n, _) = qap.get_size();
        for i in 0..64 {
            if 2_usize.pow(i) >= n {
                self._m = 2_usize.pow(i);
                break;
            }
        }


    }
}