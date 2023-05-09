pub mod constraintsystem;
pub mod srs;
pub mod dory;
pub mod protocol;

use constraintsystem::*;
use dory::DoryPrecomputation;
use srs::*;
use std::time::Instant;
use protocol::*;

#[allow(non_snake_case)]
fn main() {
    for m in 3..=8 {
        let scale = 2_usize.pow(m as u32);
        let qap = ConstraintSystem::create_simple(2*scale, scale, 3);
        println!("Crs size is {} in polynomial.", scale);
        let mut crs = SRS::create_from_qap(&qap);
        let precom1 = DoryPrecomputation::generate_with_size(2*scale);
        let precom2 = DoryPrecomputation::generate_with_size(4*scale);

        let mut ceremony = Groth16Ceremony::create_from_origin(&crs, &precom1);

        let trapdoor = Trapdoor::get_random_trapdoors();
        let contri = Groth16Ceremony::update_1(&ceremony, &mut crs, &trapdoor, &precom1);
    
        let start = Instant::now();
        let output = Groth16Ceremony::verify_1(&mut ceremony, &contri, &precom1);
        let duration = start.elapsed();
        println!("Phase 1: Verification result: {}. Time cost: {:#?}", output, duration);
        //Start the second phase
        ceremony.simple_change_phase(&mut crs, &precom2);
        let contri2 = Groth16Ceremony::update_2(&mut ceremony, &mut crs, &trapdoor, &precom2);
    
        let start = Instant::now();
        let output2 = Groth16Ceremony::verify_2(&mut ceremony, &contri2, &precom2);
        let duration = start.elapsed();
        println!("Phase 2: Verification result: {}. Time cost: {:#?}", output2, duration);
    }
}