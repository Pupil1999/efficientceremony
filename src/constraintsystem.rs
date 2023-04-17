use ark_bls12_381::Fr;
pub struct Univariate{
    pub coeffs: Vec<Fr>,
    pub degree: u64
}

impl Univariate {
    pub fn create_from_u64(coeffs: &Vec<u64>) -> Self {
        let mut coe_fr: Vec<Fr> = Vec::new();
        let mut degree = 0;
        for coe_64 in coeffs {
            coe_fr.push(Fr::from(*coe_64));
            degree += 1;
        }
        
        Self{ coeffs: coe_fr, degree: degree - 1 }
    }
}

pub struct ConstraintSystem {
    pub m: usize,
    pub n: usize,
    pub l: usize,
    pub u: Vec<Univariate>,
    pub v: Vec<Univariate>,
    pub w: Vec<Univariate>,
    pub t: Univariate,
}

impl ConstraintSystem {
    pub fn create_simple(m: usize, n: usize, l: usize) -> Self {
        let mut u = Vec::<Univariate>::with_capacity(m + 1);
        let mut v = Vec::<Univariate>::with_capacity(m + 1);
        let mut w = Vec::<Univariate>::with_capacity(m + 1);

        let gate_polys = vec![1; n - 1];

        for _i in 0..m - l {
            u.push(Univariate::create_from_u64(&gate_polys));
            v.push(Univariate::create_from_u64(&gate_polys));
            w.push(Univariate::create_from_u64(&gate_polys));
        }

        let vanishing_poly = vec![1; n];
        let t = Univariate::create_from_u64(&vanishing_poly);

        // Universal part of the CRS will be padding to be exponential size
        let actual_n = 2_usize.pow(ark_std::log2(n));

        Self {
            m,
            n: actual_n,
            l,
            u,
            v,
            w,
            t,
        }
    }

    pub fn get_size(&self) -> (usize, usize, usize) {
        let m = self.m;
        let n = self.n;
        let l = self.l;
        (m, n, l)
    }
}