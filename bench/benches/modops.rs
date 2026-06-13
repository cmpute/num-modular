#[macro_use]
extern crate criterion;
use criterion::Criterion;
use num_modular::{
    FixedMersenne64, FixedProth64, FixedTrinomialSolinas64, ModularCoreOps, Montgomery,
    PreMulInv2by1, Reducer,
};
use rand::random;

pub fn bench_u128(c: &mut Criterion) {
    const N: usize = 256;
    let mut cases: [(u128, u128, u128); N] = [(0, 0, 0); N];
    for i in 0..N {
        let (a, b, mut m) = (random(), random(), random());
        if m == 0 {
            m = 1;
        }
        cases[i] = (a, b, m);
    }

    let mut group = c.benchmark_group("u128 modular ops");
    group.bench_function("addm", |b| {
        b.iter(|| {
            cases
                .iter()
                .map(|&(a, b, m)| a.addm(b, &m))
                .reduce(|a, b| a.wrapping_add(b))
        })
    });
    group.bench_function("mulm", |b| {
        b.iter(|| {
            cases
                .iter()
                .map(|&(a, b, m)| a.mulm(b, &m))
                .reduce(|a, b| a.wrapping_add(b))
        })
    });

    group.finish();
}

pub fn bench_mod_inv_goldilocks(c: &mut Criterion) {
    const N: usize = 256;

    // Goldilocks prime: 2^64 - 2^32 + 1
    type Sol = FixedTrinomialSolinas64<64, 32, 1>;
    type Mer = FixedMersenne64<64, { (1u64 << 32) - 1 }>;
    const MOD: u64 = <Sol>::MODULUS;

    let mer = Mer::new(&MOD);
    let sol = Sol::new(&MOD);
    let monty = Montgomery::<u64>::new(MOD);
    let premul = PreMulInv2by1::<u64>::new(MOD);

    let mut bases_mer = [0u64; N];
    let mut bases_sol = [0u64; N];
    for i in 0..N {
        let v = random::<u64>() % MOD;
        bases_mer[i] = mer.transform(if v == 0 { 1 } else { v });
        bases_sol[i] = sol.transform(if v == 0 { 1 } else { v });
    }

    let mut group = c.benchmark_group("mod_inv (2^64 - 2^32 + 1)");
    group.bench_function("Mersenne", |b| {
        b.iter(|| {
            bases_mer.iter()
                .map(|&v| mer.inv(v))
                .reduce(|a, b| Some(a?.wrapping_add(b?)))
        })
    });
    group.bench_function("Solinas", |b| {
        b.iter(|| {
            bases_sol.iter()
                .map(|&v| sol.inv(v))
                .reduce(|a, b| Some(a?.wrapping_add(b?)))
        })
    });
    group.bench_function("Montgomery", |b| {
        b.iter(|| {
            bases_sol.iter()
                .map(|&v| monty.inv(v))
                .reduce(|a, b| Some(a?.wrapping_add(b?)))
        })
    });
    group.bench_function("PreMulInv2by1", |b| {
        b.iter(|| {
            bases_sol.iter()
                .map(|&v| premul.inv(v))
                .reduce(|a, b| Some(a?.wrapping_add(b?)))
        })
    });
    group.finish();
}

pub fn bench_mod_inv_proth(c: &mut Criterion) {
    const N: usize = 256;

    type Pro = FixedProth64<32, 3>;
    const MOD: u64 = <Pro>::MODULUS;

    let pro = Pro::new(&MOD);
    let monty = Montgomery::<u64>::new(MOD);
    let premul = PreMulInv2by1::<u64>::new(MOD);

    let mut bases = [0u64; N];
    for i in 0..N {
        let v = random::<u64>() % MOD;
        bases[i] = pro.transform(if v == 0 { 1 } else { v });
    }

    let mut group = c.benchmark_group("mod_inv proth 3*2^32+1");
    group.bench_function("FixedProth64", |b| {
        b.iter(|| {
            bases.iter()
                .map(|&v| pro.inv(v))
                .reduce(|a, b| Some(a?.wrapping_add(b?)))
        })
    });
    group.bench_function("Montgomery", |b| {
        b.iter(|| {
            bases.iter()
                .map(|&v| monty.inv(v))
                .reduce(|a, b| Some(a?.wrapping_add(b?)))
        })
    });
    group.bench_function("PreMulInv2by1", |b| {
        b.iter(|| {
            bases.iter()
                .map(|&v| premul.inv(v))
                .reduce(|a, b| Some(a?.wrapping_add(b?)))
        })
    });
    group.finish();
}

pub fn bench_mod_pow_goldilocks(c: &mut Criterion) {
    const N: usize = 256;

    // Goldilocks prime: 2^64 - 2^32 + 1
    // As Mersenne: 2^64 - K where K = 2^32 - 1
    type Sol = FixedTrinomialSolinas64<64, 32, 1>;
    type Mer = FixedMersenne64<64, { (1u64 << 32) - 1 }>;
    const MOD: u64 = <Sol>::MODULUS;

    let mer = Mer::new(&MOD);
    let sol = Sol::new(&MOD);
    let monty = Montgomery::<u64>::new(MOD);
    let premul = PreMulInv2by1::<u64>::new(MOD);

    // Generate random bases, pre-transform into each reducer's form
    let mut bases_mer = [0u64; N];
    let mut bases_sol = [0u64; N];
    for i in 0..N {
        bases_mer[i] = mer.transform(random::<u64>() % MOD);
        bases_sol[i] = sol.transform(random::<u64>() % MOD);
    }

    let exp = MOD - 2;

    let mut group = c.benchmark_group("mod_pow");
    group.bench_function("Mersenne (2^64 - 2^32 + 1)", |b| {
        b.iter(|| {
            bases_mer
                .iter()
                .map(|&v| mer.pow(v, &exp))
                .reduce(|a, b| a.wrapping_add(b))
        })
    });
    group.bench_function("Solinas (2^64 - 2^32 + 1)", |b| {
        b.iter(|| {
            bases_sol
                .iter()
                .map(|&v| sol.pow(v, &exp))
                .reduce(|a, b| a.wrapping_add(b))
        })
    });
    group.bench_function("Montgomery (2^64 - 2^32 + 1)", |b| {
        b.iter(|| {
            bases_sol
                .iter()
                .map(|&v| monty.pow(v, &exp))
                .reduce(|a, b| a.wrapping_add(b))
        })
    });
    group.bench_function("PreMulInv2by1 (2^64 - 2^32 + 1)", |b| {
        b.iter(|| {
            bases_sol
                .iter()
                .map(|&v| premul.pow(v, &exp))
                .reduce(|a, b| a.wrapping_add(b))
        })
    });
    group.finish();
}

pub fn bench_mod_pow_proth(c: &mut Criterion) {
    const N: usize = 256;

    // Proth prime: 3 * 2^32 + 1
    type Pro = FixedProth64<32, 3>;
    const MOD: u64 = <Pro>::MODULUS;

    let pro = Pro::new(&MOD);
    let monty = Montgomery::<u64>::new(MOD);
    let premul = PreMulInv2by1::<u64>::new(MOD);

    // Generate random bases, pre-transform into each reducer's form
    let mut bases_pro = [0u64; N];
    for i in 0..N {
        bases_pro[i] = pro.transform(random::<u64>() % MOD);
    }

    let exp = MOD - 2;

    let mut group = c.benchmark_group("mod_pow proth 3*2^32+1");
    group.bench_function("FixedProth64", |b| {
        b.iter(|| {
            bases_pro
                .iter()
                .map(|&v| pro.pow(v, &exp))
                .reduce(|a, b| a.wrapping_add(b))
        })
    });
    group.bench_function("Montgomery", |b| {
        b.iter(|| {
            bases_pro
                .iter()
                .map(|&v| monty.pow(v, &exp))
                .reduce(|a, b| a.wrapping_add(b))
        })
    });
    group.bench_function("PreMulInv2by1", |b| {
        b.iter(|| {
            bases_pro
                .iter()
                .map(|&v| premul.pow(v, &exp))
                .reduce(|a, b| a.wrapping_add(b))
        })
    });
    group.finish();
}

pub fn bench_mod_add_goldilocks(c: &mut Criterion) {
    const N: usize = 256;

    // Goldilocks prime: 2^64 - 2^32 + 1
    type Sol = FixedTrinomialSolinas64<64, 32, 1>;
    type Mer = FixedMersenne64<64, { (1u64 << 32) - 1 }>;
    const MOD: u64 = <Sol>::MODULUS;

    let mer = Mer::new(&MOD);
    let sol = Sol::new(&MOD);
    let monty = Montgomery::<u64>::new(MOD);
    let premul = PreMulInv2by1::<u64>::new(MOD);

    let mut lhs_mer = [0u64; N];
    let mut rhs_mer = [0u64; N];
    let mut lhs_sol = [0u64; N];
    let mut rhs_sol = [0u64; N];
    for i in 0..N {
        lhs_mer[i] = mer.transform(random::<u64>() % MOD);
        rhs_mer[i] = mer.transform(random::<u64>() % MOD);
        lhs_sol[i] = sol.transform(random::<u64>() % MOD);
        rhs_sol[i] = sol.transform(random::<u64>() % MOD);
    }

    let mut group = c.benchmark_group("mod_add (2^64 - 2^32 + 1)");
    group.bench_function("Mersenne", |b| {
        b.iter(|| {
            lhs_mer.iter().zip(rhs_mer.iter())
                .map(|(a, b)| mer.add(a, b))
                .reduce(|a, b| a.wrapping_add(b))
        })
    });
    group.bench_function("Solinas", |b| {
        b.iter(|| {
            lhs_sol.iter().zip(rhs_sol.iter())
                .map(|(a, b)| sol.add(a, b))
                .reduce(|a, b| a.wrapping_add(b))
        })
    });
    group.bench_function("Montgomery", |b| {
        b.iter(|| {
            lhs_sol.iter().zip(rhs_sol.iter())
                .map(|(a, b)| monty.add(a, b))
                .reduce(|a, b| a.wrapping_add(b))
        })
    });
    group.bench_function("PreMulInv2by1", |b| {
        b.iter(|| {
            lhs_sol.iter().zip(rhs_sol.iter())
                .map(|(a, b)| premul.add(a, b))
                .reduce(|a, b| a.wrapping_add(b))
        })
    });
    group.finish();
}

pub fn bench_mod_add_proth(c: &mut Criterion) {
    const N: usize = 256;

    type Pro = FixedProth64<32, 3>;
    const MOD: u64 = <Pro>::MODULUS;

    let pro = Pro::new(&MOD);
    let monty = Montgomery::<u64>::new(MOD);
    let premul = PreMulInv2by1::<u64>::new(MOD);

    let mut lhs = [0u64; N];
    let mut rhs = [0u64; N];
    for i in 0..N {
        lhs[i] = pro.transform(random::<u64>() % MOD);
        rhs[i] = pro.transform(random::<u64>() % MOD);
    }

    let mut group = c.benchmark_group("mod_add proth 3*2^32+1");
    group.bench_function("FixedProth64", |b| {
        b.iter(|| {
            lhs.iter().zip(rhs.iter())
                .map(|(a, b)| pro.add(a, b))
                .reduce(|a, b| a.wrapping_add(b))
        })
    });
    group.bench_function("Montgomery", |b| {
        b.iter(|| {
            lhs.iter().zip(rhs.iter())
                .map(|(a, b)| monty.add(a, b))
                .reduce(|a, b| a.wrapping_add(b))
        })
    });
    group.bench_function("PreMulInv2by1", |b| {
        b.iter(|| {
            lhs.iter().zip(rhs.iter())
                .map(|(a, b)| premul.add(a, b))
                .reduce(|a, b| a.wrapping_add(b))
        })
    });
    group.finish();
}

pub fn bench_mod_mul_goldilocks(c: &mut Criterion) {
    const N: usize = 256;

    type Sol = FixedTrinomialSolinas64<64, 32, 1>;
    type Mer = FixedMersenne64<64, { (1u64 << 32) - 1 }>;
    const MOD: u64 = <Sol>::MODULUS;

    let mer = Mer::new(&MOD);
    let sol = Sol::new(&MOD);
    let monty = Montgomery::<u64>::new(MOD);
    let premul = PreMulInv2by1::<u64>::new(MOD);

    let mut lhs_mer = [0u64; N];
    let mut rhs_mer = [0u64; N];
    let mut lhs_sol = [0u64; N];
    let mut rhs_sol = [0u64; N];
    for i in 0..N {
        lhs_mer[i] = mer.transform(random::<u64>() % MOD);
        rhs_mer[i] = mer.transform(random::<u64>() % MOD);
        lhs_sol[i] = sol.transform(random::<u64>() % MOD);
        rhs_sol[i] = sol.transform(random::<u64>() % MOD);
    }

    let mut group = c.benchmark_group("mod_mul (2^64 - 2^32 + 1)");
    group.bench_function("Mersenne", |b| {
        b.iter(|| {
            lhs_mer.iter().zip(rhs_mer.iter())
                .map(|(a, b)| mer.mul(a, b))
                .reduce(|a, b| a.wrapping_add(b))
        })
    });
    group.bench_function("Solinas", |b| {
        b.iter(|| {
            lhs_sol.iter().zip(rhs_sol.iter())
                .map(|(a, b)| sol.mul(a, b))
                .reduce(|a, b| a.wrapping_add(b))
        })
    });
    group.bench_function("Montgomery", |b| {
        b.iter(|| {
            lhs_sol.iter().zip(rhs_sol.iter())
                .map(|(a, b)| monty.mul(a, b))
                .reduce(|a, b| a.wrapping_add(b))
        })
    });
    group.bench_function("PreMulInv2by1", |b| {
        b.iter(|| {
            lhs_sol.iter().zip(rhs_sol.iter())
                .map(|(a, b)| premul.mul(a, b))
                .reduce(|a, b| a.wrapping_add(b))
        })
    });
    group.finish();
}

pub fn bench_mod_mul_proth(c: &mut Criterion) {
    const N: usize = 256;

    type Pro = FixedProth64<32, 3>;
    const MOD: u64 = <Pro>::MODULUS;

    let pro = Pro::new(&MOD);
    let monty = Montgomery::<u64>::new(MOD);
    let premul = PreMulInv2by1::<u64>::new(MOD);

    let mut lhs = [0u64; N];
    let mut rhs = [0u64; N];
    for i in 0..N {
        lhs[i] = pro.transform(random::<u64>() % MOD);
        rhs[i] = pro.transform(random::<u64>() % MOD);
    }

    let mut group = c.benchmark_group("mod_mul proth 3*2^32+1");
    group.bench_function("FixedProth64", |b| {
        b.iter(|| {
            lhs.iter().zip(rhs.iter())
                .map(|(a, b)| pro.mul(a, b))
                .reduce(|a, b| a.wrapping_add(b))
        })
    });
    group.bench_function("Montgomery", |b| {
        b.iter(|| {
            lhs.iter().zip(rhs.iter())
                .map(|(a, b)| monty.mul(a, b))
                .reduce(|a, b| a.wrapping_add(b))
        })
    });
    group.bench_function("PreMulInv2by1", |b| {
        b.iter(|| {
            lhs.iter().zip(rhs.iter())
                .map(|(a, b)| premul.mul(a, b))
                .reduce(|a, b| a.wrapping_add(b))
        })
    });
    group.finish();
}

criterion_group!(
    benches,
    bench_mod_inv_goldilocks,
    bench_mod_inv_proth,
    bench_u128,
    bench_mod_add_goldilocks,
    bench_mod_add_proth,
    bench_mod_mul_goldilocks,
    bench_mod_mul_proth,
    bench_mod_pow_goldilocks,
    bench_mod_pow_proth
);
criterion_main!(benches);
