#[macro_use]
extern crate criterion;
use criterion::Criterion;
use num_modular::{
    FixedMersenne64, FixedMersenneInt, FixedTrinomialSolinas64, ModularCoreOps,
    ModularPow, ModularUnaryOps, Montgomery, PreMulInv2by1, Reducer,
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

pub fn bench_modinv(c: &mut Criterion) {
    const M1: u64 = (1 << 56) - 5;
    let mut group = c.benchmark_group("modular inverse (small operands)");

    group.bench_function("extended gcd", |b| {
        b.iter(|| {
            (100u64..400u64)
                .map(|n| n.invm(&M1).unwrap())
                .reduce(|a, b| a.addm(b, &M1))
        })
    });
    group.bench_function("fermat theorem", |b| {
        b.iter(|| {
            (100u64..400u64)
                .map(|n| n.powm(M1 - 2, &M1))
                .reduce(|a, b| a.addm(b, &M1))
        })
    });
    group.bench_function("mersenne + extended gcd", |b| {
        b.iter(|| {
            (100u64..400u64)
                .map(|n| {
                    FixedMersenneInt::<56, 5>::new(n as u128, &(M1 as u128))
                        .inv()
                        .unwrap()
                })
                .reduce(|a, b| a + b)
        })
    });
    group.bench_function("mersenne + fermat theorem", |b| {
        b.iter(|| {
            (100u64..400u64)
                .map(|n| {
                    FixedMersenneInt::<56, 5>::new(n as u128, &(M1 as u128)).pow(&(M1 as u128 - 2))
                })
                .reduce(|a, b| a + b)
        })
    });

    group.finish();

    const M2: u128 = (1 << 94) - 3;
    let mut group = c.benchmark_group("modular inverse (large operands)");

    group.bench_function("extended gcd", |b| {
        b.iter(|| {
            (1_000_000_000u128..1_000_000_300u128)
                .map(|n| n.invm(&M2).unwrap())
                .reduce(|a, b| a.addm(b, &M2))
        })
    });
    group.bench_function("fermat theorem", |b| {
        b.iter(|| {
            (1_000_000_000u128..1_000_000_300u128)
                .map(|n| n.powm(M2 - 2, &M2))
                .reduce(|a, b| a.addm(b, &M2))
        })
    });
    group.bench_function("mersenne + extended gcd", |b| {
        b.iter(|| {
            (1_000_000_000u128..1_000_000_300u128)
                .map(|n| {
                    FixedMersenneInt::<94, 3>::new(n, &(M2 as u128))
                        .inv()
                        .unwrap()
                })
                .reduce(|a, b| a + b)
        })
    });
    group.bench_function("mersenne + fermat theorem", |b| {
        b.iter(|| {
            (1_000_000_000u128..1_000_000_300u128)
                .map(|n| FixedMersenneInt::<94, 3>::new(n, &(M2 as u128)).pow(&(M2 - 2)))
                .reduce(|a, b| a + b)
        })
    });

    group.finish();
}

pub fn bench_mod_pow(c: &mut Criterion) {
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

criterion_group!(benches, bench_modinv, bench_u128, bench_mod_pow);
criterion_main!(benches);
