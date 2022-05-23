#[macro_use]
extern crate criterion;
use rand::random;
use criterion::Criterion;
use num_modular::{FixedMersenneInt, ModularCoreOps, ModularPow, ModularUnaryOps};
use num_traits::{Inv, Pow};

pub fn bench_u128(c: &mut Criterion) {
    const N: usize = 256;
    let mut cases: [(u128, u128, u128); N] = [(0, 0, 0); N];
    for i in 0..N {
        cases[i] = (random(), random(), random());
    }

    let mut group = c.benchmark_group("u128 modular ops");
    group.bench_function("addm", |b| {
        b.iter(|| {
            cases.iter()
                .map(|&(a, b, m)| a.addm(b, &m))
                .reduce(|a, b| a.wrapping_add(b))
        })
    });
    group.bench_function("mulm", |b| {
        b.iter(|| {
            cases.iter()
                .map(|&(a, b, m)| a.mulm(b, &m))
                .reduce(|a, b| a.wrapping_add(b))
        })
    });
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
                .map(|n| FixedMersenneInt::<56, 5>::new(n as u128, &()).inv())
                .reduce(|a, b| a + b)
        })
    });
    group.bench_function("mersenne + fermat theorem", |b| {
        b.iter(|| {
            (100u64..400u64)
                .map(|n| FixedMersenneInt::<56, 5>::new(n as u128, &()).pow(M1 as u128 - 2))
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
                .map(|n| FixedMersenneInt::<94, 3>::new(n, &()).inv())
                .reduce(|a, b| a + b)
        })
    });
    group.bench_function("mersenne + fermat theorem", |b| {
        b.iter(|| {
            (1_000_000_000u128..1_000_000_300u128)
                .map(|n| FixedMersenneInt::<94, 3>::new(n, &()).pow(M2 - 2))
                .reduce(|a, b| a + b)
        })
    });

    group.finish();
}

criterion_group!(benches, bench_modinv, bench_u128);
criterion_main!(benches);
