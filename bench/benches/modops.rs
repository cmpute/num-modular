#[macro_use]
extern crate criterion;
use criterion::Criterion;
use num_traits::{Pow, Inv};
use num_modular::{ModularCoreOps, ModularOps, MersenneInt};

pub fn bench_modinv(c: &mut Criterion) {
    const M1: u64 = (1 << 56) - 5;
    let mut group = c.benchmark_group("modular inverse (small operands)");

    group.bench_function("extended gcd", |b| {
        b.iter(|| {
            (100u64..400u64)
                .map(|n| ModularOps::<u64>::invm(&n, &M1).unwrap())
                .reduce(|a, b| a.addm(b, &M1))
        })
    });
    group.bench_function("fermat theorem", |b| {
        b.iter(|| {
            (100u64..400u64)
                .map(|n| n.powm(M1-2, &M1))
                .reduce(|a, b| a.addm(b, &M1))
        })
    });
    group.bench_function("mersenne + extended gcd", |b| {
        b.iter(|| {
            (100u64..400u64)
                .map(|n| MersenneInt::<56, 5>::new(n as u128).inv())
                .reduce(|a, b| a + b)
        })
    });
    group.bench_function("mersenne + fermat theorem", |b| {
        b.iter(|| {
            (100u64..400u64)
                .map(|n| MersenneInt::<56, 5>::new(n as u128).pow(M1 as u128 - 2))
                .reduce(|a, b| a + b)
        })
    });

    group.finish();

    const M2: u128 = (1 << 94) - 3;
    let mut group = c.benchmark_group("modular inverse (large operands)");

    group.bench_function("extended gcd", |b| {
        b.iter(|| {
            (1_000_000_000u128..1_000_000_300u128)
                .map(|n| ModularOps::<u128>::invm(&n, &M2).unwrap())
                .reduce(|a, b| a.addm(b, &M2))
        })
    });
    group.bench_function("fermat theorem", |b| {
        b.iter(|| {
            (1_000_000_000u128..1_000_000_300u128)
                .map(|n| n.powm(M2-2, &M2))
                .reduce(|a, b| a.addm(b, &M2))
        })
    });
    group.bench_function("mersenne +  + extended gcd", |b| {
        b.iter(|| {
            (1_000_000_000u128..1_000_000_300u128)
                .map(|n| MersenneInt::<94, 3>::new(n).inv())
                .reduce(|a, b| a + b)
        })
    });
    group.bench_function("mersenne + fermat theorem", |b| {
        b.iter(|| {
            (1_000_000_000u128..1_000_000_300u128)
                .map(|n| MersenneInt::<94, 3>::new(n).pow(M2-2))
                .reduce(|a, b| a + b)
        })
    });

    group.finish();
}

criterion_group!(benches, bench_modinv);
criterion_main!(benches);
