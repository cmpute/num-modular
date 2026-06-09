#[macro_use]
extern crate criterion;
use criterion::Criterion;
use num_modular::{FixedTrinomialSolinas, FixedTrinomialSolinas32, FixedTrinomialSolinas64, Reducer};
use rand::random;

/// Solinas-like prime 2^6 - 2^2 + 1 = 61, fits in u32, u64, and u128.
const P1: u8 = 6;
const P2: u8 = 2;
const K32: i32 = 1;
const K64: i64 = 1;
const K128: i128 = 1;

const MOD32: u32 = (1u32 << P1) - (1u32 << P2) + (K32 as u32);
const MOD64: u64 = (1u64 << P1) - (1u64 << P2) + (K64 as u64);
const MOD128: u128 = (1u128 << P1) - (1u128 << P2) + (K128 as u128);

const N: usize = 256;

pub fn bench_transform(c: &mut Criterion) {
    let mut inputs: [u32; N] = [0; N];
    for i in 0..N {
        inputs[i] = random::<u32>() % MOD32;
    }

    let reducer32 = FixedTrinomialSolinas32::<P1, P2, K32>::new(&MOD32);
    let reducer64 = FixedTrinomialSolinas64::<P1, P2, K64>::new(&MOD64);
    let reducer128 = FixedTrinomialSolinas::<P1, P2, K128>::new(&MOD128);

    let mut group = c.benchmark_group("solinas transform (2^6 - 2^2 + 1)");
    group.bench_function("FixedTrinomialSolinas32", |b| {
        b.iter(|| {
            inputs
                .iter()
                .map(|&v| reducer32.transform(v))
                .reduce(|a, b| a.wrapping_add(b))
        })
    });
    group.bench_function("FixedTrinomialSolinas64", |b| {
        b.iter(|| {
            inputs
                .iter()
                .map(|&v| reducer64.transform(v as u64))
                .reduce(|a, b| a.wrapping_add(b))
        })
    });
    group.bench_function("FixedTrinomialSolinas", |b| {
        b.iter(|| {
            inputs
                .iter()
                .map(|&v| reducer128.transform(v as u128))
                .reduce(|a, b| a.wrapping_add(b))
        })
    });
    group.finish();
}

pub fn bench_mul(c: &mut Criterion) {
    let mut lhs: [u32; N] = [0; N];
    let mut rhs: [u32; N] = [0; N];
    for i in 0..N {
        lhs[i] = random::<u32>() % MOD32;
        rhs[i] = random::<u32>() % MOD32;
    }

    let reducer32 = FixedTrinomialSolinas32::<P1, P2, K32>::new(&MOD32);
    let reducer64 = FixedTrinomialSolinas64::<P1, P2, K64>::new(&MOD64);
    let reducer128 = FixedTrinomialSolinas::<P1, P2, K128>::new(&MOD128);

    let lhs32: Vec<u32> = lhs.iter().map(|&v| reducer32.transform(v)).collect();
    let rhs32: Vec<u32> = rhs.iter().map(|&v| reducer32.transform(v)).collect();
    let lhs64: Vec<u64> = lhs.iter().map(|&v| reducer64.transform(v as u64)).collect();
    let rhs64: Vec<u64> = rhs.iter().map(|&v| reducer64.transform(v as u64)).collect();
    let lhs128: Vec<u128> = lhs.iter().map(|&v| reducer128.transform(v as u128)).collect();
    let rhs128: Vec<u128> = rhs.iter().map(|&v| reducer128.transform(v as u128)).collect();

    let mut group = c.benchmark_group("solinas mul (2^6 - 2^2 + 1)");
    group.bench_function("FixedTrinomialSolinas32", |b| {
        b.iter(|| {
            lhs32
                .iter()
                .zip(rhs32.iter())
                .map(|(a, b)| reducer32.mul(a, b))
                .reduce(|a, b| a.wrapping_add(b))
        })
    });
    group.bench_function("FixedTrinomialSolinas64", |b| {
        b.iter(|| {
            lhs64
                .iter()
                .zip(rhs64.iter())
                .map(|(a, b)| reducer64.mul(a, b))
                .reduce(|a, b| a.wrapping_add(b))
        })
    });
    group.bench_function("FixedTrinomialSolinas", |b| {
        b.iter(|| {
            lhs128
                .iter()
                .zip(rhs128.iter())
                .map(|(a, b)| reducer128.mul(a, b))
                .reduce(|a, b| a.wrapping_add(b))
        })
    });
    group.finish();
}

pub fn bench_sqr(c: &mut Criterion) {
    let mut inputs: [u32; N] = [0; N];
    for i in 0..N {
        inputs[i] = random::<u32>() % MOD32;
    }

    let reducer32 = FixedTrinomialSolinas32::<P1, P2, K32>::new(&MOD32);
    let reducer64 = FixedTrinomialSolinas64::<P1, P2, K64>::new(&MOD64);
    let reducer128 = FixedTrinomialSolinas::<P1, P2, K128>::new(&MOD128);

    let sqr32: Vec<u32> = inputs.iter().map(|&v| reducer32.transform(v)).collect();
    let sqr64: Vec<u64> = inputs.iter().map(|&v| reducer64.transform(v as u64)).collect();
    let sqr128: Vec<u128> = inputs.iter().map(|&v| reducer128.transform(v as u128)).collect();

    let mut group = c.benchmark_group("solinas sqr (2^6 - 2^2 + 1)");
    group.bench_function("FixedTrinomialSolinas32", |b| {
        b.iter(|| {
            sqr32
                .iter()
                .map(|&v| reducer32.sqr(v))
                .reduce(|a, b| a.wrapping_add(b))
        })
    });
    group.bench_function("FixedTrinomialSolinas64", |b| {
        b.iter(|| {
            sqr64
                .iter()
                .map(|&v| reducer64.sqr(v))
                .reduce(|a, b| a.wrapping_add(b))
        })
    });
    group.bench_function("FixedTrinomialSolinas", |b| {
        b.iter(|| {
            sqr128
                .iter()
                .map(|&v| reducer128.sqr(v))
                .reduce(|a, b| a.wrapping_add(b))
        })
    });
    group.finish();
}

pub fn bench_add(c: &mut Criterion) {
    let mut lhs: [u32; N] = [0; N];
    let mut rhs: [u32; N] = [0; N];
    for i in 0..N {
        lhs[i] = random::<u32>() % MOD32;
        rhs[i] = random::<u32>() % MOD32;
    }

    let reducer32 = FixedTrinomialSolinas32::<P1, P2, K32>::new(&MOD32);
    let reducer64 = FixedTrinomialSolinas64::<P1, P2, K64>::new(&MOD64);
    let reducer128 = FixedTrinomialSolinas::<P1, P2, K128>::new(&MOD128);

    let lhs32: Vec<u32> = lhs.iter().map(|&v| reducer32.transform(v)).collect();
    let rhs32: Vec<u32> = rhs.iter().map(|&v| reducer32.transform(v)).collect();
    let lhs64: Vec<u64> = lhs.iter().map(|&v| reducer64.transform(v as u64)).collect();
    let rhs64: Vec<u64> = rhs.iter().map(|&v| reducer64.transform(v as u64)).collect();
    let lhs128: Vec<u128> = lhs.iter().map(|&v| reducer128.transform(v as u128)).collect();
    let rhs128: Vec<u128> = rhs.iter().map(|&v| reducer128.transform(v as u128)).collect();

    let mut group = c.benchmark_group("solinas add (2^6 - 2^2 + 1)");
    group.bench_function("FixedTrinomialSolinas32", |b| {
        b.iter(|| {
            lhs32
                .iter()
                .zip(rhs32.iter())
                .map(|(a, b)| reducer32.add(a, b))
                .reduce(|a, b| a.wrapping_add(b))
        })
    });
    group.bench_function("FixedTrinomialSolinas64", |b| {
        b.iter(|| {
            lhs64
                .iter()
                .zip(rhs64.iter())
                .map(|(a, b)| reducer64.add(a, b))
                .reduce(|a, b| a.wrapping_add(b))
        })
    });
    group.bench_function("FixedTrinomialSolinas", |b| {
        b.iter(|| {
            lhs128
                .iter()
                .zip(rhs128.iter())
                .map(|(a, b)| reducer128.add(a, b))
                .reduce(|a, b| a.wrapping_add(b))
        })
    });
    group.finish();
}

pub fn bench_inv(c: &mut Criterion) {
    // Generate values coprime to 61 (guaranteed since modulus is prime)
    let mut inputs: [u32; N] = [0; N];
    for i in 0..N {
        let v = random::<u32>() % MOD32;
        inputs[i] = if v == 0 { 1 } else { v };
    }

    let reducer32 = FixedTrinomialSolinas32::<P1, P2, K32>::new(&MOD32);
    let reducer64 = FixedTrinomialSolinas64::<P1, P2, K64>::new(&MOD64);
    let reducer128 = FixedTrinomialSolinas::<P1, P2, K128>::new(&MOD128);

    let inv32: Vec<u32> = inputs.iter().map(|&v| reducer32.transform(v)).collect();
    let inv64: Vec<u64> = inputs.iter().map(|&v| reducer64.transform(v as u64)).collect();
    let inv128: Vec<u128> = inputs.iter().map(|&v| reducer128.transform(v as u128)).collect();

    let mut group = c.benchmark_group("solinas inv (2^6 - 2^2 + 1)");
    group.bench_function("FixedTrinomialSolinas32", |b| {
        b.iter(|| {
            inv32
                .iter()
                .map(|&v| reducer32.inv(v))
                .reduce(|a, b| Some(a?.wrapping_add(b?)))
        })
    });
    group.bench_function("FixedTrinomialSolinas64", |b| {
        b.iter(|| {
            inv64
                .iter()
                .map(|&v| reducer64.inv(v))
                .reduce(|a, b| Some(a?.wrapping_add(b?)))
        })
    });
    group.bench_function("FixedTrinomialSolinas", |b| {
        b.iter(|| {
            inv128
                .iter()
                .map(|&v| reducer128.inv(v))
                .reduce(|a, b| Some(a?.wrapping_add(b?)))
        })
    });
    group.finish();
}

criterion_group!(benches, bench_transform, bench_add, bench_mul, bench_sqr, bench_inv);
criterion_main!(benches);
