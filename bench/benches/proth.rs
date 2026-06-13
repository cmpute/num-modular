#[macro_use]
extern crate criterion;
use criterion::Criterion;
use num_modular::{FixedProth32, FixedProth64, Reducer};
use rand::random;

const N_ENTRIES: usize = 256;

// N=8, K=1 → m=257 — small modulus, common path
const N_SMALL: u8 = 8;
const K_SMALL: u8 = 1;

const MOD32_SMALL: u32 = (K_SMALL as u32) * (1u32 << N_SMALL) + 1;
const MOD64_SMALL: u64 = (K_SMALL as u64) * (1u64 << N_SMALL) + 1;

// N=31, K=1 → m=2^31+1 — near u32 width
const N_WIDE32: u8 = 31;

const MOD32_WIDE: u32 = (K_SMALL as u32) * (1u32 << N_WIDE32) + 1;

// N=63, K=1 → m=2^63+1 — near u64 width
const N_WIDE64: u8 = 63;

const MOD64_WIDE: u64 = (K_SMALL as u64) * (1u64 << N_WIDE64) + 1;

pub fn bench_transform(c: &mut Criterion) {
    let mut inputs: [u32; N_ENTRIES] = [0; N_ENTRIES];
    for i in 0..N_ENTRIES {
        inputs[i] = random::<u32>() % MOD32_SMALL;
    }

    let reducer32 = FixedProth32::<N_SMALL, K_SMALL>::new(&MOD32_SMALL);
    let reducer64 = FixedProth64::<N_SMALL, K_SMALL>::new(&MOD64_SMALL);

    let mut group = c.benchmark_group("proth transform (2^8 + 1)");
    group.bench_function("FixedProth32", |b| {
        b.iter(|| {
            inputs
                .iter()
                .map(|&v| reducer32.transform(v))
                .reduce(|a, b| a.wrapping_add(b))
        })
    });
    group.bench_function("FixedProth64", |b| {
        b.iter(|| {
            inputs
                .iter()
                .map(|&v| reducer64.transform(v as u64))
                .reduce(|a, b| a.wrapping_add(b))
        })
    });
    group.finish();
}

pub fn bench_mul(c: &mut Criterion) {
    let mut lhs: [u32; N_ENTRIES] = [0; N_ENTRIES];
    let mut rhs: [u32; N_ENTRIES] = [0; N_ENTRIES];
    for i in 0..N_ENTRIES {
        lhs[i] = random::<u32>() % MOD32_SMALL;
        rhs[i] = random::<u32>() % MOD32_SMALL;
    }

    let reducer32 = FixedProth32::<N_SMALL, K_SMALL>::new(&MOD32_SMALL);
    let reducer64 = FixedProth64::<N_SMALL, K_SMALL>::new(&MOD64_SMALL);

    let lhs32: Vec<u32> = lhs.iter().map(|&v| reducer32.transform(v)).collect();
    let rhs32: Vec<u32> = rhs.iter().map(|&v| reducer32.transform(v)).collect();
    let lhs64: Vec<u64> = lhs.iter().map(|&v| reducer64.transform(v as u64)).collect();
    let rhs64: Vec<u64> = rhs.iter().map(|&v| reducer64.transform(v as u64)).collect();

    let mut group = c.benchmark_group("proth mul (2^8 + 1)");
    group.bench_function("FixedProth32", |b| {
        b.iter(|| {
            lhs32
                .iter()
                .zip(rhs32.iter())
                .map(|(a, b)| reducer32.mul(a, b))
                .reduce(|a, b| a.wrapping_add(b))
        })
    });
    group.bench_function("FixedProth64", |b| {
        b.iter(|| {
            lhs64
                .iter()
                .zip(rhs64.iter())
                .map(|(a, b)| reducer64.mul(a, b))
                .reduce(|a, b| a.wrapping_add(b))
        })
    });
    group.finish();
}

pub fn bench_sqr(c: &mut Criterion) {
    let mut inputs: [u32; N_ENTRIES] = [0; N_ENTRIES];
    for i in 0..N_ENTRIES {
        inputs[i] = random::<u32>() % MOD32_SMALL;
    }

    let reducer32 = FixedProth32::<N_SMALL, K_SMALL>::new(&MOD32_SMALL);
    let reducer64 = FixedProth64::<N_SMALL, K_SMALL>::new(&MOD64_SMALL);

    let sqr32: Vec<u32> = inputs.iter().map(|&v| reducer32.transform(v)).collect();
    let sqr64: Vec<u64> = inputs
        .iter()
        .map(|&v| reducer64.transform(v as u64))
        .collect();

    let mut group = c.benchmark_group("proth sqr (2^8 + 1)");
    group.bench_function("FixedProth32", |b| {
        b.iter(|| {
            sqr32
                .iter()
                .map(|&v| reducer32.sqr(v))
                .reduce(|a, b| a.wrapping_add(b))
        })
    });
    group.bench_function("FixedProth64", |b| {
        b.iter(|| {
            sqr64
                .iter()
                .map(|&v| reducer64.sqr(v))
                .reduce(|a, b| a.wrapping_add(b))
        })
    });
    group.finish();
}

pub fn bench_add(c: &mut Criterion) {
    let mut lhs: [u32; N_ENTRIES] = [0; N_ENTRIES];
    let mut rhs: [u32; N_ENTRIES] = [0; N_ENTRIES];
    for i in 0..N_ENTRIES {
        lhs[i] = random::<u32>() % MOD32_SMALL;
        rhs[i] = random::<u32>() % MOD32_SMALL;
    }

    let reducer32 = FixedProth32::<N_SMALL, K_SMALL>::new(&MOD32_SMALL);
    let reducer64 = FixedProth64::<N_SMALL, K_SMALL>::new(&MOD64_SMALL);

    let lhs32: Vec<u32> = lhs.iter().map(|&v| reducer32.transform(v)).collect();
    let rhs32: Vec<u32> = rhs.iter().map(|&v| reducer32.transform(v)).collect();
    let lhs64: Vec<u64> = lhs.iter().map(|&v| reducer64.transform(v as u64)).collect();
    let rhs64: Vec<u64> = rhs.iter().map(|&v| reducer64.transform(v as u64)).collect();

    let mut group = c.benchmark_group("proth add (2^8 + 1)");
    group.bench_function("FixedProth32", |b| {
        b.iter(|| {
            lhs32
                .iter()
                .zip(rhs32.iter())
                .map(|(a, b)| reducer32.add(a, b))
                .reduce(|a, b| a.wrapping_add(b))
        })
    });
    group.bench_function("FixedProth64", |b| {
        b.iter(|| {
            lhs64
                .iter()
                .zip(rhs64.iter())
                .map(|(a, b)| reducer64.add(a, b))
                .reduce(|a, b| a.wrapping_add(b))
        })
    });
    group.finish();
}

pub fn bench_inv(c: &mut Criterion) {
    let mut inputs: [u32; N_ENTRIES] = [0; N_ENTRIES];
    for i in 0..N_ENTRIES {
        let v = random::<u32>() % MOD32_SMALL;
        inputs[i] = if v == 0 { 1 } else { v };
    }

    let reducer32 = FixedProth32::<N_SMALL, K_SMALL>::new(&MOD32_SMALL);
    let reducer64 = FixedProth64::<N_SMALL, K_SMALL>::new(&MOD64_SMALL);

    let inv32: Vec<u32> = inputs.iter().map(|&v| reducer32.transform(v)).collect();
    let inv64: Vec<u64> = inputs
        .iter()
        .map(|&v| reducer64.transform(v as u64))
        .collect();

    let mut group = c.benchmark_group("proth inv (2^8 + 1)");
    group.bench_function("FixedProth32", |b| {
        b.iter(|| {
            inv32
                .iter()
                .map(|&v| reducer32.inv(v))
                .reduce(|a, b| Some(a?.wrapping_add(b?)))
        })
    });
    group.bench_function("FixedProth64", |b| {
        b.iter(|| {
            inv64
                .iter()
                .map(|&v| reducer64.inv(v))
                .reduce(|a, b| Some(a?.wrapping_add(b?)))
        })
    });
    group.finish();
}

pub fn bench_mul_near_width_u32(c: &mut Criterion) {
    let mut lhs: [u32; N_ENTRIES] = [0; N_ENTRIES];
    let mut rhs: [u32; N_ENTRIES] = [0; N_ENTRIES];
    for i in 0..N_ENTRIES {
        lhs[i] = random::<u32>() % MOD32_WIDE;
        rhs[i] = random::<u32>() % MOD32_WIDE;
    }

    let reducer32 = FixedProth32::<N_WIDE32, K_SMALL>::new(&MOD32_WIDE);
    let lhs32: Vec<u32> = lhs.iter().map(|&v| reducer32.transform(v)).collect();
    let rhs32: Vec<u32> = rhs.iter().map(|&v| reducer32.transform(v)).collect();

    let mut group = c.benchmark_group("proth mul (2^31 + 1)");
    group.bench_function("FixedProth32", |b| {
        b.iter(|| {
            lhs32
                .iter()
                .zip(rhs32.iter())
                .map(|(a, b)| reducer32.mul(a, b))
                .reduce(|a, b| a.wrapping_add(b))
        })
    });
    group.finish();
}

pub fn bench_sqr_near_width_u32(c: &mut Criterion) {
    let mut inputs: [u32; N_ENTRIES] = [0; N_ENTRIES];
    for i in 0..N_ENTRIES {
        inputs[i] = random::<u32>() % MOD32_WIDE;
    }

    let reducer32 = FixedProth32::<N_WIDE32, K_SMALL>::new(&MOD32_WIDE);
    let sqr32: Vec<u32> = inputs.iter().map(|&v| reducer32.transform(v)).collect();

    let mut group = c.benchmark_group("proth sqr (2^31 + 1)");
    group.bench_function("FixedProth32", |b| {
        b.iter(|| {
            sqr32
                .iter()
                .map(|&v| reducer32.sqr(v))
                .reduce(|a, b| a.wrapping_add(b))
        })
    });
    group.finish();
}

pub fn bench_mul_near_width_u64(c: &mut Criterion) {
    let mut lhs: [u64; N_ENTRIES] = [0; N_ENTRIES];
    let mut rhs: [u64; N_ENTRIES] = [0; N_ENTRIES];
    for i in 0..N_ENTRIES {
        lhs[i] = random::<u64>() % MOD64_WIDE;
        rhs[i] = random::<u64>() % MOD64_WIDE;
    }

    let reducer64 = FixedProth64::<N_WIDE64, K_SMALL>::new(&MOD64_WIDE);
    let lhs64: Vec<u64> = lhs.iter().map(|&v| reducer64.transform(v)).collect();
    let rhs64: Vec<u64> = rhs.iter().map(|&v| reducer64.transform(v)).collect();

    let mut group = c.benchmark_group("proth mul (2^63 + 1)");
    group.bench_function("FixedProth64", |b| {
        b.iter(|| {
            lhs64
                .iter()
                .zip(rhs64.iter())
                .map(|(a, b)| reducer64.mul(a, b))
                .reduce(|a, b| a.wrapping_add(b))
        })
    });
    group.finish();
}

pub fn bench_sqr_near_width_u64(c: &mut Criterion) {
    let mut inputs: [u64; N_ENTRIES] = [0; N_ENTRIES];
    for i in 0..N_ENTRIES {
        inputs[i] = random::<u64>() % MOD64_WIDE;
    }

    let reducer64 = FixedProth64::<N_WIDE64, K_SMALL>::new(&MOD64_WIDE);
    let sqr64: Vec<u64> = inputs.iter().map(|&v| reducer64.transform(v)).collect();

    let mut group = c.benchmark_group("proth sqr (2^63 + 1)");
    group.bench_function("FixedProth64", |b| {
        b.iter(|| {
            sqr64
                .iter()
                .map(|&v| reducer64.sqr(v))
                .reduce(|a, b| a.wrapping_add(b))
        })
    });
    group.finish();
}

criterion_group!(
    benches,
    bench_transform,
    bench_add,
    bench_mul,
    bench_sqr,
    bench_inv,
    bench_mul_near_width_u32,
    bench_sqr_near_width_u32,
    bench_mul_near_width_u64,
    bench_sqr_near_width_u64,
);
criterion_main!(benches);
