#[macro_use]
extern crate criterion;
use criterion::Criterion;
use rand::random;

use num_modular::udouble;

pub fn bench_mul(c: &mut Criterion) {
    let mut group = c.benchmark_group("single mul single");

    const N: usize = 12;
    let mut lhs: [u128; N] = [0; N];
    let mut rhs: [u128; N] = [0; N];
    for i in 0..N {
        lhs[i] = random();
        rhs[i] = random();
    }

    group.bench_function("ours", |b| {
        b.iter(|| {
            lhs.iter()
                .zip(rhs.iter())
                .map(|(&a, &b)| udouble::widening_mul(a, b))
                .reduce(|a, b| udouble::from(a.lo.wrapping_add(b.lo)))
        })
    });

    #[cfg(feature = "ethnum")]
    {
        use ethnum::U256;
        group.bench_function("ethnum", |b| {
            b.iter(|| {
                lhs.iter()
                    .zip(rhs.iter())
                    .map(|(&a, &b)| U256::from(a) * U256::from(b))
                    .reduce(|a, b| U256::from(a.0[0].wrapping_add(b.0[0])))
            })
        });
    }

    #[cfg(feature = "primitive-types")]
    {
        use primitive_types::U256;
        group.bench_function("primitive-types", |b| {
            b.iter(|| {
                lhs.iter()
                    .zip(rhs.iter())
                    .map(|(&a, &b)| U256::from(a) * U256::from(b))
                    .reduce(|a, b| U256::from(a.0[0].wrapping_add(b.0[0])))
            })
        });
    }

    #[cfg(feature = "uint")]
    {
        use uint::construct_uint;
        construct_uint! {
            pub struct U256(4);
        }
        group.bench_function("uint", |b| {
            b.iter(|| {
                lhs.iter()
                    .zip(rhs.iter())
                    .map(|(&a, &b)| U256::from(a) * U256::from(b))
                    .reduce(|a, b| U256::from(a.0[0].wrapping_add(b.0[0])))
            })
        });
    }

    #[cfg(feature = "zkp-u256")]
    {
        use zkp_u256::U256;
        group.bench_function("zkp-u256", |b| {
            b.iter(|| {
                lhs.iter()
                    .zip(rhs.iter())
                    .map(|(&a, &b)| U256::from(a) * U256::from(b))
                    .reduce(|a, b| U256::from(a.as_u128().wrapping_add(b.as_u128())))
            })
        });
    }

    #[cfg(feature = "crypto-bigint")]
    {
        use crypto_bigint::{Split, U128, U256};
        group.bench_function("crypto-bigint", |b| {
            b.iter(|| {
                lhs.iter()
                    .zip(rhs.iter())
                    .map(|(&a, &b)| U256::from(a).saturating_mul(&U256::from(b)))
                    .reduce(|a, b| U256::from((U128::ZERO, a.split().0.wrapping_add(&b.split().0))))
            })
        });
    }

    group.finish();
}

pub fn bench_div(c: &mut Criterion) {
    let mut group = c.benchmark_group("double div single");

    const N: usize = 12;
    let mut lhs: [(u128, u128); N] = [(0, 0); N];
    let mut rhs: [u128; N] = [0; N];
    for i in 0..N {
        lhs[i] = (random(), random());
        rhs[i] = random();
    }

    group.bench_function("ours", |b| {
        b.iter(|| {
            lhs.iter()
                .zip(rhs.iter())
                .map(|(&a, &b)| udouble { lo: a.0, hi: a.1 } / b)
                .reduce(|a, b| udouble::from(a.lo.wrapping_add(b.lo)))
        })
    });

    #[cfg(feature = "ethnum")]
    {
        use ethnum::U256;
        group.bench_function("ethnum", |b| {
            b.iter(|| {
                lhs.iter()
                    .zip(rhs.iter())
                    .map(|(&a, &b)| U256([a.0, a.1]) / b)
                    .reduce(|a, b| U256::from(a.0[0].wrapping_add(b.0[0])))
            })
        });
    }

    #[cfg(feature = "primitive-types")]
    {
        use primitive_types::U256;
        const MASK: u128 = (1 << 64) - 1;
        group.bench_function("primitive-types", |b| {
            b.iter(|| {
                lhs.iter()
                    .zip(rhs.iter())
                    .map(|(&a, &b)| {
                        U256([
                            (a.0 >> 64) as u64,
                            (a.0 & MASK) as u64,
                            (a.1 >> 64) as u64,
                            (a.1 & MASK) as u64,
                        ]) / b
                    })
                    .reduce(|a, b| U256::from(a.0[0].wrapping_add(b.0[0])))
            })
        });
    }

    #[cfg(feature = "uint")]
    {
        use uint::construct_uint;
        construct_uint! {
            pub struct U256(4);
        }
        const MASK: u128 = (1 << 64) - 1;
        group.bench_function("uint", |b| {
            b.iter(|| {
                lhs.iter()
                    .zip(rhs.iter())
                    .map(|(&a, &b)| {
                        U256([
                            (a.0 >> 64) as u64,
                            (a.0 & MASK) as u64,
                            (a.1 >> 64) as u64,
                            (a.1 & MASK) as u64,
                        ]) / b
                    })
                    .reduce(|a, b| U256::from(a.0[0].wrapping_add(b.0[0])))
            })
        });
    }

    #[cfg(feature = "zkp-u256")]
    {
        use zkp_u256::U256;
        group.bench_function("zkp-u256", |b| {
            b.iter(|| {
                lhs.iter()
                    .zip(rhs.iter())
                    .map(|(&a, &b)| ((U256::from(a.0) << 128) + U256::from(a.1)) / U256::from(b))
                    .reduce(|a, b| U256::from(a.as_u128().wrapping_add(b.as_u128())))
            })
        });
    }

    #[cfg(feature = "crypto-bigint")]
    {
        use crypto_bigint::{Split, U128, U256};
        group.bench_function("crypto-bigint", |b| {
            b.iter(|| {
                lhs.iter()
                    .zip(rhs.iter())
                    .map(|(&a, &b)| {
                        ((U256::from(a.0) << 128).wrapping_add(&U256::from(a.1)))
                            .div_rem(&U256::from(b))
                            .unwrap()
                            .0
                    })
                    .reduce(|a, b| U256::from((U128::ZERO, a.split().0.wrapping_add(&b.split().0))))
            })
        });
    }

    group.finish();
}

criterion_group!(benches, bench_mul, bench_div);
criterion_main!(benches);
