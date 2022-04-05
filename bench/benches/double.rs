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
            lhs.iter().zip(rhs.iter())
                .map(|(&a, &b)| udouble::widening_mul(a, b))
                .reduce(|a, b| udouble::from(a.lo.wrapping_add(b.lo)))
        })
    });

    #[cfg(feature = "ethnum")]
    {
        use ethnum::U256;
        group.bench_function("ethnum", |b| {
            b.iter(|| {
                lhs.iter().zip(rhs.iter())
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
                lhs.iter().zip(rhs.iter())
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
                lhs.iter().zip(rhs.iter())
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
                lhs.iter().zip(rhs.iter())
                    .map(|(&a, &b)| U256::from(a) * U256::from(b))
                    .reduce(|a, b| U256::from(a.as_u128().wrapping_add(b.as_u128())))
            })
        });
    }

    group.finish();
}

criterion_group!(benches, bench_mul);
criterion_main!(benches);
