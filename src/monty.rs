use crate::reduced::impl_reduced_binary_pow;
use crate::{udouble, ModularUnaryOps, Reducer, Vanilla};

/// Negated modular inverse on binary bases
trait NegModInv {
    /// Calculate `-(m^-1) mod R`, `R = 2^k. If m is odd, then result of m + 1 will be returned.
    fn neginv(m: &Self) -> Self;
}

// Entry i contains (2i+1)^(-1) mod 256.
const BINV_TABLE: [u8; 128] = [
    0x01, 0xAB, 0xCD, 0xB7, 0x39, 0xA3, 0xC5, 0xEF, 0xF1, 0x1B, 0x3D, 0xA7, 0x29, 0x13, 0x35, 0xDF,
    0xE1, 0x8B, 0xAD, 0x97, 0x19, 0x83, 0xA5, 0xCF, 0xD1, 0xFB, 0x1D, 0x87, 0x09, 0xF3, 0x15, 0xBF,
    0xC1, 0x6B, 0x8D, 0x77, 0xF9, 0x63, 0x85, 0xAF, 0xB1, 0xDB, 0xFD, 0x67, 0xE9, 0xD3, 0xF5, 0x9F,
    0xA1, 0x4B, 0x6D, 0x57, 0xD9, 0x43, 0x65, 0x8F, 0x91, 0xBB, 0xDD, 0x47, 0xC9, 0xB3, 0xD5, 0x7F,
    0x81, 0x2B, 0x4D, 0x37, 0xB9, 0x23, 0x45, 0x6F, 0x71, 0x9B, 0xBD, 0x27, 0xA9, 0x93, 0xB5, 0x5F,
    0x61, 0x0B, 0x2D, 0x17, 0x99, 0x03, 0x25, 0x4F, 0x51, 0x7B, 0x9D, 0x07, 0x89, 0x73, 0x95, 0x3F,
    0x41, 0xEB, 0x0D, 0xF7, 0x79, 0xE3, 0x05, 0x2F, 0x31, 0x5B, 0x7D, 0xE7, 0x69, 0x53, 0x75, 0x1F,
    0x21, 0xCB, 0xED, 0xD7, 0x59, 0xC3, 0xE5, 0x0F, 0x11, 0x3B, 0x5D, 0xC7, 0x49, 0x33, 0x55, 0xFF,
];

impl NegModInv for u8 {
    fn neginv(m: &Self) -> Self {
        let i = BINV_TABLE[((m >> 1) & 0x7F) as usize];
        i.wrapping_neg()
    }
}

impl NegModInv for u16 {
    fn neginv(m: &Self) -> Self {
        let mut i = BINV_TABLE[((m >> 1) & 0x7F) as usize] as u16;
        // hensel lifting
        i = 2u16.wrapping_sub(i.wrapping_mul(*m)).wrapping_mul(i);
        i.wrapping_neg()
    }
}

impl NegModInv for u32 {
    fn neginv(m: &Self) -> Self {
        let mut i = BINV_TABLE[((m >> 1) & 0x7F) as usize] as u32;
        i = 2u32.wrapping_sub(i.wrapping_mul(*m)).wrapping_mul(i);
        i = 2u32.wrapping_sub(i.wrapping_mul(*m)).wrapping_mul(i);
        i.wrapping_neg()
    }
}

impl NegModInv for u64 {
    fn neginv(m: &Self) -> Self {
        let mut i = BINV_TABLE[((m >> 1) & 0x7F) as usize] as u64;
        i = 2u64.wrapping_sub(i.wrapping_mul(*m)).wrapping_mul(i);
        i = 2u64.wrapping_sub(i.wrapping_mul(*m)).wrapping_mul(i);
        i = 2u64.wrapping_sub(i.wrapping_mul(*m)).wrapping_mul(i);
        i.wrapping_neg()
    }
}

impl NegModInv for u128 {
    fn neginv(m: &Self) -> Self {
        let mut i = BINV_TABLE[((m >> 1) & 0x7F) as usize] as u128;
        i = 2u128.wrapping_sub(i.wrapping_mul(*m)).wrapping_mul(i);
        i = 2u128.wrapping_sub(i.wrapping_mul(*m)).wrapping_mul(i);
        i = 2u128.wrapping_sub(i.wrapping_mul(*m)).wrapping_mul(i);
        i = 2u128.wrapping_sub(i.wrapping_mul(*m)).wrapping_mul(i);
        i.wrapping_neg()
    }
}

/// A modular reducer based on [Montgomery form](https://en.wikipedia.org/wiki/Montgomery_modular_multiplication#Montgomery_form), only supports odd modulus.
///
/// The generic type T represents the underlying integer representation for modular inverse `-m^-1 mod R`,
/// and `R=2^B` will be used as the auxiliary modulus, where B is automatically selected
/// based on the size of T.
#[derive(Debug, Clone, Copy)]
pub struct Montgomery<I, M>(I, M); // modular inverse of the modulus, modulus itself

macro_rules! impl_uprim_montgomery_reduce {
    ($t:ty, $double:ty) => {
        impl Montgomery<$t, $t> {
            fn reduce(&self, monty: $double) -> $t {
                debug_assert!(monty < ((self.1 as $double) << <$t>::BITS));

                // REDC algorithm
                let m = self.1;
                let tm = (monty as $t).wrapping_mul(self.0);
                let (t, overflow) = monty.overflowing_add((tm as $double) * (m as $double));
                let t = (t >> <$t>::BITS) as $t;

                if overflow {
                    t + m.wrapping_neg()
                } else if t >= m {
                    t - m
                } else {
                    t
                }
            }
        }
    };
}

macro_rules! impl_uprim_montgomery_core {
    ($single:ty) => {
        #[inline]
        fn new(m: &$single) -> Self {
            if m & 1 == 0 {
                panic!("Only odd modulus are supported by the Montgomery form");
            }
            Self(<$single>::neginv(m), *m)
        }
        #[inline(always)]
        fn modulus(&self) -> $single {
            self.1
        }
        #[inline(always)]
        fn is_zero(&self, target: &$single) -> bool {
            *target == 0
        }

        #[inline(always)]
        fn add(&self, lhs: $single, rhs: $single) -> $single {
            Vanilla::<$single>::new(&self.1).add(lhs, rhs)
        }

        #[inline(always)]
        fn double(&self, target: $single) -> $single {
            Vanilla::<$single>::new(&self.1).double(target)
        }

        #[inline(always)]
        fn sub(&self, lhs: $single, rhs: $single) -> $single {
            Vanilla::<$single>::new(&self.1).sub(lhs, rhs)
        }

        #[inline(always)]
        fn neg(&self, target: $single) -> $single {
            Vanilla::<$single>::new(&self.1).neg(target)
        }

        #[inline(always)]
        fn inv(&self, target: $single) -> Option<$single> {
            // TODO: support direct montgomery inverse
            // REF: http://cetinkayakoc.net/docs/j82.pdf
            self.residue(target)
                .invm(&self.1)
                .map(|v| self.transform(v))
        }

        impl_reduced_binary_pow!($single, $single);
    };
}

macro_rules! impl_uprim_montgomery {
    ($single:ty, $double:ty) => {
        impl_uprim_montgomery_reduce!($single, $double);

        impl Reducer<$single> for Montgomery<$single, $single> {
            impl_uprim_montgomery_core!($single);

            #[inline]
            fn transform(&self, target: $single) -> $single {
                (((target as $double) << <$single>::BITS) % (self.1 as $double)) as _
            }

            #[inline]
            fn residue(&self, target: $single) -> $single {
                self.reduce(target as $double)
            }

            #[inline]
            fn mul(&self, lhs: $single, rhs: $single) -> $single {
                self.reduce((lhs as $double) * (rhs as $double))
            }

            #[inline]
            fn square(&self, target: $single) -> $single {
                let d = target as $double;
                self.reduce(d * d)
            }
        }
    };
}

impl_uprim_montgomery!(u8, u16);
impl_uprim_montgomery!(u16, u32);
impl_uprim_montgomery!(u32, u64);
impl_uprim_montgomery!(u64, u128);

impl Montgomery<u128, u128> {
    fn reduce(&self, monty: udouble) -> u128 {
        debug_assert!(monty < udouble { hi: self.1, lo: 0 });

        // REDC algorithm
        let m = self.1;
        let tm = monty.lo.wrapping_mul(self.0);
        let (t, overflow) = monty.overflowing_add(udouble::widening_mul(tm, m));

        if overflow {
            t.hi + m.wrapping_neg()
        } else if t.hi >= m {
            t.hi - m
        } else {
            t.hi
        }
    }
}

impl Reducer<u128> for Montgomery<u128, u128> {
    #[inline]
    fn transform(&self, target: u128) -> u128 {
        if target == 0 {
            return 0;
        }
        udouble { hi: target, lo: 0 } % self.1
    }

    #[inline]
    fn residue(&self, target: u128) -> u128 {
        self.reduce(target.into())
    }

    #[inline]
    fn mul(&self, lhs: u128, rhs: u128) -> u128 {
        self.reduce(udouble::widening_mul(lhs, rhs))
    }

    #[inline]
    fn square(&self, target: u128) -> u128 {
        self.reduce(udouble::widening_square(target))
    }

    impl_uprim_montgomery_core!(u128);
}

// TODO: support bigints, use relaxed form described in https://cetinkayakoc.net/docs/j56.pdf and https://eprint.iacr.org/2011/239.pdf
// impl Reducer<BigUInt> for Montgomery<usize>
// or impl Reducer<&[usize]> for Montgomery<usize> ??
//
// We could directly base the operations on a specific bigint library (ibig is a good candidate), and convert all other bigint types to this type for arithmetics (rug::Integer::as_limbs, num_bigint::BigUint::to_u32_digits/to_u64_digits)
//     but there're two problems for ibig-rs now: it's unable to construct and deconstruct big integer by moving, and it has not implemented the Integer trait yet.
// So the better option is to create a standalone "reduce" function that takes &[usize] as input, and then call this reduce function for each big integer backend.
// And consider create separate ReducedInt type for small and bigints, as we can provide convenient interface for multi-by-single operators

// TODO(v0.5.x): accept even numbers by removing 2 factors from m and store the exponent
// Requirement: 1. A separate class to perform modular arithmetics with 2^n as modulus
//              2. Algorithm for construct residue from two components (see http://koclab.cs.ucsb.edu/teaching/cs154/docx/Notes7-Montgomery.pdf)
// Or we can just provide crt function, and let the implementation of monty int with full modulus support as an example code.

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{ModularCoreOps, ModularPow, ModularUnaryOps};
    use rand::random;

    const NRANDOM: u32 = 10;

    #[test]
    fn creation_test() {
        // a deterministic test case for u128
        let a = (0x81u128 << 120) - 1;
        let m = (0x81u128 << 119) - 1;
        let m = m >> m.trailing_zeros();
        let r = Montgomery::new(&m);
        assert_eq!(r.residue(r.transform(a)), a % m);

        // is_zero test
        let r = Montgomery::new(&11u8);
        assert!(r.is_zero(&r.transform(0)));
        let five = r.transform(5u8);
        let six = r.transform(6u8);
        assert!(r.is_zero(&r.add(five, six)));

        // random creation test
        for _ in 0..NRANDOM {
            let a = random::<u8>();
            let m = random::<u8>() | 1;
            let r = Montgomery::new(&m);
            assert_eq!(r.residue(r.transform(a)), a % m);

            let a = random::<u16>();
            let m = random::<u16>() | 1;
            let r = Montgomery::new(&m);
            assert_eq!(r.residue(r.transform(a)), a % m);

            let a = random::<u32>();
            let m = random::<u32>() | 1;
            let r = Montgomery::new(&m);
            assert_eq!(r.residue(r.transform(a)), a % m);

            let a = random::<u64>();
            let m = random::<u64>() | 1;
            let r = Montgomery::new(&m);
            assert_eq!(r.residue(r.transform(a)), a % m);

            let a = random::<u128>();
            let m = random::<u128>() | 1;
            let r = Montgomery::new(&m);
            assert_eq!(r.residue(r.transform(a)), a % m);
        }
    }

    #[test]
    fn test_against_prim() {
        macro_rules! tests_for {
            ($($T:ty)*) => ($(
                let m = random::<$T>() | 1;
                let r = Montgomery::new(&m);
                let e = random::<$T>() as $T;
                let (a, b) = (random::<$T>(), random::<$T>());
                let am = r.transform(a);
                let bm = r.transform(b);
                assert_eq!(r.residue(r.add(am, bm)), a.addm(b, &m));
                assert_eq!(r.residue(r.sub(am, bm)), a.subm(b, &m));
                assert_eq!(r.residue(r.mul(am, bm)), a.mulm(b, &m));
                assert_eq!(r.residue(r.neg(am)), a.negm(&m));
                assert_eq!(r.inv(am).map(|v| r.residue(v)), a.invm(&m));
                assert_eq!(r.residue(r.double(am)), a.dblm(&m));
                assert_eq!(r.residue(r.square(am)), a.sqm(&m));
                assert_eq!(r.residue(r.pow(am, e)), a.powm(e, &m));
            )*);
        }

        for _ in 0..NRANDOM {
            tests_for!(u8 u16 u32 u64 u128);
        }
    }
}
