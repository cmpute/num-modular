use crate::reduced::impl_reduced_binary_pow;
use crate::{udouble, umax, ModularUnaryOps, Reducer, Vanilla};

// Proth primes: m = K * 2^N + 1 (K odd, K < 2^N)
//
// Montgomery REDC with R = 2^N, m' = R - 1 (since m ≡ 1 mod R):
//   REDC(v) = let v0 = v & mask;  hi = v >> N
//     if v0 == 0: result = hi
//     else:       result = hi + 1 + (R - v0) * K
//   All shift / multiply / add — no division in the fold.
//
// Values are stored in Montgomery form (a·R mod m).  One REDC + one
// final % normalises the result.

macro_rules! impl_fixed_proth {
    (
        $TypeName:ident,
        $T:ty,
        $D:ty,
        $max_N:expr,
        $kind:ident
    ) => {
        impl<const N: u8, const K: $T> $TypeName<N, K> {
            const BITMASK: $T = match (1 as $T).checked_shl(N as u32) {
                Some(v) => v.wrapping_sub(1),
                None => <$T>::MAX,
            };
            pub const MODULUS: $T = {
                let p2n = match (1 as $T).checked_shl(N as u32) {
                    Some(v) => v,
                    None => 0,
                };
                K.wrapping_mul(p2n).wrapping_add(1)
            };

            fn compute_r2() -> $T {
                impl_fixed_proth!(@compute_r2_body, $kind, $T, $D)
            }

            impl_fixed_proth!(@reduce_double, $kind, $T, $D);
        }

        impl<const N: u8, const K: $T> Reducer<$T> for $TypeName<N, K> {
            #[inline]
            fn new(m: &$T) -> Self {
                assert!(
                    *m == Self::MODULUS,
                    "the given modulus doesn't match with the generic params"
                );
                debug_assert!(N <= $max_N);
                debug_assert!(N > 0);
                debug_assert!(K > 0);
                debug_assert!(K % 2 == 1);
                debug_assert!((K as u128) < (1u128 << (N as u32)));
                debug_assert!(
                    (Self::MODULUS == 3 || Self::MODULUS % 3 != 0)
                        && (Self::MODULUS == 5 || Self::MODULUS % 5 != 0)
                        && (Self::MODULUS == 7 || Self::MODULUS % 7 != 0)
                        && (Self::MODULUS == 11 || Self::MODULUS % 11 != 0)
                        && (Self::MODULUS == 13 || Self::MODULUS % 13 != 0)
                );
                let r2 = Self::compute_r2();
                Self { r2 }
            }
            #[inline]
            fn transform(&self, target: $T) -> $T {
                if target == 0 {
                    return 0;
                }
                Self::reduce_double(impl_fixed_proth!(@widen_mul,
                    $kind, $T, $D, &target, &self.r2))
            }
            #[inline]
            fn check(&self, target: &$T) -> bool {
                *target < Self::MODULUS
            }
            #[inline]
            fn residue(&self, target: $T) -> $T {
                if target == 0 {
                    return 0;
                }
                Self::reduce_double(impl_fixed_proth!(@to_double, $kind, $D, target))
            }
            #[inline]
            fn modulus(&self) -> $T {
                Self::MODULUS
            }
            #[inline]
            fn is_zero(&self, target: &$T) -> bool {
                target == &0
            }

            #[inline]
            fn add(&self, lhs: &$T, rhs: &$T) -> $T {
                Vanilla::<$T>::add(&Self::MODULUS, *lhs, *rhs)
            }
            #[inline]
            fn sub(&self, lhs: &$T, rhs: &$T) -> $T {
                Vanilla::<$T>::sub(&Self::MODULUS, *lhs, *rhs)
            }
            #[inline]
            fn dbl(&self, target: $T) -> $T {
                Vanilla::<$T>::dbl(&Self::MODULUS, target)
            }
            #[inline]
            fn neg(&self, target: $T) -> $T {
                Vanilla::<$T>::neg(&Self::MODULUS, target)
            }
            #[inline]
            fn mul(&self, lhs: &$T, rhs: &$T) -> $T {
                Self::reduce_double(impl_fixed_proth!(@widen_mul, $kind, $T, $D, lhs, rhs))
            }
            #[inline]
            fn inv(&self, target: $T) -> Option<$T> {
                let plain = if target == 0 {
                    0
                } else {
                    Self::reduce_double(impl_fixed_proth!(@to_double, $kind, $D, target))
                };
                let inv_plain = if (N as u32) < usize::BITS {
                    (plain as usize)
                        .invm(&(Self::MODULUS as usize))
                        .map(|v| v as $T)
                } else {
                    plain.invm(&Self::MODULUS)
                }?;
                if inv_plain == 0 {
                    return Some(0);
                }
                Some(Self::reduce_double(impl_fixed_proth!(@widen_mul,
                    $kind, $T, $D, &inv_plain, &self.r2)))
            }
            #[inline]
            fn sqr(&self, target: $T) -> $T {
                Self::reduce_double(impl_fixed_proth!(@widen_sqr, $kind, $T, $D, target))
            }

            impl_reduced_binary_pow!($T);
        }
    };

    // Internal: compute R^2 mod m
    (@compute_r2_body, primitive, $T:ty, $D:ty) => {{
        let r = (1 as $D) << N;
        let r2 = r * r;
        (r2 % Self::MODULUS as $D) as $T
    }};
    (@compute_r2_body, udouble, $T:ty, $D:ty) => {{
        let r = udouble { hi: 0, lo: 1 } << N;
        let r2 = udouble::widening_square(r.lo);
        let m_ud = udouble { hi: 0, lo: Self::MODULUS };
        (r2 % m_ud).lo
    }};

    // Internal: reduce_double — one REDC (shift-based fold) + normalisation.
    // After REDC the result is < (K+1)·m.  For typical Proth primes K is
    // small (e.g. 1, 3, 5), so a few conditional subtractions outperform a
    // hardware division.  When K is large we fall back to %.
    (@reduce_double, primitive, $T:ty, $D:ty) => {
        pub fn reduce_double(v: $D) -> $T {
            let v0 = (v as $T) & Self::BITMASK;
            let v1 = v >> N;
            let mut acc: $D = if v0 == 0 {
                v1
            } else {
                let t = ((1 as $D) << N) - (v0 as $D);
                v1 + 1 + t * (K as $D)
            };
            // Normalise: result < (K+1)·m, at most K+1 subtractions of m.
            let m = Self::MODULUS as $D;
            let limit = (K as usize) + 1;
            for _ in 0..limit {
                if acc < m {
                    return acc as $T;
                }
                acc -= m;
            }
            (acc % m) as $T
        }
    };

    // Internal: reduce_double for udouble
    (@reduce_double, udouble, $T:ty, $D:ty) => {
        pub fn reduce_double(v: $D) -> $T {
            let v0 = v.lo & Self::BITMASK;
            let v1 = v >> N;
            let acc = if v0 == 0 {
                v1
            } else {
                let t = ((1 as umax) << N) - v0;
                let mut sum = v1;
                sum.lo += 1;
                if sum.lo < 1 {
                    sum.hi += 1;
                }
                let (tp, _) = udouble::widening_mul(t, K).overflowing_add(sum);
                tp
            };
            let m_ud = udouble { hi: 0, lo: Self::MODULUS };
            if acc.hi > 0 || acc.lo >= Self::MODULUS {
                (acc % m_ud).lo
            } else {
                acc.lo
            }
        }
    };

    // Convert T to D (primitive: as-cast, udouble: wrap in lo)
    (@to_double, primitive, $D:ty, $v:expr) => { $v as $D };
    (@to_double, udouble, $D:ty, $v:expr) => { udouble { hi: 0, lo: $v } };

    // Widening multiplication
    (@widen_mul, primitive, $T:ty, $D:ty, $lhs:expr, $rhs:expr) => {
        (*$lhs as $D) * (*$rhs as $D)
    };
    (@widen_mul, udouble, $T:ty, $D:ty, $lhs:expr, $rhs:expr) => {
        <$D>::widening_mul(*$lhs, *$rhs)
    };

    // Widening square
    (@widen_sqr, primitive, $T:ty, $D:ty, $target:expr) => {
        ($target as $D) * ($target as $D)
    };
    (@widen_sqr, udouble, $T:ty, $D:ty, $target:expr) => {
        <$D>::widening_square($target)
    };
}

/// A modular reducer for Proth primes `K * 2^N + 1` with 32-bit operands.
///
/// Supports `N` up to 31, `K` odd with `K < 2^N`. All inputs and outputs are `u32`.
/// The modulus `K * 2^N + 1` must be prime for modular inverse and Fermat-based
/// operations to be valid.  Uses Montgomery form with `R = 2^N` internally.
///
/// # Example
///
/// ```rust
/// use num_modular::{FixedProth32, Reducer};
///
/// const N: u8 = 4;
/// const K: u32 = 1;
/// let modulus = K * (1u32 << N) + 1; // 1*2^4 + 1 = 17
/// let reducer = FixedProth32::<N, K>::new(&modulus);
/// let a = reducer.transform(3);
/// let b = reducer.transform(5);
/// assert_eq!(reducer.residue(reducer.add(&a, &b)), 8);
/// assert_eq!(reducer.residue(reducer.mul(&a, &b)), 15);
/// ```
#[derive(Debug, Clone, Copy)]
pub struct FixedProth32<const N: u8, const K: u32> {
    r2: u32,
}

impl_fixed_proth!(FixedProth32, u32, u64, 31, primitive);

/// A modular reducer for Proth primes `K * 2^N + 1` with 64-bit operands.
///
/// Supports `N` up to 63, `K` odd with `K < 2^N`. All inputs and outputs are `u64`.
/// Uses Montgomery form with `R = 2^N` internally.
///
/// # Example
///
/// ```rust
/// use num_modular::{FixedProth64, Reducer};
///
/// const N: u8 = 5;
/// const K: u64 = 3;
/// let modulus = K * (1u64 << N) + 1; // 3*2^5 + 1 = 97
/// let reducer = FixedProth64::<N, K>::new(&modulus);
/// let a = reducer.transform(10);
/// let b = reducer.transform(20);
/// assert_eq!(reducer.residue(reducer.mul(&a, &b)), (10u64 * 20) % 97);
/// ```
#[derive(Debug, Clone, Copy)]
pub struct FixedProth64<const N: u8, const K: u64> {
    r2: u64,
}

impl_fixed_proth!(FixedProth64, u64, u128, 63, primitive);

/// A modular reducer for Proth primes `K * 2^N + 1`.
///
/// Supports `N` up to 127, `K` odd with `K < 2^N`. All inputs and outputs are [umax] (currently `u128`).
/// Uses Montgomery form with `R = 2^N` internally.
///
/// # Example
///
/// ```rust
/// use num_modular::{FixedProth, Reducer};
///
/// const N: u8 = 16;
/// const K: u128 = 1;
/// let modulus = K * (1u128 << N) + 1; // 2^16 + 1 = 65537
/// let reducer = FixedProth::<N, K>::new(&modulus);
/// let a = reducer.transform(1000);
/// let b = reducer.transform(2000);
/// assert_eq!(reducer.residue(reducer.mul(&a, &b)), (1000u128 * 2000) % modulus);
/// ```
#[derive(Debug, Clone, Copy)]
pub struct FixedProth<const N: u8, const K: umax> {
    r2: umax,
}

impl_fixed_proth!(FixedProth, umax, udouble, 127, udouble);

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{ModularCoreOps, ModularPow};
    use rand::random;

    // u128 types
    type P128_1 = FixedProth<2, 1>;   // m = 5
    type P128_2 = FixedProth<4, 1>;   // m = 17
    type P128_3 = FixedProth<5, 3>;   // m = 97
    type P128_4 = FixedProth<8, 3>;   // m = 769
    type P128_5 = FixedProth<16, 1>;  // m = 65537

    // u64 types
    type P64_1 = FixedProth64<4, 1>;   // m = 17
    type P64_2 = FixedProth64<5, 3>;   // m = 97
    type P64_3 = FixedProth64<8, 1>;   // m = 257
    type P64_4 = FixedProth64<16, 1>;  // m = 65537

    // u32 types
    type P32_1 = FixedProth32<2, 1>;  // m = 5
    type P32_2 = FixedProth32<2, 3>;  // m = 13
    type P32_3 = FixedProth32<4, 1>;  // m = 17
    type P32_4 = FixedProth32<3, 5>;  // m = 41

    const NRANDOM: u32 = 10;

    #[test]
    fn creation_test_u128() {
        for _ in 0..NRANDOM {
            let a = random::<u128>();

            const M1: u128 = <P128_1>::MODULUS;
            let r1 = P128_1::new(&M1);
            assert_eq!(r1.residue(r1.transform(a % M1)), a % M1);

            const M2: u128 = <P128_2>::MODULUS;
            let r2 = P128_2::new(&M2);
            assert_eq!(r2.residue(r2.transform(a % M2)), a % M2);

            const M3: u128 = <P128_3>::MODULUS;
            let r3 = P128_3::new(&M3);
            assert_eq!(r3.residue(r3.transform(a % M3)), a % M3);

            const M4: u128 = <P128_4>::MODULUS;
            let r4 = P128_4::new(&M4);
            assert_eq!(r4.residue(r4.transform(a % M4)), a % M4);

            const M5: u128 = <P128_5>::MODULUS;
            let r5 = P128_5::new(&M5);
            assert_eq!(r5.residue(r5.transform(a % M5)), a % M5);
        }
    }

    #[test]
    fn creation_test_u64() {
        for _ in 0..NRANDOM {
            let a = random::<u64>();

            const M1: u64 = <P64_1>::MODULUS;
            let r1 = P64_1::new(&M1);
            assert_eq!(r1.residue(r1.transform(a % M1)), a % M1);

            const M2: u64 = <P64_2>::MODULUS;
            let r2 = P64_2::new(&M2);
            assert_eq!(r2.residue(r2.transform(a % M2)), a % M2);

            const M3: u64 = <P64_3>::MODULUS;
            let r3 = P64_3::new(&M3);
            assert_eq!(r3.residue(r3.transform(a % M3)), a % M3);

            const M4: u64 = <P64_4>::MODULUS;
            let r4 = P64_4::new(&M4);
            assert_eq!(r4.residue(r4.transform(a % M4)), a % M4);
        }
    }

    #[test]
    fn creation_test_u32() {
        for _ in 0..NRANDOM {
            let a = random::<u32>();

            const M1: u32 = <P32_1>::MODULUS;
            let r1 = P32_1::new(&M1);
            assert_eq!(r1.residue(r1.transform(a % M1)), a % M1);

            const M2: u32 = <P32_2>::MODULUS;
            let r2 = P32_2::new(&M2);
            assert_eq!(r2.residue(r2.transform(a % M2)), a % M2);

            const M3: u32 = <P32_3>::MODULUS;
            let r3 = P32_3::new(&M3);
            assert_eq!(r3.residue(r3.transform(a % M3)), a % M3);

            const M4: u32 = <P32_4>::MODULUS;
            let r4 = P32_4::new(&M4);
            assert_eq!(r4.residue(r4.transform(a % M4)), a % M4);
        }
    }

    #[test]
    fn test_against_modops_u128() {
        macro_rules! tests_for {
            ($a:ident, $b:ident, $e:ident; $($M:ty)*) => ($({
                const P: u128 = <$M>::MODULUS;
                let r = <$M>::new(&P);
                let am = r.transform($a);
                let bm = r.transform($b);
                assert_eq!(r.residue(r.add(&am, &bm)), $a.addm($b, &P));
                assert_eq!(r.residue(r.sub(&am, &bm)), $a.subm($b, &P));
                assert_eq!(r.residue(r.mul(&am, &bm)), $a.mulm($b, &P));
                assert_eq!(r.residue(r.neg(am)), $a.negm(&P));
                assert_eq!(r.residue(r.dbl(am)), $a.dblm(&P));
                assert_eq!(r.residue(r.sqr(am)), $a.sqm(&P));
                assert_eq!(r.residue(r.pow(am, &$e)), $a.powm($e, &P));
                if let (Some(inv), Some(ref_inv)) = (r.inv(am), $a.invm(&P)) {
                    assert_eq!(r.residue(inv), ref_inv);
                }
            })*);
        }

        for _ in 0..NRANDOM {
            let a = random::<u128>();
            let b = random::<u128>();
            let e = random::<u8>() as u128;
            tests_for!(a, b, e; P128_1 P128_2 P128_3 P128_4 P128_5);
        }
    }

    #[test]
    fn test_against_modops_u64() {
        macro_rules! tests_for {
            ($a:ident, $b:ident, $e:ident; $($M:ty)*) => ($({
                const P: u64 = <$M>::MODULUS;
                let r = <$M>::new(&P);
                let am = r.transform($a);
                let bm = r.transform($b);
                assert_eq!(r.residue(r.add(&am, &bm)), $a.addm($b, &P));
                assert_eq!(r.residue(r.sub(&am, &bm)), $a.subm($b, &P));
                assert_eq!(r.residue(r.mul(&am, &bm)), $a.mulm($b, &P));
                assert_eq!(r.residue(r.neg(am)), $a.negm(&P));
                assert_eq!(r.residue(r.dbl(am)), $a.dblm(&P));
                assert_eq!(r.residue(r.sqr(am)), $a.sqm(&P));
                assert_eq!(r.residue(r.pow(am, &$e)), $a.powm($e, &P));
                if let (Some(inv), Some(ref_inv)) = (r.inv(am), $a.invm(&P)) {
                    assert_eq!(r.residue(inv), ref_inv);
                }
            })*);
        }

        for _ in 0..NRANDOM {
            let a = random::<u64>();
            let b = random::<u64>();
            let e = random::<u8>() as u64;
            tests_for!(a, b, e; P64_1 P64_2 P64_3 P64_4);
        }
    }

    #[test]
    fn test_against_modops_u32() {
        macro_rules! tests_for {
            ($a:ident, $b:ident, $e:ident; $($M:ty)*) => ($({
                const P: u32 = <$M>::MODULUS;
                let r = <$M>::new(&P);
                let am = r.transform($a);
                let bm = r.transform($b);
                assert_eq!(r.residue(r.add(&am, &bm)), $a.addm($b, &P));
                assert_eq!(r.residue(r.sub(&am, &bm)), $a.subm($b, &P));
                assert_eq!(r.residue(r.mul(&am, &bm)), $a.mulm($b, &P));
                assert_eq!(r.residue(r.neg(am)), $a.negm(&P));
                assert_eq!(r.residue(r.dbl(am)), $a.dblm(&P));
                assert_eq!(r.residue(r.sqr(am)), $a.sqm(&P));
                assert_eq!(r.residue(r.pow(am, &$e)), $a.powm($e, &P));
                if let (Some(inv), Some(ref_inv)) = (r.inv(am), $a.invm(&P)) {
                    assert_eq!(r.residue(inv), ref_inv);
                }
            })*);
        }

        for _ in 0..NRANDOM {
            let a = random::<u32>();
            let b = random::<u32>();
            let e = random::<u8>() as u32;
            tests_for!(a, b, e; P32_1 P32_2 P32_3 P32_4);
        }
    }

    #[test]
    fn test_add_near_overflow_u64() {
        type S = FixedProth64<32, 3>;
        const M: u64 = <S>::MODULUS;
        let r = S::new(&M);

        let a = M - 1;
        let b = M - 2;
        let am = r.transform(a);
        let bm = r.transform(b);
        let sum = r.add(&am, &bm);
        assert_eq!(r.residue(sum), a.addm(b, &M));

        let a2 = M - 1;
        let a2m = r.transform(a2);
        let dbl = r.dbl(a2m);
        assert_eq!(r.residue(dbl), a2.dblm(&M));
    }
}
