use crate::reduced::{impl_reduced_binary_pow, impl_reduced_ops};
use crate::{imax, udouble, umax, ModularUnaryOps, Reducer};

// REF: Handbook of Cryptography 14.3.4

macro_rules! impl_fixed_trinomial_solinas {
    (
        $TypeName:ident,
        $T:ty,
        $K:ty,
        $D:ty,
        $half_bits:expr,
        $max_P1:expr,
        $kind:ident
    ) => {
        impl<const P1: u8, const P2: u8, const K: $K> $TypeName<P1, P2, K> {
            const BITMASK: $T = match (1 as $T).checked_shl(P1 as u32) {
                Some(v) => v.wrapping_sub(1),
                None => <$T>::MAX,
            };
            pub const MODULUS: $T = {
                let p1 = match (1 as $T).checked_shl(P1 as u32) {
                    Some(v) => v,
                    None => 0,
                };
                let p2 = match (1 as $T).checked_shl(P2 as u32) {
                    Some(v) => v,
                    None => panic!("P2 exceeds type width"),
                };
                if K >= 0 {
                    p1.wrapping_sub(p2).wrapping_add(K as $T)
                } else {
                    p1.wrapping_sub(p2).wrapping_sub((-K) as $T)
                }
            };

            /// Worst-case fold count for `reduce_double`.
            /// Each fold removes roughly (P1−P2) bits; ⌈P1/(P1−P2)⌉ folds
            /// shrink from 2·P1 bits to ≤ P1, plus 1 (K>0) or 2 (K<0) for the carry tail.
            const FOLDS: u32 = {
                let gap = (P1 - P2) as u32;
                let folds_ceil = ((P1 as u32) + gap - 1) / gap;
                if K > 0 {
                    folds_ceil + 1
                } else if K < 0 {
                    folds_ceil + 2
                } else {
                    1 // K == 0: trivial reduction, single fold
                }
            };

            impl_fixed_trinomial_solinas!(@reduce_single, $kind, $T, $D);
            impl_fixed_trinomial_solinas!(@reduce_double, $kind, $T, $D);
        }

        impl<const P1: u8, const P2: u8, const K: $K> Reducer<$T> for $TypeName<P1, P2, K> {
            #[inline]
            fn new(m: &$T) -> Self {
                assert!(
                    *m == Self::MODULUS,
                    "the given modulus doesn't match with the generic params"
                );
                debug_assert!(P1 <= $max_P1);
                debug_assert!(P2 > 0 && P1 > P2);
                debug_assert!(K % 2 != 0); // modulus must be odd
                // |K| < 2^P2 keeps each reduction step non-negative in Z (required for unsigned arithmetic)
                debug_assert!((K.unsigned_abs() as u128) < (1u128 << (P2 as u32)));
                debug_assert!(
                    (Self::MODULUS == 3 || Self::MODULUS % 3 != 0)
                        && (Self::MODULUS == 5 || Self::MODULUS % 5 != 0)
                        && (Self::MODULUS == 7 || Self::MODULUS % 7 != 0)
                        && (Self::MODULUS == 11 || Self::MODULUS % 11 != 0)
                        && (Self::MODULUS == 13 || Self::MODULUS % 13 != 0)
                ); // error on easy composites
                Self {}
            }
            #[inline]
            fn transform(&self, target: $T) -> $T {
                Self::reduce_single(target)
            }
            #[inline]
            fn residue(&self, target: $T) -> $T {
                target
            }

            impl_reduced_ops!($T);

            #[inline]
            fn mul(&self, lhs: &$T, rhs: &$T) -> $T {
                if (P1 as u32) < $half_bits {
                    Self::reduce_single(lhs * rhs)
                } else {
                    Self::reduce_double(impl_fixed_trinomial_solinas!(@widen_mul, $kind, $T, $D, lhs, rhs))
                }
            }
            #[inline]
            fn inv(&self, target: $T) -> Option<$T> {
                // TODO: inv can be specialized
                // REF: https://xn--2-umb.com/22/goldilocks/
                if (P1 as u32) < usize::BITS {
                    (target as usize)
                        .invm(&(Self::MODULUS as usize))
                        .map(|v| v as $T)
                } else {
                    target.invm(&Self::MODULUS)
                }
            }
            #[inline]
            fn sqr(&self, target: $T) -> $T {
                if (P1 as u32) < $half_bits {
                    Self::reduce_single(target * target)
                } else {
                    Self::reduce_double(impl_fixed_trinomial_solinas!(@widen_sqr, $kind, $T, $D, target))
                }
            }

            impl_reduced_binary_pow!($T);
        }
    };

    // Internal: reduce_single for primitive double-width types (u32→u64, u64→u128)
    (@reduce_single, primitive, $T:ty, $D:ty) => {
        /// Reduces a single-width value `v` modulo `2^P1 - 2^P2 + K`.
        ///
        /// For the result of a widening multiplication or square, use
        /// [`reduce_double`](Self::reduce_double) instead.
        pub const fn reduce_single(v: $T) -> $T {
            let mut v: $D = v as $D;
            while v >> P1 > 0 {
                let lo = (v as $T) & Self::BITMASK;
                let hi = v >> P1;
                let mut sum: $D = (hi << (P2 as u32)) + (lo as $D);
                if K > 0 {
                    sum -= hi * (K as $D);
                } else if K < 0 {
                    sum += hi * ((-K) as $D);
                }
                v = sum;
            }
            let v = v as $T;
            if v >= Self::MODULUS {
                v - Self::MODULUS
            } else {
                v
            }
        }
    };

    // Internal: reduce_single for udouble (umax→udouble). Stays in udouble for the same reason
    // as reduce_double below: `hi << P2` can exceed `umax` during the tail.
    (@reduce_single, udouble, $T:ty, $D:ty) => {
        /// Reduces a single-width value `v` modulo `2^P1 - 2^P2 + K`.
        ///
        /// For the result of a widening multiplication or square, use
        /// [`reduce_double`](Self::reduce_double) instead.
        pub fn reduce_single(v: $T) -> $T {
            let mut v: $D = udouble { hi: 0, lo: v };
            while v.hi > 0 || v.lo >> P1 > 0 {
                let lo = v.lo & Self::BITMASK;
                let hi = v >> P1;
                let mut sum = (hi << (P2 as u32)) + lo;
                if K > 0 {
                    sum -= hi * (K as umax);
                } else if K < 0 {
                    sum += hi * ((-K) as umax);
                }
                v = sum;
            }
            let v = v.lo;
            if v >= Self::MODULUS {
                v - Self::MODULUS
            } else {
                v
            }
        }
    };

    // Internal: reduce_double for primitive double-width types (u32→u64, u64→u128)
    //
    // When the worst-case fold count is small, replace the while loop with
    // straight-line unconditional folds. Each fold is a no-op once hi reaches 0.
    // FOLDS from the expert formula: ⌈P1/(P1−P2)⌉ + 1 (K>0) or +2 (K<0).
    // Unrolling condition: P2 ≤ ⌊2·P1/3⌋  ⇔  FOLDS ≤ 4.
    (@reduce_double, primitive, $T:ty, $D:ty) => {
        /// Reduces a double-width value `v` modulo `2^P1 - 2^P2 + K`.
        ///
        /// This handles widening-multiplication or widening-square results.
        /// For single-width values, use [`reduce_single`](Self::reduce_single).
        pub fn reduce_double(v: $D) -> $T {
            let mut lo = (v as $T) & Self::BITMASK;
            let mut hi = v >> P1;
            macro_rules! solinas_fold {
                () => {
                    let mut sum: $D = (hi << (P2 as u32)) + (lo as $D);
                    if K > 0 { sum -= hi * (K as $D); }
                    else if K < 0 { sum += hi * ((-K) as $D); }
                    lo = (sum as $T) & Self::BITMASK;
                    hi = sum >> P1;
                };
            }
            if Self::FOLDS <= 3 {
                #[allow(unused_assignments)] { solinas_fold!(); }
                #[allow(unused_assignments)] { solinas_fold!(); }
                #[allow(unused_assignments)] { solinas_fold!(); }
            } else if Self::FOLDS == 4 {
                #[allow(unused_assignments)] { solinas_fold!(); }
                #[allow(unused_assignments)] { solinas_fold!(); }
                #[allow(unused_assignments)] { solinas_fold!(); }
                #[allow(unused_assignments)] { solinas_fold!(); }
            } else {
                while hi > 0 { solinas_fold!(); }
            }
            if lo >= Self::MODULUS {
                lo - Self::MODULUS
            } else {
                lo
            }
        }
    };

    // Internal: reduce_double for udouble (u128→udouble)
    //
    // Unlike [Mersenne](crate::FixedMersenne)'s two-phase loop (udouble while `hi.hi > 0`, then
    // `umax` while `hi.lo > 0`), Solinas keeps `hi` as [udouble] until fully zero. Mersenne's
    // tail step is `hi * K + lo`, which stays within `umax` when `K < 2^(P-1)`. Solinas uses
    // `hi << P2`, which can exceed `umax` even when `hi` fits in one word (e.g. `hi * 2^P2`), so
    // the tail must stay in double-width arithmetic.
    (@reduce_double, udouble, $T:ty, $D:ty) => {
        /// Reduces a double-width value `v` modulo `2^P1 - 2^P2 + K`.
        ///
        /// This handles widening-multiplication or widening-square results.
        /// For single-width values, use [`reduce_single`](Self::reduce_single).
        pub fn reduce_double(v: $D) -> $T {
            let mut lo = v.lo & Self::BITMASK;
            let mut hi = v >> P1;
            macro_rules! udouble_fold {
                () => {
                    let mut sum = (hi << (P2 as u32)) + lo;
                    if K > 0 { sum -= hi * (K as umax); }
                    else if K < 0 { sum += hi * ((-K) as umax); }
                    lo = sum.lo & Self::BITMASK;
                    hi = sum >> P1;
                };
            }
            if Self::FOLDS <= 3 {
                #[allow(unused_assignments)] { udouble_fold!(); }
                #[allow(unused_assignments)] { udouble_fold!(); }
                #[allow(unused_assignments)] { udouble_fold!(); }
            } else if Self::FOLDS == 4 {
                #[allow(unused_assignments)] { udouble_fold!(); }
                #[allow(unused_assignments)] { udouble_fold!(); }
                #[allow(unused_assignments)] { udouble_fold!(); }
                #[allow(unused_assignments)] { udouble_fold!(); }
            } else {
                while hi.hi > 0 || hi.lo > 0 { udouble_fold!(); }
            }
            if lo >= Self::MODULUS {
                lo - Self::MODULUS
            } else {
                lo
            }
        }
    };

    // Internal: widening multiplication for primitive types
    (@widen_mul, primitive, $T:ty, $D:ty, $lhs:expr, $rhs:expr) => {
        (*$lhs as $D) * (*$rhs as $D)
    };

    // Internal: widening multiplication for udouble
    (@widen_mul, udouble, $T:ty, $D:ty, $lhs:expr, $rhs:expr) => {
        <$D>::widening_mul(*$lhs, *$rhs)
    };

    // Internal: widening square for primitive types
    (@widen_sqr, primitive, $T:ty, $D:ty, $target:expr) => {
        ($target as $D) * ($target as $D)
    };

    // Internal: widening square for udouble
    (@widen_sqr, udouble, $T:ty, $D:ty, $target:expr) => {
        <$D>::widening_square($target)
    };
}

/// A modular reducer for trinomial Solinas numbers `2^P1 - 2^P2 + K` as modulus with 32-bit operands.
///
/// Supports `P1` up to 32, `P2 < P1`, and odd signed `K` with `|K| < 2^P2`. All inputs and outputs are `u32`.
/// The modulus `2^P1 - 2^P2 + K` must be prime for modular inverse and Fermat-based operations to be valid.
///
/// # Example
///
/// ```rust
/// use num_modular::{FixedTrinomialSolinas32, Reducer};
///
/// const P1: u8 = 4;
/// const P2: u8 = 2;
/// const K: i32 = 1;
/// let modulus = (1u32 << P1) - (1u32 << P2) + (K as u32); // 2^4 - 2^2 + 1 = 13
/// let reducer = FixedTrinomialSolinas32::<P1, P2, K>::new(&modulus);
/// let a = reducer.transform(3);
/// let b = reducer.transform(5);
/// assert_eq!(reducer.residue(reducer.add(&a, &b)), 8);
/// ```
#[must_use]
#[derive(Debug, Clone, Copy)]
pub struct FixedTrinomialSolinas32<const P1: u8, const P2: u8, const K: i32>();

impl_fixed_trinomial_solinas!(FixedTrinomialSolinas32, u32, i32, u64, 16, 32, primitive);

/// A modular reducer for trinomial Solinas numbers `2^P1 - 2^P2 + K` as modulus with 64-bit operands.
///
/// Supports `P1` up to 64, `P2 < P1`, and odd signed `K` with `|K| < 2^P2`. All inputs and outputs are `u64`.
/// Uses `u128` as the double-width intermediate for multiplication and reduction.
/// The modulus `2^P1 - 2^P2 + K` must be prime for modular inverse and Fermat-based operations to be valid.
///
/// # Example
///
/// ```rust
/// use num_modular::{FixedTrinomialSolinas64, Reducer};
///
/// const P1: u8 = 6;
/// const P2: u8 = 2;
/// const K: i64 = 1;
/// let modulus = (1u64 << P1) - (1u64 << P2) + (K as u64); // 2^6 - 2^2 + 1 = 61
/// let reducer = FixedTrinomialSolinas64::<P1, P2, K>::new(&modulus);
/// let a = reducer.transform(10);
/// let b = reducer.transform(20);
/// assert_eq!(reducer.residue(reducer.mul(&a, &b)), (10u64 * 20) % 61);
/// ```
#[must_use]
#[derive(Debug, Clone, Copy)]
pub struct FixedTrinomialSolinas64<const P1: u8, const P2: u8, const K: i64>();

impl_fixed_trinomial_solinas!(FixedTrinomialSolinas64, u64, i64, u128, 32, 64, primitive);

/// A modular reducer for trinomial Solinas numbers `2^P1 - 2^P2 + K` as modulus.
///
/// Supports `P1` up to 127, `P2 < P1`, and odd signed `K` with `|K| < 2^P2`. All inputs and outputs are [umax] (currently `u128`).
///
/// The `P1` is limited to 127 so that overflow checks aren't necessary. This covers all trinomial
/// Solinas primes within the range of [umax] (i.e. `u128`).
///
/// # Example
///
/// ```rust
/// use num_modular::{FixedTrinomialSolinas, Reducer};
///
/// const P1: u8 = 31;
/// const P2: u8 = 13;
/// const K: i128 = 1;
/// let modulus = (1u128 << P1) - (1u128 << P2) + (K as u128);
/// let reducer = FixedTrinomialSolinas::<P1, P2, K>::new(&modulus);
/// let a = reducer.transform(1000);
/// let b = reducer.transform(2000);
/// assert_eq!(reducer.residue(reducer.mul(&a, &b)), (1000u128 * 2000) % modulus);
/// ```
#[must_use]
#[derive(Debug, Clone, Copy)]
pub struct FixedTrinomialSolinas<const P1: u8, const P2: u8, const K: imax>();

impl_fixed_trinomial_solinas!(FixedTrinomialSolinas, umax, imax, udouble, 64, 127, udouble);

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{ModularCoreOps, ModularPow};
    use rand::random;

    // u128 types
    type S1 = FixedTrinomialSolinas<31, 13, 1>;
    type S2 = FixedTrinomialSolinas<61, 30, 1>;
    type S3 = FixedTrinomialSolinas<127, 64, 1>;
    type S4 = FixedTrinomialSolinas<32, 16, 1>;
    type S5 = FixedTrinomialSolinas<56, 28, -1>;
    type S6 = FixedTrinomialSolinas<122, 61, -3>;

    // u64 types
    type S64_1 = FixedTrinomialSolinas64<31, 13, 1>;
    type S64_2 = FixedTrinomialSolinas64<61, 30, 1>;
    type S64_3 = FixedTrinomialSolinas64<32, 16, 1>;
    type S64_4 = FixedTrinomialSolinas64<64, 32, 1>; // 2^64 - 2^32 + 1

    // u32 types
    type S32_1 = FixedTrinomialSolinas32<4, 2, 1>;
    type S32_2 = FixedTrinomialSolinas32<5, 3, -1>;
    type S32_3 = FixedTrinomialSolinas32<6, 2, 1>;
    type S32_4 = FixedTrinomialSolinas32<32, 20, 1>;

    const NRANDOM: u32 = 10;

    #[test]
    fn creation_test_u128() {
        const P: umax = <S1>::MODULUS;
        let m = S1::new(&P);
        assert_eq!(m.residue(m.transform(0)), 0);
        assert_eq!(m.residue(m.transform(1)), 1);
        assert_eq!(m.residue(m.transform(P)), 0);
        assert_eq!(m.residue(m.transform(P - 1)), P - 1);
        assert_eq!(m.residue(m.transform(P + 1)), 1);

        for _ in 0..NRANDOM {
            let a = random::<umax>();

            const P1: umax = <S1>::MODULUS;
            let m1 = S1::new(&P1);
            assert_eq!(m1.residue(m1.transform(a)), a % P1);
            const P2: umax = <S2>::MODULUS;
            let m2 = S2::new(&P2);
            assert_eq!(m2.residue(m2.transform(a)), a % P2);
            const P3: umax = <S3>::MODULUS;
            let m3 = S3::new(&P3);
            assert_eq!(m3.residue(m3.transform(a)), a % P3);
            const P4: umax = <S4>::MODULUS;
            let m4 = S4::new(&P4);
            assert_eq!(m4.residue(m4.transform(a)), a % P4);
            const P5: umax = <S5>::MODULUS;
            let m5 = S5::new(&P5);
            assert_eq!(m5.residue(m5.transform(a)), a % P5);
            const P6: umax = <S6>::MODULUS;
            let m6 = S6::new(&P6);
            assert_eq!(m6.residue(m6.transform(a)), a % P6);
        }
    }

    #[test]
    fn creation_test_u64() {
        for _ in 0..NRANDOM {
            let a = random::<u64>();

            const P1: u64 = <S64_1>::MODULUS;
            let m1 = S64_1::new(&P1);
            assert_eq!(m1.residue(m1.transform(a)), a % P1);
            const P2: u64 = <S64_2>::MODULUS;
            let m2 = S64_2::new(&P2);
            assert_eq!(m2.residue(m2.transform(a)), a % P2);
            const P3: u64 = <S64_3>::MODULUS;
            let m3 = S64_3::new(&P3);
            assert_eq!(m3.residue(m3.transform(a)), a % P3);
            const P4: u64 = <S64_4>::MODULUS;
            let m4 = S64_4::new(&P4);
            assert_eq!(m4.residue(m4.transform(a)), a % P4);
        }
    }

    #[test]
    fn creation_test_u32() {
        for _ in 0..NRANDOM {
            let a = random::<u32>();

            const P1: u32 = <S32_1>::MODULUS;
            let m1 = S32_1::new(&P1);
            assert_eq!(m1.residue(m1.transform(a)), a % P1);
            const P2: u32 = <S32_2>::MODULUS;
            let m2 = S32_2::new(&P2);
            assert_eq!(m2.residue(m2.transform(a)), a % P2);
            const P3: u32 = <S32_3>::MODULUS;
            let m3 = S32_3::new(&P3);
            assert_eq!(m3.residue(m3.transform(a)), a % P3);
            const P4: u32 = <S32_4>::MODULUS;
            let m4 = S32_4::new(&P4);
            assert_eq!(m4.residue(m4.transform(a)), a % P4);
        }
    }

    #[test]
    fn test_against_modops_u128() {
        macro_rules! tests_for {
            ($a:tt, $b:tt, $e:tt; $($M:ty)*) => ($({
                const P: umax = <$M>::MODULUS;
                let r = <$M>::new(&P);
                let am = r.transform($a);
                let bm = r.transform($b);
                assert_eq!(r.add(&am, &bm), $a.addm($b, &P));
                assert_eq!(r.sub(&am, &bm), $a.subm($b, &P));
                assert_eq!(r.mul(&am, &bm), $a.mulm($b, &P));
                assert_eq!(r.neg(am), $a.negm(&P));
                assert_eq!(r.inv(am), $a.invm(&P));
                assert_eq!(r.dbl(am), $a.dblm(&P));
                assert_eq!(r.sqr(am), $a.sqm(&P));
                assert_eq!(r.pow(am, &$e), $a.powm($e, &P));
            })*);
        }

        for _ in 0..NRANDOM {
            let (a, b) = (random::<u128>(), random::<u128>());
            let e = random::<u8>() as umax;
            tests_for!(a, b, e; S1 S2 S3 S4 S5 S6);
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
                assert_eq!(r.add(&am, &bm), $a.addm($b, &P));
                assert_eq!(r.sub(&am, &bm), $a.subm($b, &P));
                assert_eq!(r.mul(&am, &bm), $a.mulm($b, &P));
                assert_eq!(r.neg(am), $a.negm(&P));
                assert_eq!(r.inv(am), $a.invm(&P));
                assert_eq!(r.dbl(am), $a.dblm(&P));
                assert_eq!(r.sqr(am), $a.sqm(&P));
                assert_eq!(r.pow(am, &$e), $a.powm($e, &P));
            })*);
        }

        for _ in 0..NRANDOM {
            let a = random::<u64>();
            let b = random::<u64>();
            let e = random::<u8>() as u64;
            tests_for!(a, b, e; S64_1 S64_2 S64_3 S64_4);
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
                assert_eq!(r.add(&am, &bm), $a.addm($b, &P));
                assert_eq!(r.sub(&am, &bm), $a.subm($b, &P));
                assert_eq!(r.mul(&am, &bm), $a.mulm($b, &P));
                assert_eq!(r.neg(am), $a.negm(&P));
                assert_eq!(r.inv(am), $a.invm(&P));
                assert_eq!(r.dbl(am), $a.dblm(&P));
                assert_eq!(r.sqr(am), $a.sqm(&P));
                assert_eq!(r.pow(am, &$e), $a.powm($e, &P));
            })*);
        }

        for _ in 0..NRANDOM {
            let a = random::<u32>();
            let b = random::<u32>();
            let e = random::<u8>() as u32;
            tests_for!(a, b, e; S32_1 S32_2 S32_3 S32_4);
        }
    }

    #[test]
    fn test_add_near_overflow_u64() {
        // 2^64 - 2^32 + 1 = 0xFFFFFFFF00000001, near u64::MAX
        type S = FixedTrinomialSolinas64<64, 32, 1>;
        const P: u64 = <S>::MODULUS;
        assert_eq!(P, 0xFFFFFFFF00000001);
        let r = S::new(&P);
        // Values near P-1; their sum exceeds u64::MAX
        // (P-1) + (P-2) = 2P-3 ≡ P-3 (mod P)
        let a = r.transform(P - 1);
        let b = r.transform(P - 2);
        assert_eq!(r.residue(r.add(&a, &b)), P - 3);
        // dbl near overflow: 2*(P-1) = 2P-2 ≡ P-2 (mod P)
        let c = r.transform(P - 1);
        assert_eq!(r.residue(r.dbl(c)), P - 2);
    }
}
