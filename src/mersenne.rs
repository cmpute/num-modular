use crate::reduced::{impl_reduced_binary_pow, impl_reduced_ops};
use crate::{udouble, umax, ModularUnaryOps, Reducer};

macro_rules! impl_fixed_mersenne {
    (
        $TypeName:ident,
        $T:ty,
        $D:ty,
        $half_bits:expr,
        $max_P:expr,
        $kind:ident
    ) => {
        impl<const P: u8, const K: $T> $TypeName<P, K> {
            const BITMASK: $T = match (1 as $T).checked_shl(P as u32) {
                Some(v) => v.wrapping_sub(1),
                None => <$T>::MAX,
            };
            pub const MODULUS: $T = {
                let p1 = match (1 as $T).checked_shl(P as u32) {
                    Some(v) => v,
                    None => 0,
                };
                p1.wrapping_sub(K)
            };

            /// Worst-case fold count for `reduce_double`.
            /// Each fold replaces V = hi·2^P + lo with hi·K + lo (since 2^P ≡ K).
            /// For K = 1: always 2 folds (the carry chain terminates in at most one
            /// extra step). For K > 1: ⌈P/(P−⌈log₂K⌉)⌉ + 1 folds.
            const FOLDS: u32 = if K == 1 {
                2
            } else {
                let s = <$T>::BITS - K.leading_zeros(); // bit-width of K
                let gap = P as u32 - s;
                let folds_ceil = (P as u32 + gap - 1) / gap;
                folds_ceil + 1
            };

            /// Reduces a single-width value `v` modulo `2^P - K`.
            ///
            /// For the result of a widening multiplication or square, use
            /// [`reduce_double`](Self::reduce_double) instead.
            pub const fn reduce_single(v: $T) -> $T {
                let mut lo = v & Self::BITMASK;
                let mut hi = match v.checked_shr(P as u32) {
                    Some(s) => s,
                    None => 0,
                };
                while hi > 0 {
                    let sum = if K == 1 { hi + lo } else { hi * K + lo };
                    lo = sum & Self::BITMASK;
                    hi = match sum.checked_shr(P as u32) {
                        Some(s) => s,
                        None => 0,
                    };
                }
                if lo >= Self::MODULUS {
                    lo - Self::MODULUS
                } else {
                    lo
                }
            }

            impl_fixed_mersenne!(@reduce_double, $kind, $T, $D);
        }

        impl<const P: u8, const K: $T> Reducer<$T> for $TypeName<P, K> {
            #[inline]
            fn new(m: &$T) -> Self {
                assert!(
                    *m == Self::MODULUS,
                    "the given modulus doesn't match with the generic params"
                );
                debug_assert!(P <= $max_P);
                debug_assert!(K > 0 && K < (2 as $T).pow(P as u32 - 1) && K % 2 == 1);
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
                if (P as u32) < $half_bits {
                    Self::reduce_single(lhs * rhs)
                } else {
                    Self::reduce_double(impl_fixed_mersenne!(@widen_mul, $kind, $T, $D, lhs, rhs))
                }
            }
            #[inline]
            fn inv(&self, target: $T) -> Option<$T> {
                if (P as u32) < usize::BITS {
                    (target as usize)
                        .invm(&(Self::MODULUS as usize))
                        .map(|v| v as $T)
                } else {
                    target.invm(&Self::MODULUS)
                }
            }
            #[inline]
            fn sqr(&self, target: $T) -> $T {
                if (P as u32) < $half_bits {
                    Self::reduce_single(target * target)
                } else {
                    Self::reduce_double(impl_fixed_mersenne!(@widen_sqr, $kind, $T, $D, target))
                }
            }

            impl_reduced_binary_pow!($T);
        }
    };

    // Internal: reduce_double for primitive double-width types (u32→u64, u64→u128)
    //
    // For real pseudo-Mersennes, FOLDS is always ≤ 3 (K=1 → 2; small K → 3).
    // Unrolling replaces the data-dependent while loop with straight-line folds.
    // Extra folds past the true count are no-ops (hi reaches 0).
    (@reduce_double, primitive, $T:ty, $D:ty) => {
        /// Reduces a double-width value `v` modulo `2^P - K`.
        ///
        /// This handles widening-multiplication or widening-square results.
        /// For single-width values, use [`reduce_single`](Self::reduce_single).
        pub fn reduce_double(v: $D) -> $T {
            let mut lo = (v as $T) & Self::BITMASK;
            let mut hi = v >> P;
            macro_rules! mersenne_fold {
                () => {
                    let sum = if K == 1 {
                        hi + lo as $D
                    } else {
                        hi * (K as $D) + lo as $D
                    };
                    lo = (sum as $T) & Self::BITMASK;
                    hi = sum >> P;
                };
            }
            if Self::FOLDS <= 2 {
                #[allow(unused_assignments)] { mersenne_fold!(); }
                #[allow(unused_assignments)] { mersenne_fold!(); }
            } else if Self::FOLDS == 3 {
                #[allow(unused_assignments)] { mersenne_fold!(); }
                #[allow(unused_assignments)] { mersenne_fold!(); }
                #[allow(unused_assignments)] { mersenne_fold!(); }
            } else {
                while hi > 0 { mersenne_fold!(); }
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
    // Phase 1 (udouble while hi.hi > 0) is unreachable for valid P ≤ 128 since
    // hi = v >> P < 2^P ≤ 2^128 always fits in one word. Phase 2 uses u128
    // arithmetic and is unrolled when FOLDS ≤ 3 (all practical pseudo-Mersennes).
    (@reduce_double, udouble, $T:ty, $D:ty) => {
        /// Reduces a double-width value `v` modulo `2^P - K`.
        ///
        /// This handles widening-multiplication or widening-square results.
        /// For single-width values, use [`reduce_single`](Self::reduce_single).
        pub fn reduce_double(v: $D) -> $T {
            let mut lo = v.lo & Self::BITMASK;
            let mut hi = v >> P;
            while hi.hi > 0 {
                let sum = if K == 1 { hi + lo } else { hi * K + lo };
                lo = sum.lo & Self::BITMASK;
                hi = sum >> P;
            }
            let mut hi = hi.lo;
            macro_rules! mersenne_u128_fold {
                () => {
                    let sum = if K == 1 { hi + lo } else { hi * K + lo };
                    lo = sum & Self::BITMASK;
                    hi = match sum.checked_shr(P as u32) {
                        Some(s) => s,
                        None => 0,
                    };
                };
            }
            if Self::FOLDS <= 2 {
                #[allow(unused_assignments)] { mersenne_u128_fold!(); }
                #[allow(unused_assignments)] { mersenne_u128_fold!(); }
            } else if Self::FOLDS == 3 {
                #[allow(unused_assignments)] { mersenne_u128_fold!(); }
                #[allow(unused_assignments)] { mersenne_u128_fold!(); }
                #[allow(unused_assignments)] { mersenne_u128_fold!(); }
            } else {
                while hi > 0 { mersenne_u128_fold!(); }
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

/// A modular reducer for (pseudo) Mersenne numbers `2^P - K` as modulus with 32-bit operands.
///
/// Supports `P` up to 32 and `K < 2^(P-1)`. All inputs and outputs are `u32`.
/// The modulus `2^P - K` must be prime for modular inverse and Fermat-based operations to be valid.
///
/// # Example
///
/// ```rust
/// use num_modular::{FixedMersenne32, Reducer};
///
/// const P: u8 = 13;
/// const K: u32 = 1;
/// let modulus = (1u32 << P) - K; // 2^13 - 1 = 8191 (Mersenne prime)
/// let reducer = FixedMersenne32::<P, K>::new(&modulus);
/// let a = reducer.transform(100);
/// let b = reducer.transform(200);
/// assert_eq!(reducer.residue(reducer.add(&a, &b)), 300 % modulus);
/// ```
#[must_use]
#[derive(Debug, Clone, Copy)]
pub struct FixedMersenne32<const P: u8, const K: u32>();

impl_fixed_mersenne!(FixedMersenne32, u32, u64, 16, 32, primitive);

/// A modular reducer for (pseudo) Mersenne numbers `2^P - K` as modulus with 64-bit operands.
///
/// Supports `P` up to 64 and `K < 2^(P-1)`. All inputs and outputs are `u64`.
/// Uses `u128` as the double-width intermediate for multiplication and reduction.
/// The modulus `2^P - K` must be prime for modular inverse and Fermat-based operations to be valid.
///
/// # Example
///
/// ```rust
/// use num_modular::{FixedMersenne64, Reducer};
///
/// const P: u8 = 61;
/// const K: u64 = 1;
/// let modulus = (1u64 << P) - K; // 2^61 - 1 (Mersenne prime)
/// let reducer = FixedMersenne64::<P, K>::new(&modulus);
/// let a = reducer.transform(1000);
/// let b = reducer.transform(2000);
/// assert_eq!(reducer.residue(reducer.mul(&a, &b)), (1000u64 * 2000) % modulus);
/// ```
#[must_use]
#[derive(Debug, Clone, Copy)]
pub struct FixedMersenne64<const P: u8, const K: u64>();

impl_fixed_mersenne!(FixedMersenne64, u64, u128, 32, 64, primitive);

/// A modular reducer for (pseudo) Mersenne numbers `2^P - K` as modulus.
///
/// Supports `P` up to 128 and `K < 2^(P-1)`. All inputs and outputs are [umax] (currently `u128`).
///
/// The `P` is limited to 128 so that overflow checks aren't necessary. This covers all Mersenne
/// primes within the range of [umax] (i.e. `u128`).
///
/// # Example
///
/// ```rust
/// use num_modular::{FixedMersenne, Reducer, umax};
///
/// const P: u8 = 31;
/// const K: umax = 1;
/// let modulus = (1 << P) - K; // 2^31 - 1 (Mersenne prime)
/// let reducer = FixedMersenne::<P, K>::new(&modulus);
/// let a = reducer.transform(1000);
/// let b = reducer.transform(2000);
/// assert_eq!(reducer.residue(reducer.mul(&a, &b)), (1000 * 2000) % modulus);
/// ```
#[must_use]
#[derive(Debug, Clone, Copy)]
pub struct FixedMersenne<const P: u8, const K: umax>();

impl_fixed_mersenne!(FixedMersenne, umax, udouble, 64, 128, udouble);

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{ModularCoreOps, ModularPow};
    use rand::random;

    // u128 tests (existing)
    type M1 = FixedMersenne<31, 1>;
    type M2 = FixedMersenne<61, 1>;
    type M3 = FixedMersenne<127, 1>;
    type M4 = FixedMersenne<32, 5>;
    type M5 = FixedMersenne<56, 5>;
    type M6 = FixedMersenne<122, 3>;
    type M7 = FixedMersenne<128, 159>;

    // u64 tests
    type M64_1 = FixedMersenne64<31, 1>;
    type M64_2 = FixedMersenne64<61, 1>;
    type M64_3 = FixedMersenne64<32, 5>;
    type M64_4 = FixedMersenne64<64, 59>;

    // u32 tests
    type M32_1 = FixedMersenne32<13, 1>;
    type M32_2 = FixedMersenne32<31, 1>;
    type M32_3 = FixedMersenne32<16, 5>;

    const NRANDOM: u32 = 10;

    #[test]
    fn creation_test_u128() {
        const P: umax = (1 << 31) - 1;
        let m = M1::new(&P);
        assert_eq!(m.residue(m.transform(0)), 0);
        assert_eq!(m.residue(m.transform(1)), 1);
        assert_eq!(m.residue(m.transform(P)), 0);
        assert_eq!(m.residue(m.transform(P - 1)), P - 1);
        assert_eq!(m.residue(m.transform(P + 1)), 1);

        for _ in 0..NRANDOM {
            let a = random::<umax>();

            const P1: umax = (1 << 31) - 1;
            let m1 = M1::new(&P1);
            assert_eq!(m1.residue(m1.transform(a)), a % P1);
            const P2: umax = (1 << 61) - 1;
            let m2 = M2::new(&P2);
            assert_eq!(m2.residue(m2.transform(a)), a % P2);
            const P3: umax = (1 << 127) - 1;
            let m3 = M3::new(&P3);
            assert_eq!(m3.residue(m3.transform(a)), a % P3);
            const P4: umax = (1 << 32) - 5;
            let m4 = M4::new(&P4);
            assert_eq!(m4.residue(m4.transform(a)), a % P4);
            const P5: umax = (1 << 56) - 5;
            let m5 = M5::new(&P5);
            assert_eq!(m5.residue(m5.transform(a)), a % P5);
            const P6: umax = (1 << 122) - 3;
            let m6 = M6::new(&P6);
            assert_eq!(m6.residue(m6.transform(a)), a % P6);
            const P7: umax = M7::MODULUS;
            let m7 = M7::new(&P7);
            assert_eq!(m7.residue(m7.transform(a)), a % P7);
        }
    }

    #[test]
    fn creation_test_u64() {
        for _ in 0..NRANDOM {
            let a = random::<u64>();

            const P1: u64 = (1 << 31) - 1;
            let m1 = M64_1::new(&P1);
            assert_eq!(m1.residue(m1.transform(a)), a % P1);
            const P2: u64 = (1 << 61) - 1;
            let m2 = M64_2::new(&P2);
            assert_eq!(m2.residue(m2.transform(a)), a % P2);
            const P3: u64 = (1 << 32) - 5;
            let m3 = M64_3::new(&P3);
            assert_eq!(m3.residue(m3.transform(a)), a % P3);
            const P4: u64 = M64_4::MODULUS;
            let m4 = M64_4::new(&P4);
            assert_eq!(m4.residue(m4.transform(a)), a % P4);
        }
    }

    #[test]
    fn creation_test_u32() {
        for _ in 0..NRANDOM {
            let a = random::<u32>();

            const P1: u32 = (1 << 13) - 1;
            let m1 = M32_1::new(&P1);
            assert_eq!(m1.residue(m1.transform(a)), a % P1);
            const P2: u32 = (1 << 31) - 1;
            let m2 = M32_2::new(&P2);
            assert_eq!(m2.residue(m2.transform(a)), a % P2);
            const P3: u32 = (1 << 16) - 5;
            let m3 = M32_3::new(&P3);
            assert_eq!(m3.residue(m3.transform(a)), a % P3);
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
            tests_for!(a, b, e; M1 M2 M3 M4 M5 M6);
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
            tests_for!(a, b, e; M64_1 M64_2 M64_3);
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
            tests_for!(a, b, e; M32_1 M32_2 M32_3);
        }
    }

    #[test]
    fn small_prime_moduli() {
        // Small primes that are themselves 3, 5, 7, 11, or 13 should not
        // trigger the composite check debug_assert.
        // 2^2 - 1 = 3, 2^3 - 1 = 7, 2^4 - 1 = 15 (not prime, skip)
        // 2^5 - 1 = 31 (Mersenne prime)
        type M3 = FixedMersenne<2, 1>; // modulus = 3
        type M7 = FixedMersenne<3, 1>; // modulus = 7

        const M3_MOD: umax = M3::MODULUS;
        let r3 = M3::new(&M3_MOD);
        assert_eq!(r3.residue(r3.transform(5 % M3_MOD)), 5 % M3_MOD);

        const M7_MOD: umax = M7::MODULUS;
        let r7 = M7::new(&M7_MOD);
        assert_eq!(r7.residue(r7.transform(10 % M7_MOD)), 10 % M7_MOD);

        // 32-bit variants
        type M32_3 = FixedMersenne32<2, 1>; // modulus = 3
        type M32_7 = FixedMersenne32<3, 1>; // modulus = 7

        const P3: u32 = M32_3::MODULUS;
        let r3 = M32_3::new(&P3);
        assert_eq!(r3.residue(r3.transform(5 % P3)), 5 % P3);

        const P7: u32 = M32_7::MODULUS;
        let r7 = M32_7::new(&P7);
        assert_eq!(r7.residue(r7.transform(10 % P7)), 10 % P7);
    }
}
