use crate::impl_fixed_monty_ops;
use crate::reduced::{impl_reduced_binary_pow, impl_reduced_ops};
use crate::{powm_u32, powm_u64, udouble, umax, ModularUnaryOps, Reducer};

// Proth primes: m = K * 2^N + 1 (K odd, K < 2^N)
//
// Montgomery REDC with R = 2^BITS (so R > m always).  The Montgomery
// constant N0 = -m⁻¹ mod R is computed at compile time via Newton
// iteration.  Because R > m the REDC result is always < 2·m, so a
// single conditional subtraction normalises.
//
// The product m·p inside REDC is expanded using the Proth form:
//   m·(K·2^N + 1) = (m·K)<<N + m
// Since K ≤ 255 (u8), the (m·K) multiply is narrow (at most 8 bits),
// replacing a full-width multiply-add with a shift and an add.

// --- macro for FixedProth32 / FixedProth64 ----------------------------------

/// Debug-only primality heuristic: checks that `m` is not divisible by
/// small primes (3, 5, 7, 11, 13), allowing `m` to be the small prime itself.
macro_rules! debug_assert_prime_candidate {
    ($m:expr) => {
        debug_assert!(
            ($m == 3 || $m % 3 != 0)
                && ($m == 5 || $m % 5 != 0)
                && ($m == 7 || $m % 7 != 0)
                && ($m == 11 || $m % 11 != 0)
                && ($m == 13 || $m % 13 != 0)
        )
    };
}

macro_rules! impl_fixed_proth_inherent {
    ($TypeName:ident, $T:ty, $D:ty, $neginv_fn:path, $powm:ident) => {
        impl<const N: u8, const K: u8> $TypeName<N, K> {
            /// Compile-time guard: N must be strictly less than the type bit-width.
            const _N_BOUND_CHECK: () = assert!((N as u32) < <$T>::BITS);

            pub const MODULUS: $T = {
                let p2n = match (1 as $T).checked_shl(N as u32) {
                    Some(v) => v,
                    None => unreachable!(),
                };
                let m = (K as $T).wrapping_mul(p2n).wrapping_add(1);
                // MODULUS ≤ φ·R guarantees `reduce` never overflows the double-word
                // sum (φ = (√5−1)/2 ≈ 0.618).
                assert!(
                    m as u128
                        <= match <$T>::BITS {
                            32 => 2_654_435_769u128,
                            64 => 11_400_714_819_323_199_485u128,
                            _ => unreachable!(),
                        },
                    "MODULUS exceeds overflow-free bound; lower N or use FixedMontgomery"
                );
                m
            };

            /// Montgomery constant:  -MODULUS⁻¹ mod 2^BITS
            const N0: $T = $neginv_fn(Self::MODULUS);

            /// R² mod MODULUS  (R = 2^BITS, so R² = 2^{2·BITS})
            const R2: $T = $powm(2, (2 * <$T>::BITS) as $T, Self::MODULUS);

            #[inline]
            pub fn reduce(&self, t: $D) -> $T {
                // Standard Montgomery REDC with Proth-optimised m·p product.
                // MODULUS ≤ φ·R (guaranteed at compile time) ensures the sum
                // t + m·MODULUS never exceeds the double-word width.
                let m = (t as $T).wrapping_mul(Self::N0);
                // m·p = m·(K·2^N + 1) = (m·K)<<N + m
                let mp = ((m as $D) * (K as $D)) << N;
                let mp = mp.wrapping_add(m as $D);
                let r = (t.wrapping_add(mp) >> <$T>::BITS) as $T;
                if r >= Self::MODULUS {
                    r - Self::MODULUS
                } else {
                    r
                }
            }
        }
    };
}

/// A modular reducer for Proth primes `K * 2^N + 1` with 32-bit operands.
///
/// Supports `N` up to 31, `K` odd with `K < 2^N`.  Montgomery form with `R = 2³²`.
///
/// # Example
///
/// ```rust
/// use num_modular::{FixedProth32, Reducer};
///
/// const N: u8 = 4;
/// const K: u8 = 1;
/// let modulus = (K as u32) * (1u32 << N) + 1; // 1*2^4 + 1 = 17
/// let reducer = FixedProth32::<N, K>::new(&modulus);
/// let a = reducer.transform(3);
/// let b = reducer.transform(5);
/// assert_eq!(reducer.residue(reducer.add(&a, &b)), 8);
/// assert_eq!(reducer.residue(reducer.mul(&a, &b)), 15);
/// ```
#[must_use]
#[derive(Debug, Clone, Copy)]
pub struct FixedProth32<const N: u8, const K: u8>;

impl_fixed_proth_inherent!(
    FixedProth32,
    u32,
    u64,
    crate::monty::neg_mod_inv::u32::neginv,
    powm_u32
);

impl<const N: u8, const K: u8> Reducer<u32> for FixedProth32<N, K> {
    #[inline]
    fn new(m: &u32) -> Self {
        assert!(
            *m == Self::MODULUS,
            "the given modulus doesn't match with the generic params"
        );
        assert!(N < 32, "N must be less than type bit width");
        assert!(N > 0, "N must be positive");
        assert!(K > 0, "K must be positive");
        assert!(K % 2 == 1, "K must be odd");
        assert!(
            (K as u64) * (1_u64 << (N as u32)) < u32::MAX as u64,
            "K·2^N + 1 exceeds type maximum"
        );
        debug_assert!((K as u32) < (1u32 << (N as u32)), "K must be less than 2^N");
        debug_assert_prime_candidate!(Self::MODULUS);
        Self {}
    }
    impl_fixed_monty_ops!(u32, u64, Self::R2, primitive);
}

/// A modular reducer for Proth primes `K * 2^N + 1` with 64-bit operands.
///
/// Supports `N` up to 63, `K` odd with `K < 2^N`.  Montgomery form with `R = 2⁶⁴`.
///
/// # Example
///
/// ```rust
/// use num_modular::{FixedProth64, Reducer};
///
/// const N: u8 = 5;
/// const K: u8 = 3;
/// let modulus = (K as u64) * (1u64 << N) + 1; // 3*2^5 + 1 = 97
/// let reducer = FixedProth64::<N, K>::new(&modulus);
/// let a = reducer.transform(10);
/// let b = reducer.transform(20);
/// assert_eq!(reducer.residue(reducer.mul(&a, &b)), (10u64 * 20) % 97);
/// ```
#[must_use]
#[derive(Debug, Clone, Copy)]
pub struct FixedProth64<const N: u8, const K: u8>;

impl_fixed_proth_inherent!(
    FixedProth64,
    u64,
    u128,
    crate::monty::neg_mod_inv::u64::neginv,
    powm_u64
);

impl<const N: u8, const K: u8> Reducer<u64> for FixedProth64<N, K> {
    #[inline]
    fn new(m: &u64) -> Self {
        assert!(
            *m == Self::MODULUS,
            "the given modulus doesn't match with the generic params"
        );
        assert!(N < 64, "N must be less than type bit width");
        assert!(N > 0, "N must be positive");
        assert!(K > 0, "K must be positive");
        assert!(K % 2 == 1, "K must be odd");
        assert!(
            (K as u128) * (1_u128 << (N as u32)) < u64::MAX as u128,
            "K·2^N + 1 exceeds type maximum"
        );
        debug_assert!((K as u64) < (1u64 << (N as u32)), "K must be less than 2^N");
        debug_assert_prime_candidate!(Self::MODULUS);
        Self {}
    }
    impl_fixed_monty_ops!(u64, u128, Self::R2, primitive);
}

// ── FixedProth (umax / udouble) ──────────────────────────────────────────────

/// A modular reducer for Proth primes `K * 2^N + 1`.
///
/// Supports `N` up to 127, `K` odd with `K < 2^N`.  Montgomery form with `R = 2¹²⁸`.
///
/// # Example
///
/// ```rust
/// use num_modular::{FixedProth, Reducer};
///
/// const N: u8 = 16;
/// const K: u8 = 1;
/// let modulus = (K as u128) * (1u128 << N) + 1; // 2^16 + 1 = 65537
/// let reducer = FixedProth::<N, K>::new(&modulus);
/// let a = reducer.transform(1000);
/// let b = reducer.transform(2000);
/// assert_eq!(reducer.residue(reducer.mul(&a, &b)), (1000u128 * 2000) % modulus);
/// ```
#[must_use]
#[derive(Debug, Clone, Copy)]
pub struct FixedProth<const N: u8, const K: u8>;

impl<const N: u8, const K: u8> FixedProth<N, K> {
    /// Compile-time guard: N must be strictly less than 128.
    const _N_BOUND_CHECK_U128: () = assert!(N < 128);

    pub const MODULUS: umax = {
        let p2n = match 1u128.checked_shl(N as u32) {
            Some(v) => v,
            None => unreachable!(),
        };
        let m = (K as u128).wrapping_mul(p2n).wrapping_add(1);
        // MODULUS ≤ φ·R guarantees `reduce` never overflows the udouble sum
        // (φ = (√5−1)/2 ≈ 0.618).
        assert!(
            m <= 210_306_068_529_402_891_650_266_558_847_000_772_608,
            "MODULUS exceeds overflow-free bound; lower N or use FixedMontgomery"
        );
        m
    };

    /// Montgomery constant:  -MODULUS⁻¹ mod 2¹²⁸
    const N0: umax = crate::monty::neg_mod_inv::u128::neginv(Self::MODULUS);

    /// R² mod MODULUS  (R = 2¹²⁸, so R² = 2²⁵⁶)
    const R2: umax = {
        let r = udouble { hi: 1, lo: 0 }.div_rem_2by1(Self::MODULUS).1; // 2¹²⁸ mod MODULUS
        udouble::widening_square(r).div_rem_2by1(Self::MODULUS).1 // 2²⁵⁶ mod MODULUS
    };

    /// Montgomery REDC with R = 2¹²⁸ and Proth-optimised m·p product.
    #[must_use]
    #[inline]
    pub fn reduce(&self, t: udouble) -> umax {
        let m = t.lo.wrapping_mul(Self::N0);
        // m·p = m·(K·2^N + 1) = (m·K)<<N + m
        // K ≤ 255, so the widening_mul is narrow (at most 8 bits).
        let mk = udouble::widening_mul(m, K as u128);
        let mp = mk.shl_u32(N as u32) + udouble { hi: 0, lo: m };
        let r = (t + mp).hi;
        if r >= Self::MODULUS {
            r - Self::MODULUS
        } else {
            r
        }
    }
}

impl<const N: u8, const K: u8> Reducer<umax> for FixedProth<N, K> {
    #[inline]
    fn new(m: &umax) -> Self {
        assert!(
            *m == Self::MODULUS,
            "the given modulus doesn't match with the generic params"
        );
        assert!(N < 128, "N must be less than type bit width");
        assert!(N > 0, "N must be positive");
        assert!(K > 0, "K must be positive");
        assert!(K % 2 == 1, "K must be odd");
        assert!(
            (K as u128) * (1u128 << (N as u32)) < u128::MAX,
            "K·2^N + 1 exceeds type maximum"
        );
        debug_assert!(
            (K as u128) < (1u128 << (N as u32)),
            "K must be less than 2^N"
        );
        debug_assert_prime_candidate!(Self::MODULUS);
        Self {}
    }
    impl_fixed_monty_ops!(umax, udouble, Self::R2, udouble);
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{ModularCoreOps, ModularPow};
    use rand::random;

    // u128 types
    type P128_1 = FixedProth<2, 1>; // m = 5
    type P128_2 = FixedProth<4, 1>; // m = 17
    type P128_3 = FixedProth<5, 3>; // m = 97
    type P128_4 = FixedProth<8, 3>; // m = 769
    type P128_5 = FixedProth<16, 1>; // m = 65537

    // u64 types
    type P64_1 = FixedProth64<4, 1>; // m = 17
    type P64_2 = FixedProth64<5, 3>; // m = 97
    type P64_3 = FixedProth64<8, 1>; // m = 257
    type P64_4 = FixedProth64<16, 1>; // m = 65537

    // u32 types
    type P32_1 = FixedProth32<2, 1>; // m = 5
    type P32_2 = FixedProth32<2, 3>; // m = 13
    type P32_3 = FixedProth32<4, 1>; // m = 17
    type P32_4 = FixedProth32<3, 5>; // m = 41

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

    /// Reduce correctness with MODULUS near the overflow-free bound.
    #[test]
    fn test_reduce_near_bound() {
        // 255·2^23 + 1 = 2,139,095,041 (close to φ·2^32 threshold 2,654,435,769)
        type S = FixedProth32<23, 255>;
        const M: u32 = <S>::MODULUS;
        let r = S::new(&M);

        for _ in 0..10 {
            let a = random::<u32>() % M;
            let b = random::<u32>() % M;
            let am = r.transform(a);
            let bm = r.transform(b);
            let result = r.residue(r.mul(&am, &bm));
            assert_eq!(result, a.mulm(b, &M));
        }
    }

    /// inv with MODULUS > usize::MAX should not truncate.
    #[test]
    fn test_inv_no_truncation_u128() {
        // N=60 < 64 but MODULUS = 31·2^60+1 > u64::MAX, so the old
        // `N < usize::BITS` gate would incorrectly take the usize path.
        type S = FixedProth<60, 31>;
        const M: u128 = <S>::MODULUS;
        assert!(
            M > u64::MAX as u128,
            "MODULUS must exceed usize for this test"
        );
        let r = S::new(&M);

        let a: u128 = 1234567890123456789 % M;
        let a_mont = r.transform(a);
        let inv = r.inv(a_mont).expect("inv should succeed");
        let result = r.residue(inv);
        assert_eq!(result.mulm(a, &M), 1u128, "inv truncation bug");
    }

    /// K·2^N exceeding type max should panic, not silently wrap.
    #[test]
    #[should_panic(expected = "exceeds type maximum")]
    fn test_modulus_overflow_panics_u32() {
        type S = FixedProth32<31, 3>; // 3·2^31+1 > 2^32
        const M: u32 = <S>::MODULUS; // wraps to 2^31+1
        let _ = S::new(&M); // should panic
    }

    /// FixedProth with N>64 should compute reduce correctly
    /// (no shift truncation in the Proth-optimised m·p product).
    #[test]
    fn test_reduce_n_gt_64() {
        type S = FixedProth<65, 3>; // MODULUS = 3·2^65 + 1
        const M: u128 = <S>::MODULUS;
        let r = S::new(&M);

        for _ in 0..10 {
            let a = random::<u128>() % M;
            let b = random::<u128>() % M;
            let am = r.transform(a);
            let bm = r.transform(b);
            let result = r.residue(r.mul(&am, &bm));
            assert_eq!(result, a.mulm(b, &M), "shift truncation bug for N>64");
        }
    }
}
