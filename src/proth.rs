use crate::reduced::impl_reduced_binary_pow;
use crate::{powm_u32, powm_u64, udouble, umax, ModularUnaryOps, Reducer, Vanilla};

// Proth primes: m = K * 2^N + 1 (K odd, K < 2^N)
//
// Montgomery REDC with R = 2^BITS (so R > m always).  The Montgomery
// constant N0 = -m⁻¹ mod R is computed at compile time via Newton
// iteration.  Because R > m the REDC result is always < 2·m, so a
// single conditional subtraction normalises.
//
// The product m·p inside REDC is expanded using the Proth form:
//   m·(K·2^N + 1) = (m·K)<<N + m
// which replaces a full-width multiply-add with a narrow multiply
// (K is small), a shift, and an add.

// --- const helpers (not pub) ------------------------------------------------

/// Compute x⁻¹ mod 2³² using Newton iteration (quadratic convergence).
const fn mod_inv_pow2_u32(x: u32) -> u32 {
    // x must be odd (gcd(x,2)=1)
    let mut inv = 1u32; // correct mod 2¹
    inv = inv.wrapping_mul(2u32.wrapping_sub(x.wrapping_mul(inv))); // mod 2²
    inv = inv.wrapping_mul(2u32.wrapping_sub(x.wrapping_mul(inv))); // mod 2⁴
    inv = inv.wrapping_mul(2u32.wrapping_sub(x.wrapping_mul(inv))); // mod 2⁸
    inv = inv.wrapping_mul(2u32.wrapping_sub(x.wrapping_mul(inv))); // mod 2¹⁶
    inv = inv.wrapping_mul(2u32.wrapping_sub(x.wrapping_mul(inv))); // mod 2³²
    inv
}

/// Compute x⁻¹ mod 2⁶⁴ using Newton iteration.
const fn mod_inv_pow2_u64(x: u64) -> u64 {
    let mut inv = 1u64;
    inv = inv.wrapping_mul(2u64.wrapping_sub(x.wrapping_mul(inv))); // mod 2²
    inv = inv.wrapping_mul(2u64.wrapping_sub(x.wrapping_mul(inv))); // mod 2⁴
    inv = inv.wrapping_mul(2u64.wrapping_sub(x.wrapping_mul(inv))); // mod 2⁸
    inv = inv.wrapping_mul(2u64.wrapping_sub(x.wrapping_mul(inv))); // mod 2¹⁶
    inv = inv.wrapping_mul(2u64.wrapping_sub(x.wrapping_mul(inv))); // mod 2³²
    inv = inv.wrapping_mul(2u64.wrapping_sub(x.wrapping_mul(inv))); // mod 2⁶⁴
    inv
}

// --- macro for FixedProth32 / FixedProth64 ----------------------------------

macro_rules! impl_fixed_proth_primitive {
    ($TypeName:ident, $T:ty, $D:ty, $max_N:expr, $inv_pow2:ident, $powm:ident) => {
        impl<const N: u8, const K: $T> $TypeName<N, K> {
            pub const MODULUS: $T = {
                let p2n = match (1 as $T).checked_shl(N as u32) {
                    Some(v) => v,
                    None => 0,
                };
                K.wrapping_mul(p2n).wrapping_add(1)
            };

            /// Montgomery constant:  -MODULUS⁻¹ mod 2^BITS
            const N0: $T = (0 as $T).wrapping_sub($inv_pow2(Self::MODULUS));

            /// R² mod MODULUS  (R = 2^BITS, so R² = 2^{2·BITS})
            const R2: $T = $powm(2, (2 * <$T>::BITS) as $T, Self::MODULUS);

            #[inline]
            pub fn reduce_single(&self, v: $T) -> $T {
                self.reduce_double(v as $D)
            }

            pub fn reduce_double(&self, t: $D) -> $T {
                // Standard Montgomery REDC with Proth-optimised m·p product.
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
                Self {}
            }
            #[inline]
            fn transform(&self, target: $T) -> $T {
                if target == 0 {
                    return 0;
                }
                self.reduce_double((target as $D) * (Self::R2 as $D))
            }
            #[inline]
            fn check(&self, target: &$T) -> bool { *target < Self::MODULUS }
            #[inline]
            fn residue(&self, target: $T) -> $T {
                if target == 0 {
                    return 0;
                }
                self.reduce_single(target)
            }
            #[inline]
            fn modulus(&self) -> $T { Self::MODULUS }
            #[inline]
            fn is_zero(&self, target: &$T) -> bool { target == &0 }

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
                self.reduce_double((*lhs as $D) * (*rhs as $D))
            }
            #[inline]
            fn inv(&self, target: $T) -> Option<$T> {
                let plain = if target == 0 { 0 } else { self.reduce_single(target) };
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
                Some(self.reduce_double((inv_plain as $D) * (Self::R2 as $D)))
            }
            #[inline]
            fn sqr(&self, target: $T) -> $T {
                self.reduce_double((target as $D) * (target as $D))
            }

            impl_reduced_binary_pow!($T);
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
/// const K: u32 = 1;
/// let modulus = K * (1u32 << N) + 1; // 1*2^4 + 1 = 17
/// let reducer = FixedProth32::<N, K>::new(&modulus);
/// let a = reducer.transform(3);
/// let b = reducer.transform(5);
/// assert_eq!(reducer.residue(reducer.add(&a, &b)), 8);
/// assert_eq!(reducer.residue(reducer.mul(&a, &b)), 15);
/// ```
#[derive(Debug, Clone, Copy)]
pub struct FixedProth32<const N: u8, const K: u32>;

impl_fixed_proth_primitive!(
    FixedProth32, u32, u64, 31,
    mod_inv_pow2_u32, powm_u32
);

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
/// const K: u64 = 3;
/// let modulus = K * (1u64 << N) + 1; // 3*2^5 + 1 = 97
/// let reducer = FixedProth64::<N, K>::new(&modulus);
/// let a = reducer.transform(10);
/// let b = reducer.transform(20);
/// assert_eq!(reducer.residue(reducer.mul(&a, &b)), (10u64 * 20) % 97);
/// ```
#[derive(Debug, Clone, Copy)]
pub struct FixedProth64<const N: u8, const K: u64>;

impl_fixed_proth_primitive!(
    FixedProth64, u64, u128, 63,
    mod_inv_pow2_u64, powm_u64
);

// ── FixedProth (umax / udouble) ──────────────────────────────────────────────

/// A modular reducer for Proth primes `K * 2^N + 1`.
///
/// Supports `N` up to 127, `K` odd with `K < 2^N`.  Montgomery form with `R = 2^N`.
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

impl<const N: u8, const K: umax> FixedProth<N, K> {
    const BITMASK: umax = match 1u128.checked_shl(N as u32) {
        Some(v) => v.wrapping_sub(1),
        None => umax::MAX,
    };
    pub const MODULUS: umax = {
        let p2n = match 1u128.checked_shl(N as u32) {
            Some(v) => v,
            None => 0,
        };
        K.wrapping_mul(p2n).wrapping_add(1)
    };

    #[inline]
    pub fn reduce_single(&self, v: umax) -> umax {
        Self::reduce_double(udouble { hi: 0, lo: v })
    }

    /// REDC with R = 2^N.  Uses the Proth identity m ≡ 1 mod R.
    pub fn reduce_double(v: udouble) -> umax {
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
        if acc.hi > 0 || acc.lo >= Self::MODULUS {
            acc % Self::MODULUS
        } else {
            acc.lo
        }
    }
}

impl<const N: u8, const K: umax> Reducer<umax> for FixedProth<N, K> {
    #[inline]
    fn new(m: &umax) -> Self {
        assert!(
            *m == Self::MODULUS,
            "the given modulus doesn't match with the generic params"
        );
        debug_assert!(N <= 127);
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
        let r = udouble { hi: 0, lo: 1 } << N;
        let r2 = udouble::widening_square(r.lo);
        Self {
            r2: r2 % Self::MODULUS,
        }
    }
    #[inline]
    fn transform(&self, target: umax) -> umax {
        if target == 0 {
            return 0;
        }
        Self::reduce_double(udouble::widening_mul(target, self.r2))
    }
    #[inline]
    fn check(&self, target: &umax) -> bool { *target < Self::MODULUS }
    #[inline]
    fn residue(&self, target: umax) -> umax {
        if target == 0 {
            return 0;
        }
        self.reduce_single(target)
    }
    #[inline]
    fn modulus(&self) -> umax { Self::MODULUS }
    #[inline]
    fn is_zero(&self, target: &umax) -> bool { target == &0 }

    #[inline]
    fn add(&self, lhs: &umax, rhs: &umax) -> umax {
        Vanilla::<umax>::add(&Self::MODULUS, *lhs, *rhs)
    }
    #[inline]
    fn sub(&self, lhs: &umax, rhs: &umax) -> umax {
        Vanilla::<umax>::sub(&Self::MODULUS, *lhs, *rhs)
    }
    #[inline]
    fn dbl(&self, target: umax) -> umax {
        Vanilla::<umax>::dbl(&Self::MODULUS, target)
    }
    #[inline]
    fn neg(&self, target: umax) -> umax {
        Vanilla::<umax>::neg(&Self::MODULUS, target)
    }
    #[inline]
    fn mul(&self, lhs: &umax, rhs: &umax) -> umax {
        Self::reduce_double(udouble::widening_mul(*lhs, *rhs))
    }
    #[inline]
    fn inv(&self, target: umax) -> Option<umax> {
        let plain = if target == 0 { 0 } else { self.reduce_single(target) };
        let inv_plain = if (N as u32) < usize::BITS {
            (plain as usize)
                .invm(&(Self::MODULUS as usize))
                .map(|v| v as umax)
        } else {
            plain.invm(&Self::MODULUS)
        }?;
        if inv_plain == 0 {
            return Some(0);
        }
        Some(Self::reduce_double(udouble::widening_mul(inv_plain, self.r2)))
    }
    #[inline]
    fn sqr(&self, target: umax) -> umax {
        Self::reduce_double(udouble::widening_square(target))
    }

    impl_reduced_binary_pow!(umax);
}

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
