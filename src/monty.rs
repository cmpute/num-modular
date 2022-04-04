use crate::{udouble, ModularInteger};
use core::ops::{Add, Mul, Neg, Sub};
use num_integer::Integer;
use num_traits::Pow;

#[cfg(std)]
use std::rc::Rc;

/// Operations of a integer represented in Montgomery form. This data type can
/// be used in place of a normal integer with regard to modular arithmetics.
///
/// The generic type T represents the underlying integer representation, and
/// R=2^B will be used as the auxiliary modulus, where B is automatically selected
/// based on the size of T.
pub trait Montgomery: Sized {
    /// The type for inversion of the modulus.
    ///
    /// This type is usually the same as Self, but it can be smaller when using
    /// Montgomery form on multi-precision integer representations.
    type Inv;

    /// The type of integer with double width. It is only used in `reduce()`,
    /// so it's okay that it's not actually doubled with
    type Double;

    /// Calculate -(m^-1) mod R
    fn neginv(m: &Self) -> Self::Inv; // TODO: return Option<Self::Inv>, return None if m is even

    /// Transform a normal integer into Montgomery form (compute `target*R mod m`)
    fn transform(target: Self, m: &Self) -> Self;

    /// Transform a montgomery form back to normal integer (compute `monty/R mod m`)
    fn reduce(monty: Self::Double, m: &Self, minv: &Self::Inv) -> Self;

    /// Calculate (lhs + rhs) mod m in Montgomery form
    fn add(lhs: &Self, rhs: &Self, m: &Self) -> Self;

    /// Calculate (lhs - rhs) mod m in Montgomery form
    fn sub(lhs: &Self, rhs: &Self, m: &Self) -> Self;

    /// Calculate -monty mod m in Montgomery form
    fn neg(monty: &Self, m: &Self) -> Self;

    /// Calculate (lhs * rhs) mod m in Montgomery form
    fn mul(lhs: &Self, rhs: &Self, m: &Self, minv: &Self::Inv) -> Self;

    /// Calculate base ^ exp mod m in Montgomery form
    fn pow(base: &Self, exp: &Self, m: &Self, minv: &Self::Inv) -> Self;
}

// Entry i contains (2i+1)^(-1) mod 256.
const BINVERT_TABLE: [u8; 128] = [
    0x01, 0xAB, 0xCD, 0xB7, 0x39, 0xA3, 0xC5, 0xEF, 0xF1, 0x1B, 0x3D, 0xA7, 0x29, 0x13, 0x35, 0xDF,
    0xE1, 0x8B, 0xAD, 0x97, 0x19, 0x83, 0xA5, 0xCF, 0xD1, 0xFB, 0x1D, 0x87, 0x09, 0xF3, 0x15, 0xBF,
    0xC1, 0x6B, 0x8D, 0x77, 0xF9, 0x63, 0x85, 0xAF, 0xB1, 0xDB, 0xFD, 0x67, 0xE9, 0xD3, 0xF5, 0x9F,
    0xA1, 0x4B, 0x6D, 0x57, 0xD9, 0x43, 0x65, 0x8F, 0x91, 0xBB, 0xDD, 0x47, 0xC9, 0xB3, 0xD5, 0x7F,
    0x81, 0x2B, 0x4D, 0x37, 0xB9, 0x23, 0x45, 0x6F, 0x71, 0x9B, 0xBD, 0x27, 0xA9, 0x93, 0xB5, 0x5F,
    0x61, 0x0B, 0x2D, 0x17, 0x99, 0x03, 0x25, 0x4F, 0x51, 0x7B, 0x9D, 0x07, 0x89, 0x73, 0x95, 0x3F,
    0x41, 0xEB, 0x0D, 0xF7, 0x79, 0xE3, 0x05, 0x2F, 0x31, 0x5B, 0x7D, 0xE7, 0x69, 0x53, 0x75, 0x1F,
    0x21, 0xCB, 0xED, 0xD7, 0x59, 0xC3, 0xE5, 0x0F, 0x11, 0x3B, 0x5D, 0xC7, 0x49, 0x33, 0x55, 0xFF,
];

macro_rules! impl_uprim_montgomery_core {
    () => {
        #[inline]
        fn transform(target: Self, m: &Self) -> Self {
            if target == 0 {
                return 0;
            }
            (((target as Self::Double) << Self::BITS) % (*m as Self::Double)) as _
        }

        // REDC algorithm
        fn reduce(monty: Self::Double, m: &Self, minv: &Self::Inv) -> Self {
            debug_assert!(monty < ((*m as Self::Double) << Self::BITS));

            let tm = (monty as Self).wrapping_mul(*minv);
            let (t, overflow) = monty.overflowing_add((tm as Self::Double) * (*m as Self::Double));
            let t = (t >> Self::BITS) as Self;

            // in case of overflow, we need to add another `R mod m` = `R - m`
            let t = if overflow { t + m.wrapping_neg() } else { t };

            if &t >= m {
                return t - m;
            } else {
                return t;
            }
        }

        #[inline]
        fn mul(lhs: &Self, rhs: &Self, m: &Self, minv: &Self::Inv) -> Self {
            Montgomery::reduce((*lhs as Self::Double) * (*rhs as Self::Double), m, minv)
        }
    };
}
macro_rules! impl_uprim_montgomery {
    () => {
        #[inline]
        fn add(lhs: &Self, rhs: &Self, m: &Self) -> Self {
            let (sum, overflow) = lhs.overflowing_add(*rhs);
            if overflow {
                sum + m.wrapping_neg()
            } else {
                sum
            }
        }

        #[inline]
        fn sub(lhs: &Self, rhs: &Self, m: &Self) -> Self {
            if lhs >= rhs {
                lhs - rhs
            } else {
                m - (rhs - lhs)
            }
        }

        #[inline]
        fn neg(monty: &Self, m: &Self) -> Self {
            if monty == &0 {
                0
            } else {
                m - monty
            }
        }

        fn pow(base: &Self, exp: &Self, m: &Self, minv: &Self::Inv) -> Self {
            match *exp {
                1 => *base,
                2 => Montgomery::mul(base, base, m, minv),
                e => {
                    let mut multi = *base;
                    let mut exp = e;
                    let mut result = Montgomery::transform(1, m);
                    while exp > 0 {
                        if exp & 1 != 0 {
                            result = Montgomery::mul(&result, &multi, m, minv);
                        }
                        multi = Montgomery::mul(&multi, &multi, m, minv);
                        exp >>= 1;
                    }
                    result
                }
            }
        }
    };
}

impl Montgomery for u8 {
    type Inv = u8;
    type Double = u16;

    impl_uprim_montgomery_core!();
    impl_uprim_montgomery!();

    fn neginv(m: &Self) -> Self {
        BINVERT_TABLE[((m >> 1) & 0x7F) as usize].wrapping_neg()
    }
}

impl Montgomery for u16 {
    type Inv = u16;
    type Double = u32;

    impl_uprim_montgomery_core!();
    impl_uprim_montgomery!();

    fn neginv(m: &Self) -> Self {
        let i = BINVERT_TABLE[((m >> 1) & 0x7F) as usize] as u16;
        // Newton-Rhapson iteration (hensel lifting), see https://arxiv.org/abs/1303.0328
        i.wrapping_mul(*m).wrapping_sub(2).wrapping_mul(i)
    }
}

impl Montgomery for u32 {
    type Inv = u32;
    type Double = u64;

    impl_uprim_montgomery_core!();
    impl_uprim_montgomery!();

    fn neginv(m: &Self) -> Self {
        let i = BINVERT_TABLE[((m >> 1) & 0x7F) as usize] as u32;
        let i = 2u32.wrapping_sub(i.wrapping_mul(*m)).wrapping_mul(i);
        i.wrapping_mul(*m).wrapping_sub(2).wrapping_mul(i)
    }
}

impl Montgomery for u64 {
    type Inv = u64;
    type Double = u128;

    impl_uprim_montgomery_core!();
    impl_uprim_montgomery!();

    fn neginv(m: &Self) -> Self {
        let i = BINVERT_TABLE[((m >> 1) & 0x7F) as usize] as u64;
        let i = 2u64.wrapping_sub(i.wrapping_mul(*m)).wrapping_mul(i);
        let i = 2u64.wrapping_sub(i.wrapping_mul(*m)).wrapping_mul(i);
        i.wrapping_mul(*m).wrapping_sub(2).wrapping_mul(i)
    }
}

impl Montgomery for u128 {
    type Inv = u128;
    type Double = udouble;

    fn neginv(m: &Self) -> Self {
        let i = BINVERT_TABLE[((m >> 1) & 0x7F) as usize] as u128;
        let i = 2u128.wrapping_sub(i.wrapping_mul(*m)).wrapping_mul(i);
        let i = 2u128.wrapping_sub(i.wrapping_mul(*m)).wrapping_mul(i);
        let i = 2u128.wrapping_sub(i.wrapping_mul(*m)).wrapping_mul(i);
        i.wrapping_mul(*m).wrapping_sub(2).wrapping_mul(i)
    }

    #[inline]
    fn transform(target: Self, m: &Self) -> Self {
        if target == 0 {
            return 0;
        }
        let r = udouble { hi: target, lo: 0 } % udouble::from(*m);
        r.lo
    }

    // REDC algorithm
    fn reduce(monty: Self::Double, m: &Self, minv: &Self::Inv) -> Self {
        debug_assert!(monty < udouble { hi: *m, lo: 0 });

        let tm = monty.lo.wrapping_mul(*minv);
        let (t, overflow) = monty.overflowing_add(udouble::widening_mul(tm, *m));

        // in case of overflow, we need to add another `R mod m` = `R - m`
        let t = if overflow {
            t.hi + m.wrapping_neg()
        } else {
            t.hi
        };

        if &t >= m {
            return t - m;
        } else {
            return t;
        }
    }

    #[inline]
    fn mul(lhs: &Self, rhs: &Self, m: &Self, minv: &Self::Inv) -> Self {
        Montgomery::reduce(udouble::widening_mul(*lhs, *rhs), m, minv)
    }

    impl_uprim_montgomery!();
}

/// An integer represented in Montgomery form, it implements [ModularInteger] interface
/// and it's generally more efficient than the vanilla integer in modular operations.
#[derive(Debug, Clone, Copy)]
pub struct MontgomeryInt<T: Integer + Montgomery> {
    /// The Montgomery representation of the integer.
    a: T,

    /// The modulus.
    m: T,

    /// The negated modular inverse of the modulus mod R
    mi: T::Inv,
}

impl<T: Integer + Montgomery> MontgomeryInt<T> {
    #[inline]
    fn check_modulus_eq(&self, rhs: &Self) {
        if self.m != rhs.m {
            panic!("The modulus of two operators should be the same!");
        }
    }
}

impl<T: Integer + Montgomery> MontgomeryInt<T>
where
    T::Double: From<T>,
{
    /// Convert n into the modulo ring ℤ/mℤ (i.e. `n % m`)
    #[inline]
    pub fn new(n: T, m: T) -> Self {
        let minv = Montgomery::neginv(&m);
        let a = Montgomery::transform(n, &m);
        MontgomeryInt { a, m, mi: minv }
    }
}

impl<T: Integer + Montgomery> PartialEq for MontgomeryInt<T> {
    #[inline]
    fn eq(&self, other: &Self) -> bool {
        self.check_modulus_eq(other);
        self.a == other.a
    }
}

impl<T: Integer + Montgomery> Add for MontgomeryInt<T> {
    type Output = Self;

    #[inline]
    fn add(self, rhs: Self) -> Self::Output {
        self.check_modulus_eq(&rhs);
        let a = Montgomery::add(&self.a, &rhs.a, &self.m);
        MontgomeryInt {
            a,
            m: self.m,
            mi: self.mi,
        }
    }
}

impl<T: Integer + Montgomery> Sub for MontgomeryInt<T> {
    type Output = Self;

    #[inline]
    fn sub(self, rhs: Self) -> Self::Output {
        self.check_modulus_eq(&rhs);
        let a = Montgomery::sub(&self.a, &rhs.a, &self.m);
        MontgomeryInt {
            a,
            m: self.m,
            mi: self.mi,
        }
    }
}

impl<T: Integer + Montgomery> Neg for MontgomeryInt<T> {
    type Output = Self;

    #[inline]
    fn neg(self) -> Self::Output {
        let a = Montgomery::neg(&self.a, &self.m);
        MontgomeryInt {
            a,
            m: self.m,
            mi: self.mi,
        }
    }
}

impl<T: Integer + Montgomery> Mul for MontgomeryInt<T> {
    type Output = Self;

    #[inline]
    fn mul(self, rhs: Self) -> Self::Output {
        self.check_modulus_eq(&rhs);
        let a = Montgomery::mul(&self.a, &rhs.a, &self.m, &self.mi);
        MontgomeryInt {
            a,
            m: self.m,
            mi: self.mi,
        }
    }
}

impl<T: Integer + Montgomery> Pow<T> for MontgomeryInt<T> {
    type Output = Self;

    #[inline]
    fn pow(self, rhs: T) -> Self::Output {
        let a = Montgomery::pow(&self.a, &rhs, &self.m, &self.mi);
        MontgomeryInt {
            a,
            m: self.m,
            mi: self.mi,
        }
    }
}

impl<T: Integer + Montgomery + Clone> ModularInteger for MontgomeryInt<T>
where
    T::Double: From<T>,
    T::Inv: Clone,
{
    type Base = T;

    #[inline]
    fn modulus(&self) -> &T {
        &self.m
    }

    #[inline]
    fn residue(&self) -> T {
        Montgomery::reduce(T::Double::from(self.a.clone()), &self.m, &self.mi)
    }

    #[inline]
    fn new(&self, n: T) -> Self { // TODO(v0.3): rename to convert
        let a = Montgomery::transform(n, &self.m);
        MontgomeryInt {
            a,
            m: self.m.clone(),
            mi: self.mi.clone(),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ModularCoreOps;
    use rand::random;

    #[test]
    fn creation_test() {
        // a deterministic test case for u128
        let a = (0x81u128 << 120) - 1;
        let m = (0x81u128 << 119) - 1;
        let m = m >> m.trailing_zeros();
        assert_eq!(MontgomeryInt::new(a, m).residue(), a % m);

        // random creation test
        for _ in 0..10 {
            let a = random::<u8>();
            let m = random::<u8>() | 1;
            assert_eq!(MontgomeryInt::new(a, m).residue(), a % m);
    
            let a = random::<u16>();
            let m = random::<u16>() | 1;
            assert_eq!(MontgomeryInt::new(a, m).residue(), a % m);
    
            let a = random::<u32>();
            let m = random::<u32>() | 1;
            assert_eq!(MontgomeryInt::new(a, m).residue(), a % m);
    
            let a = random::<u64>();
            let m = random::<u64>() | 1;
            assert_eq!(MontgomeryInt::new(a, m).residue(), a % m);
    
            let a = random::<u128>();
            let m = random::<u128>() | 1;
            assert_eq!(MontgomeryInt::new(a, m).residue(), a % m);
        }
    }

    #[test]
    fn binary_op_tests() {
        // TODO(v0.3): test more operations
        for _ in 0..10 {
            let m = random::<u8>() | 1;
            let (a, b) = (random::<u8>(), random::<u8>());
            let am = MontgomeryInt::new(a, m);
            let bm = am.new(b);
            assert_eq!((am + bm).residue(), a.addm(b, &m));
            assert_eq!((am - bm).residue(), a.subm(b, &m));
            assert_eq!((am * bm).residue(), a.mulm(b, &m));

            let m = random::<u16>() | 1;
            let (a, b) = (random::<u16>(), random::<u16>());
            let am = MontgomeryInt::new(a, m);
            let bm = am.new(b);
            assert_eq!((am + bm).residue(), a.addm(b, &m));
            assert_eq!((am - bm).residue(), a.subm(b, &m));
            assert_eq!((am * bm).residue(), a.mulm(b, &m));
    
            let m = random::<u32>() | 1;
            let (a, b) = (random::<u32>(), random::<u32>());
            let am = MontgomeryInt::new(a, m);
            let bm = am.new(b);
            assert_eq!((am + bm).residue(), a.addm(b, &m));
            assert_eq!((am - bm).residue(), a.subm(b, &m));
            assert_eq!((am * bm).residue(), a.mulm(b, &m));

            let m = random::<u64>() | 1;
            let (a, b) = (random::<u64>(), random::<u64>());
            let am = MontgomeryInt::new(a, m);
            let bm = am.new(b);
            assert_eq!((am + bm).residue(), a.addm(b, &m));
            assert_eq!((am - bm).residue(), a.subm(b, &m));
            assert_eq!((am * bm).residue(), a.mulm(b, &m));
    
            let m = random::<u128>() | 1;
            let (a, b) = (random::<u128>(), random::<u128>());
            let am = MontgomeryInt::new(a, m);
            let bm = am.new(b);
            assert_eq!((am + bm).residue(), a.addm(b, &m));
            assert_eq!((am - bm).residue(), a.subm(b, &m));
            assert_eq!((am * bm).residue(), a.mulm(b, &m));
        }
    }
}
