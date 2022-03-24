use crate::ModularInteger;
use num_integer::Integer;
use num_traits::Pow;
use core::borrow::Borrow;
use core::ops::{Add, Mul, Neg, Sub};

#[cfg(std)]
use std::rc::Rc;

/// Operations of a integer represented in Montgomery form. This data type can
/// be used in place of a normal integer with regard to modular arithmetics.
///
/// The generic type T represents the underlying integer representation, and
/// R=2^B will be used as the auxiliary modulus, where B is automatically selected
/// based on the size of T.
// TODO: finalize API after we got the big integer Montgomery implemented
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
    fn neginv(m: &Self) -> Self::Inv;

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

macro_rules! impl_uprim_montgomery {
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

        #[inline]
        fn mul(lhs: &Self, rhs: &Self, m: &Self, minv: &Self::Inv) -> Self {
            Montgomery::reduce((*lhs as Self::Double) * (*rhs as Self::Double), m, minv)
        }

        fn pow(base: &Self, exp: &Self, m: &Self, minv: &Self::Inv) -> Self {
            match *exp {
                1 => *base,
                2 => Montgomery::mul(base, base, m, minv),
                e => {
                    let mut multi = *base;
                    let mut exp = e;
                    let mut result = Montgomery::transform(1, m); // TODO: subtract exp by 1 and use base?
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

    fn neginv(m: &Self) -> Self {
        BINVERT_TABLE[((m >> 1) & 0x7F) as usize].wrapping_neg()
    }

    impl_uprim_montgomery!();
}

impl Montgomery for u16 {
    type Inv = u16;
    type Double = u32;

    fn neginv(m: &Self) -> Self {
        let i = BINVERT_TABLE[((m >> 1) & 0x7F) as usize] as u16;
        // Newton-Rhapson iteration
        // See: https://arxiv.org/abs/1303.0328
        i.wrapping_mul(*m).wrapping_sub(2).wrapping_mul(i)
    }

    impl_uprim_montgomery!();
}

impl Montgomery for u32 {
    type Inv = u32;
    type Double = u64;

    fn neginv(m: &Self) -> Self {
        let i = BINVERT_TABLE[((m >> 1) & 0x7F) as usize] as u32;
        let i = 2u32.wrapping_sub(i.wrapping_mul(*m)).wrapping_mul(i);
        i.wrapping_mul(*m).wrapping_sub(2).wrapping_mul(i)
    }

    impl_uprim_montgomery!();
}

impl Montgomery for u64 {
    type Inv = u64;
    type Double = u128;

    fn neginv(m: &Self) -> Self {
        let i = BINVERT_TABLE[((m >> 1) & 0x7F) as usize] as u64;
        let i = 2u64.wrapping_sub(i.wrapping_mul(*m)).wrapping_mul(i);
        let i = 2u64.wrapping_sub(i.wrapping_mul(*m)).wrapping_mul(i);
        i.wrapping_mul(*m).wrapping_sub(2).wrapping_mul(i)
    }

    impl_uprim_montgomery!();
}

// XXX: implement Montgomery for u128 (double type is also u128), which requires efficient implementation of dual word mul_mod.
// REF: https://github.com/coreutils/coreutils/blob/master/src/factor.c (mulredc2)
// We can implement it efficiently with carrying_mul and widening_mul implemented (rust#85532)

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

/// A big integer represented in Montgomery form, it implements [ModularInteger] interface
/// and it's generally more efficient than the vanilla integer in modular operations.
///
/// The modulus is stored in heap to prevent frequent copying
#[derive(Debug, Clone)]
#[cfg(std)]
pub struct MontgomeryBigint<T: Integer + Montgomery> {
    /// The Montgomery representation of the integer.
    a: T,

    /// The modulus and its negated modular inverse.
    ///
    /// It's stored as a pointer to prevent frequent copying. It also allows
    /// quick checking of the equity of two moduli.
    minv: Rc<(T, T::Inv)>,
}

/// A word-size integer in Montgomery form with fixed modulus
// TODO: implement after we have const implementation of invm
#[derive(Debug, Clone, Copy)]
struct MontgomeryWord<const M: usize> (usize);

// XXX: we can also implement MontgomeryMersenne<const M: usize> to support Montgomery form
// with (Pseudo) Mersenne prime as modulo. REF: https://eprint.iacr.org/2018/1038.pdf

impl<T: Integer + Montgomery> MontgomeryInt<T> {
    #[inline]
    fn check_modulus_eq(&self, rhs: &Self) {
        if self.m != rhs.m {
            panic!("The modulus of two operators should be the same!");
        }
    }
}

#[cfg(std)]
impl<T: Integer + Montgomery> MontgomeryBigint<T> {
    #[inline]
    fn check_modulus_eq(&self, rhs: &Self) {
        if Rc::ptr_eq(&self.minv, &rhs.minv) {
            if self.minv.0 != rhs.minv.0 {
                panic!("The modulus of two operators should be the same!");
            }
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

#[cfg(std)]
impl<T: Integer + Montgomery> MontgomeryBigint<T>
where
    T::Double: From<T>,
{
    /// Convert n into the modulo ring ℤ/mℤ (i.e. `n % m`)
    #[inline]
    pub fn new(n: T, m: T) -> Self {
        let inv = Montgomery::neginv(&m);
        let a = Montgomery::transform(n, &m);
        MontgomeryBigint {
            a,
            minv: Rc::new((m, inv)),
        }
    }
}

impl<T: Integer + Montgomery> PartialEq for MontgomeryInt<T> {
    #[inline]
    fn eq(&self, other: &Self) -> bool {
        self.check_modulus_eq(other);
        self.a == other.a
    }
}

#[cfg(std)]
impl<T: Integer + Montgomery> PartialEq for MontgomeryBigint<T> {
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

#[cfg(std)]
impl<T: Integer + Montgomery> Add for MontgomeryBigint<T> {
    type Output = Self;

    #[inline]
    fn add(self, rhs: Self) -> Self::Output {
        self.check_modulus_eq(&rhs);
        let m = &Borrow::<(T, T::Inv)>::borrow(&self.minv).0;
        let a = Montgomery::add(&self.a, &rhs.a, m);
        MontgomeryBigint { a, minv: self.minv }
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

#[cfg(std)]
impl<T: Integer + Montgomery> Sub for MontgomeryBigint<T> {
    type Output = Self;

    #[inline]
    fn sub(self, rhs: Self) -> Self::Output {
        self.check_modulus_eq(&rhs);
        let m = &Borrow::<(T, T::Inv)>::borrow(&self.minv).0;
        let a = Montgomery::sub(&self.a, &rhs.a, m);
        MontgomeryBigint { a, minv: self.minv }
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

#[cfg(std)]
impl<T: Integer + Montgomery> Neg for MontgomeryBigint<T> {
    type Output = Self;

    #[inline]
    fn neg(self) -> Self::Output {
        let m = &Borrow::<(T, T::Inv)>::borrow(&self.minv).0;
        let a = Montgomery::neg(&self.a, m);
        MontgomeryBigint { a, minv: self.minv }
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

#[cfg(std)]
impl<T: Integer + Montgomery> Mul for MontgomeryBigint<T> {
    type Output = Self;

    #[inline]
    fn mul(self, rhs: Self) -> Self::Output {
        self.check_modulus_eq(&rhs);
        let minv = Borrow::<(T, T::Inv)>::borrow(&self.minv);
        let a = Montgomery::mul(&self.a, &rhs.a, &minv.0, &minv.1);
        MontgomeryBigint { a, minv: self.minv }
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

#[cfg(std)]
impl<T: Integer + Montgomery> Pow<T> for MontgomeryBigint<T> {
    type Output = Self;

    #[inline]
    fn pow(self, rhs: T) -> Self::Output {
        let minv = Borrow::<(T, T::Inv)>::borrow(&self.minv);
        let a = Montgomery::pow(&self.a, &rhs, &minv.0, &minv.1);
        MontgomeryBigint { a, minv: self.minv }
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
    fn new(&self, n: T) -> Self {
        let a = Montgomery::transform(n, &self.m);
        MontgomeryInt {
            a,
            m: self.m.clone(),
            mi: self.mi.clone(),
        }
    }
}

#[cfg(std)]
impl<T: Integer + Montgomery + Clone> ModularInteger for MontgomeryBigint<T>
where
    T::Double: From<T>,
{
    type Base = T;

    #[inline]
    fn modulus(&self) -> &T {
        &Borrow::<(T, T::Inv)>::borrow(&self.minv).0
    }

    #[inline]
    fn residue(&self) -> T {
        let minv = Borrow::<(T, T::Inv)>::borrow(&self.minv);
        Montgomery::reduce(T::Double::from(self.a.clone()), &minv.0, &minv.1)
    }

    #[inline]
    fn new(&self, n: T) -> Self {
        let m = &Borrow::<(T, T::Inv)>::borrow(&self.minv).0;
        let a = Montgomery::transform(n, &m);
        MontgomeryBigint {
            a,
            minv: self.minv.clone(),
        }
    }
}
