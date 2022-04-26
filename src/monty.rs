use crate::{udouble, ModularInteger};
use core::ops::{Add, Mul, Neg, Sub};
use num_integer::Integer;
use num_traits::Pow;

// TODO(v0.next): refactor the Montgomery trait into a "Reducer" trait?
//
// We just need to ensure that the behavior for multi-precision integers also follow this form, so this refactorization should
// be done after we have a implementation for big integers.
//
// This trait might look like 
trait Reducer {
    type Elem;

    fn transform(target: Self::Elem, m: &Self::Elem) -> Self::Elem;
    fn residue(&self, target: Self::Elem, m: &Self::Elem) -> Self::Elem;

    fn add(&self, lhs: &Self::Elem, rhs: &Self::Elem, m: &Self::Elem) -> Self::Elem;
    fn double(&self, target: &Self::Elem, m: &Self::Elem) -> Self::Elem;
    fn sub(&self, lhs: &Self::Elem, rhs: &Self::Elem, m: &Self::Elem) -> Self::Elem;
    fn neg(&self, target: &Self::Elem, m: &Self::Elem) -> Self::Elem;
    fn mul(&self, lhs: &Self::Elem, rhs: &Self::Elem, m: &Self::Elem) -> Self::Elem;
    fn square(&self, target: &Self::Elem, m: &Self::Elem) -> Self::Elem;
    fn pow(&self, base: &Self, exp: &Self, m: &Self::Elem) -> Self::Elem;
}
// Then
//
// struct Vanilla<T>(T): Reducer (trivial modular ring)
// struct Montgomery<T>(T, TInv): Reducer, implement for primitive integers
// struct MontgomeryMP<T>(Rc<(T, TInv)>): Reducer, implement for each bigint type, use relaxed form described in https://cetinkayakoc.net/docs/j56.pdf and https://eprint.iacr.org/2011/239.pdf
// struct Barret<T>(T, TInv): Reducer
// struct BarretMP<T>(Rc<(T, TInv)>): Reducer, implement for Vec[usize]
// struct ReducedInt<T, Reducer<Elem = T>>: ModularInteger
// type MontgomeryInt<T> = ReducedInt<T, Montgomery<T>>
// type BarretInt<T> = ReducedInt<T, Barret<T>>
//
// Besides, we could directly base the operations on a specific bigint library (ibig is a good candidate), and convert all other bigint types to this type for arithmetics (rug::Integer::as_limbs, num_bigint::BigUint::to_u32_digits/to_u64_digits)
//     but there're two problems for ibig-rs now: it's unable to construct and deconstruct big integer by moving, and it has not implemented the Integer trait yet.
// So the better option is to create a standalone "reduce" function that takes &[usize] as input, and then call this reduce function for each big integer backend.
// And consider create separate ReducedInt type for small and bigints, as we can provide convenient interface for multi-by-single operators

/// Operations of a integer represented in Montgomery form. Types implementing this
/// trait can be used to construct a [MontgomeryInt].
///
/// The generic type T represents the underlying integer representation, and
/// R=2^B will be used as the auxiliary modulus, where B is automatically selected
/// based on the size of T.
pub trait Montgomery {
    /// The type for inversion of the modulus.
    ///
    /// This type is usually the same as Self, but it can be smaller when using
    /// Montgomery form on multi-precision integer representations.
    type Inv;

    /// The type of integer with double width. It is only used in `reduce()`,
    /// so it's okay that it's not actually doubled with
    type Double;

    /// Calculate -(m^-1) mod R, return [None] if the inverse doesn't exist.
    fn neginv(m: &Self) -> Option<Self::Inv>;

    /// Transform a normal integer into Montgomery form (compute `target*R mod m`)
    fn transform(target: Self, m: &Self) -> Self;

    /// Transform a montgomery form back to normal integer (compute `monty/R mod m`)
    fn reduce(monty: Self::Double, m: &Self, minv: &Self::Inv) -> Self;

    /// Calculate (lhs + rhs) mod m in Montgomery form
    fn add(lhs: &Self, rhs: &Self, m: &Self) -> Self;

    /// Calculate (lhs - rhs) mod m in Montgomery form
    fn sub(lhs: &Self, rhs: &Self, m: &Self) -> Self;

    /// Calculate 2*monty mod m
    fn double(monty: &Self, m: &Self) -> Self where Self: Sized {
        Montgomery::add(monty, monty, m)
    }

    /// Calculate -monty mod m in Montgomery form
    fn neg(monty: &Self, m: &Self) -> Self;

    /// Calculate (lhs * rhs) mod m in Montgomery form
    fn mul(lhs: &Self, rhs: &Self, m: &Self, minv: &Self::Inv) -> Self;

    /// Calculate monty^2 mod m in Montgomery form
    fn square(monty: &Self, m: &Self, minv: &Self::Inv) -> Self where Self: Sized {
        Montgomery::mul(monty, monty, m, minv)
    }

    /// Calculate base ^ exp mod m in Montgomery form
    fn pow(base: &Self, exp: &Self, m: &Self, minv: &Self::Inv) -> Self;

    // TODO: support montgomery inverse, see http://cetinkayakoc.net/docs/j82.pdf
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

            if overflow {
                t + m.wrapping_neg()
            } else if &t >= m {
                t - m
            } else {
                t
            }
        }

        #[inline]
        fn mul(lhs: &Self, rhs: &Self, m: &Self, minv: &Self::Inv) -> Self {
            Montgomery::reduce((*lhs as Self::Double) * (*rhs as Self::Double), m, minv)
        }

        #[inline]
        fn square(monty: &Self, m: &Self, minv: &Self::Inv) -> Self {
            let d = *monty as Self::Double;
            Montgomery::reduce(d * d, m, minv)
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
            } else if &sum > m {
                sum - m
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

    fn neginv(m: &Self) -> Option<Self> {
        if m & 1 == 0 {
            return None;
        }
        let i = BINVERT_TABLE[((m >> 1) & 0x7F) as usize].wrapping_neg();
        Some(i)
    }
}

impl Montgomery for u16 {
    type Inv = u16;
    type Double = u32;

    impl_uprim_montgomery_core!();
    impl_uprim_montgomery!();

    fn neginv(m: &Self) -> Option<Self> {
        if m & 1 == 0 {
            return None;
        }
        let i = BINVERT_TABLE[((m >> 1) & 0x7F) as usize] as u16;
        // Newton-Rhapson iteration (hensel lifting), see https://arxiv.org/abs/1303.0328
        let i = i.wrapping_mul(*m).wrapping_sub(2).wrapping_mul(i);
        Some(i)
    }
}

impl Montgomery for u32 {
    type Inv = u32;
    type Double = u64;

    impl_uprim_montgomery_core!();
    impl_uprim_montgomery!();

    fn neginv(m: &Self) -> Option<Self> {
        if m & 1 == 0 {
            return None;
        }
        let i = BINVERT_TABLE[((m >> 1) & 0x7F) as usize] as u32;
        let i = 2u32.wrapping_sub(i.wrapping_mul(*m)).wrapping_mul(i);
        let i = i.wrapping_mul(*m).wrapping_sub(2).wrapping_mul(i);
        Some(i)
    }
}

impl Montgomery for u64 {
    type Inv = u64;
    type Double = u128;

    impl_uprim_montgomery_core!();
    impl_uprim_montgomery!();

    fn neginv(m: &Self) -> Option<Self> {
        if m & 1 == 0 {
            return None;
        }
        let i = BINVERT_TABLE[((m >> 1) & 0x7F) as usize] as u64;
        let i = 2u64.wrapping_sub(i.wrapping_mul(*m)).wrapping_mul(i);
        let i = 2u64.wrapping_sub(i.wrapping_mul(*m)).wrapping_mul(i);
        let i = i.wrapping_mul(*m).wrapping_sub(2).wrapping_mul(i);
        Some(i)
    }
}

impl Montgomery for u128 {
    type Inv = u128;
    type Double = udouble;

    fn neginv(m: &Self) -> Option<Self> {
        if m & 1 == 0 {
            return None;
        }
        let i = BINVERT_TABLE[((m >> 1) & 0x7F) as usize] as u128;
        let i = 2u128.wrapping_sub(i.wrapping_mul(*m)).wrapping_mul(i);
        let i = 2u128.wrapping_sub(i.wrapping_mul(*m)).wrapping_mul(i);
        let i = 2u128.wrapping_sub(i.wrapping_mul(*m)).wrapping_mul(i);
        let i = i.wrapping_mul(*m).wrapping_sub(2).wrapping_mul(i);
        Some(i)
    }

    #[inline]
    fn transform(target: Self, m: &Self) -> Self {
        if target == 0 {
            return 0;
        }
        udouble { hi: target, lo: 0 } % *m
    }

    // REDC algorithm
    fn reduce(monty: Self::Double, m: &Self, minv: &Self::Inv) -> Self {
        debug_assert!(monty < udouble { hi: *m, lo: 0 });

        let tm = monty.lo.wrapping_mul(*minv);
        let (t, overflow) = monty.overflowing_add(udouble::widening_mul(tm, *m));

        if overflow {
            t.hi + m.wrapping_neg()
        } else if &t.hi >= m {
            t.hi - m
        } else {
            t.hi
        }
    }

    #[inline]
    fn mul(lhs: &Self, rhs: &Self, m: &Self, minv: &Self::Inv) -> Self {
        Montgomery::reduce(udouble::widening_mul(*lhs, *rhs), m, minv)
    }

    #[inline]
    fn square(monty: &Self, m: &Self, minv: &Self::Inv) -> Self {
        Montgomery::reduce(udouble::widening_square(*monty), m, minv)
    }

    impl_uprim_montgomery!();
}

/// An integer represented in [Montgomery form](https://en.wikipedia.org/wiki/Montgomery_modular_multiplication#Montgomery_form).
/// 
/// It supports common operators just like a common integer, but it requires that the operands have the same modulus. Note that this condition
/// is only checked in the debug build.
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
    #[inline(always)]
    fn check_modulus_eq(&self, rhs: &Self) {
        if cfg!(debug_assertions) && self.m != rhs.m {
            panic!("The modulus of two operators should be the same!");
        }
    }
    #[inline(always)]
    pub fn is_zero(&self) -> bool {
        self.a.is_zero()
    }
    #[inline(always)]
    pub fn repr(&self) -> &T {
        &self.a
    }
}

// TODO(v0.4.x): accept even numbers by removing 2 factors from m and store the exponent
// Requirement: 1. A separate class to perform modular arithmetics with 2^n as modulus
//              2. Algorithm for construct residue from two components (see http://koclab.cs.ucsb.edu/teaching/cs154/docx/Notes7-Montgomery.pdf)
// Or we can just provide crt function, and let the implementation of monty int with full modulus support as an example code.

impl<T: Integer + Montgomery> MontgomeryInt<T>
where
    T::Double: From<T>,
{
    /// Convert n into the modulo ring ℤ/mℤ (i.e. `n % m`)
    #[inline]
    pub fn new(n: T, m: T) -> Self {
        let minv =
            Montgomery::neginv(&m).expect("the modulus has to be odd for 2^n based Montgomery");
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
        let Self { a, m, mi } = self;
        let a = Montgomery::add(&a, &rhs.a, &m);
        MontgomeryInt { a, m, mi }
    }
}

impl<T: Integer + Montgomery> Add<&Self> for MontgomeryInt<T> {
    type Output = Self;
    #[inline]
    fn add(self, rhs: &Self) -> Self::Output {
        self.check_modulus_eq(rhs);
        let Self { a, m, mi } = self;
        let a = Montgomery::add(&a, &rhs.a, &m);
        MontgomeryInt { a, m, mi }
    }
}
impl<T: Integer + Montgomery> Add<MontgomeryInt<T>> for &MontgomeryInt<T> {
    type Output = MontgomeryInt<T>;
    #[inline]
    fn add(self, rhs: MontgomeryInt<T>) -> Self::Output {
        self.check_modulus_eq(&rhs);
        let MontgomeryInt { a, m, mi } = rhs;
        let a = Montgomery::add(&self.a, &a, &m);
        MontgomeryInt { a, m, mi }
    }
}
impl<T: Integer + Montgomery + Clone> Add for &MontgomeryInt<T> 
where T::Inv: Clone {
    type Output = MontgomeryInt<T>;
    #[inline]
    fn add(self, rhs: Self) -> Self::Output {
        self.check_modulus_eq(rhs);
        let MontgomeryInt { a, m, mi } = self;
        let a = Montgomery::add(a, &rhs.a, m);
        MontgomeryInt { a, m: m.clone(), mi: mi.clone() }
    }
}

impl<T: Integer + Montgomery> Sub for MontgomeryInt<T> {
    type Output = Self;
    #[inline]
    fn sub(self, rhs: Self) -> Self::Output {
        self.check_modulus_eq(&rhs);
        let Self { a, m, mi } = self;
        let a = Montgomery::sub(&a, &rhs.a, &m);
        MontgomeryInt { a, m, mi }
    }
}
impl<T: Integer + Montgomery> Sub<&Self> for MontgomeryInt<T> {
    type Output = Self;
    #[inline]
    fn sub(self, rhs: &Self) -> Self::Output {
        self.check_modulus_eq(rhs);
        let Self { a, m, mi } = self;
        let a = Montgomery::sub(&a, &rhs.a, &m);
        MontgomeryInt { a, m, mi }
    }
}
impl<T: Integer + Montgomery> Sub<MontgomeryInt<T>> for &MontgomeryInt<T> {
    type Output = MontgomeryInt<T>;
    #[inline]
    fn sub(self, rhs: MontgomeryInt<T>) -> Self::Output {
        self.check_modulus_eq(&rhs);
        let MontgomeryInt { a, m, mi } = rhs;
        let a = Montgomery::sub(&self.a, &a, &m);
        MontgomeryInt { a, m, mi }
    }
}
impl<T: Integer + Montgomery + Clone> Sub for &MontgomeryInt<T> 
where T::Inv: Clone {
    type Output = MontgomeryInt<T>;
    #[inline]
    fn sub(self, rhs: Self) -> Self::Output {
        self.check_modulus_eq(rhs);
        let MontgomeryInt { a, m, mi } = self;
        let a = Montgomery::sub(a, &rhs.a, m);
        MontgomeryInt { a, m: m.clone(), mi: mi.clone() }
    }
}

impl<T: Integer + Montgomery> Neg for MontgomeryInt<T> {
    type Output = Self;
    #[inline]
    fn neg(self) -> Self::Output {
        let Self { a, m, mi } = self;
        let a = Montgomery::neg(&a, &m);
        MontgomeryInt { a, m, mi }
    }
}

impl<T: Integer + Montgomery> Mul for MontgomeryInt<T> {
    type Output = Self;
    #[inline]
    fn mul(self, rhs: Self) -> Self::Output {
        self.check_modulus_eq(&rhs);
        let Self { a, m, mi } = self;
        let a = Montgomery::mul(&a, &rhs.a, &m, &mi);
        MontgomeryInt { a, m, mi }
    }
}
impl<T: Integer + Montgomery> Mul<&Self> for MontgomeryInt<T> {
    type Output = Self;
    #[inline]
    fn mul(self, rhs: &Self) -> Self::Output {
        self.check_modulus_eq(rhs);
        let Self { a, m, mi } = self;
        let a = Montgomery::mul(&a, &rhs.a, &m, &mi);
        MontgomeryInt { a, m, mi }
    }
}
impl<T: Integer + Montgomery> Mul<MontgomeryInt<T>> for &MontgomeryInt<T> {
    type Output = MontgomeryInt<T>;
    #[inline]
    fn mul(self, rhs: MontgomeryInt<T>) -> Self::Output {
        self.check_modulus_eq(&rhs);
        let MontgomeryInt { a, m, mi } = rhs;
        let a = Montgomery::mul(&self.a, &a, &m, &mi);
        MontgomeryInt { a, m, mi }
    }
}
impl<T: Integer + Montgomery + Clone> Mul for &MontgomeryInt<T> 
where T::Inv: Clone {
    type Output = MontgomeryInt<T>;
    #[inline]
    fn mul(self, rhs: Self) -> Self::Output {
        self.check_modulus_eq(rhs);
        let MontgomeryInt { a, m, mi } = self;
        let a = Montgomery::mul(a, &rhs.a, m, mi);
        MontgomeryInt { a, m: m.clone(), mi: mi.clone() }
    }
}

impl<T: Integer + Montgomery> Pow<T> for MontgomeryInt<T> {
    type Output = Self;

    #[inline]
    fn pow(self, rhs: T) -> Self::Output {
        let Self { a, m, mi } = self;
        let a = Montgomery::pow(&a, &rhs, &m, &mi);
        MontgomeryInt { a, m, mi }
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
    fn convert(&self, n: T) -> Self {
        let a = Montgomery::transform(n, &self.m);
        MontgomeryInt {
            a,
            m: self.m.clone(),
            mi: self.mi.clone(),
        }
    }

    #[inline]
    fn double(self) -> Self {
        let Self { a, m, mi } = self;
        let a = Montgomery::double(&a, &m);
        MontgomeryInt { a, m, mi }
    }

    #[inline]
    fn square(self) -> Self {
        let Self { a, m, mi } = self;
        let a = Montgomery::square(&a, &m, &mi);
        MontgomeryInt { a, m, mi }
    }
}

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
        assert_eq!(MontgomeryInt::new(a, m).residue(), a % m);

        // random creation test
        for _ in 0..NRANDOM {
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
    fn test_against_prim() {
        for _ in 0..NRANDOM {
            let m = random::<u8>() | 1;
            let e = random::<u8>();
            let (a, b) = (random::<u8>(), random::<u8>());
            let am = MontgomeryInt::new(a, m);
            let bm = am.convert(b);
            assert_eq!((am + bm).residue(), a.addm(b, &m));
            assert_eq!((am - bm).residue(), a.subm(b, &m));
            assert_eq!((am * bm).residue(), a.mulm(b, &m));
            assert_eq!((-am).residue(), a.negm(&m));
            assert_eq!(am.pow(e).residue(), a.powm(e, &m));
            assert_eq!(am.double().residue(), a.dblm(&m));
            assert_eq!(am.square().residue(), a.sqm(&m));

            let m = random::<u16>() | 1;
            let e = e as u16;
            let (a, b) = (random::<u16>(), random::<u16>());
            let am = MontgomeryInt::new(a, m);
            let bm = am.convert(b);
            assert_eq!((am + bm).residue(), a.addm(b, &m));
            assert_eq!((am - bm).residue(), a.subm(b, &m));
            assert_eq!((am * bm).residue(), a.mulm(b, &m));
            assert_eq!((-am).residue(), a.negm(&m));
            assert_eq!(am.pow(e).residue(), a.powm(e, &m));
            assert_eq!(am.double().residue(), a.dblm(&m));
            assert_eq!(am.square().residue(), a.sqm(&m));

            let m = random::<u32>() | 1;
            let e = e as u32;
            let (a, b) = (random::<u32>(), random::<u32>());
            let am = MontgomeryInt::new(a, m);
            let bm = am.convert(b);
            assert_eq!((am + bm).residue(), a.addm(b, &m));
            assert_eq!((am - bm).residue(), a.subm(b, &m));
            assert_eq!((am * bm).residue(), a.mulm(b, &m));
            assert_eq!((-am).residue(), a.negm(&m));
            assert_eq!(am.pow(e).residue(), a.powm(e, &m));
            assert_eq!(am.double().residue(), a.dblm(&m));
            assert_eq!(am.square().residue(), a.sqm(&m));

            let m = random::<u64>() | 1;
            let e = e as u64;
            let (a, b) = (random::<u64>(), random::<u64>());
            let am = MontgomeryInt::new(a, m);
            let bm = am.convert(b);
            assert_eq!((am + bm).residue(), a.addm(b, &m));
            assert_eq!((am - bm).residue(), a.subm(b, &m));
            assert_eq!((am * bm).residue(), a.mulm(b, &m));
            assert_eq!((-am).residue(), a.negm(&m));
            assert_eq!(am.pow(e).residue(), a.powm(e, &m));
            assert_eq!(am.double().residue(), a.dblm(&m));
            assert_eq!(am.square().residue(), a.sqm(&m));

            let m = random::<u128>() | 1;
            let e = e as u128;
            let (a, b) = (random::<u128>(), random::<u128>());
            let am = MontgomeryInt::new(a, m);
            let bm = am.convert(b);
            assert_eq!((am + bm).residue(), a.addm(b, &m));
            assert_eq!((am - bm).residue(), a.subm(b, &m));
            assert_eq!((am * bm).residue(), a.mulm(b, &m));
            assert_eq!((-am).residue(), a.negm(&m));
            assert_eq!(am.pow(e).residue(), a.powm(e, &m));
            assert_eq!(am.double().residue(), a.dblm(&m));
            assert_eq!(am.square().residue(), a.sqm(&m));
        }
    }
}
