use crate::{Reducer, ModularInteger};
use core::ops::{Add, Mul, Neg, Sub};
use num_traits::Pow;

/// This 
#[derive(Debug, Clone, Copy)]
pub struct ReducedInt<T, R: Reducer<T>> {
    /// The reduced representation of the integer in a modulo ring.
    a: T,

    /// The modulus.
    m: T,

    /// The reducer for the integer
    r: R,
}

impl<T, R: Reducer<T>> ReducedInt<T, R> {
    /// Convert n into the modulo ring ℤ/mℤ (i.e. `n % m`)
    #[inline]
    pub fn new(n: T, m: T) -> Self {
        let r = R::new(&m);
        let a = R::transform(n, &m);
        ReducedInt { a, m, r }
    }

    #[inline(always)]
    fn check_modulus_eq(&self, rhs: &Self) where T: Eq {
        if cfg!(debug_assertions) && self.m != rhs.m {
            panic!("The modulus of two operators should be the same!");
        }
    }

    #[inline(always)]
    pub fn is_zero(&self) -> bool {
        self.r.is_zero(&self.a, &self.m)
    }

    #[inline(always)]
    pub fn repr(&self) -> &T {
        &self.a
    }
}

impl<T: Eq, R: Reducer<T>> PartialEq for ReducedInt<T, R> {
    #[inline]
    fn eq(&self, other: &Self) -> bool {
        self.check_modulus_eq(other);
        self.a == other.a
    }
}

macro_rules! impl_binops {
    ($method:ident, impl $op:ident) => {        
        impl<T: Eq, R: Reducer<T>> $op for ReducedInt<T, R> {
            type Output = Self;
            fn $method(self, rhs: Self) -> Self::Output {
                self.check_modulus_eq(&rhs);
                let Self { a, m, r } = self;
                let a = r.$method(a, rhs.a, &m);
                Self { a, m, r }
            }
        }

        impl<T: Eq + Clone, R: Reducer<T>> $op<&Self> for ReducedInt<T, R> {
            type Output = Self;
            #[inline]
            fn $method(self, rhs: &Self) -> Self::Output {
                self.check_modulus_eq(&rhs);
                let Self { a, m, r } = self;
                let a = r.$method(a, rhs.a.clone(), &m);
                Self { a, m, r }
            }
        }

        impl<T: Eq + Clone, R: Reducer<T>> $op<ReducedInt<T, R>> for &ReducedInt<T, R> {
            type Output = ReducedInt<T, R>;
            #[inline]
            fn $method(self, rhs: ReducedInt<T, R>) -> Self::Output {
                self.check_modulus_eq(&rhs);
                let ReducedInt { a, m, r } = rhs;
                let a = r.$method(self.a.clone(), a, &m);
                ReducedInt { a, m, r }
            }
        }

        impl<T: Eq + Clone, R: Reducer<T> + Clone> $op<&ReducedInt<T, R>> for &ReducedInt<T, R> {
            type Output = ReducedInt<T, R>;
            #[inline]
            fn $method(self, rhs: &ReducedInt<T, R>) -> Self::Output {
                self.check_modulus_eq(&rhs);
                let a = self.r.$method(self.a.clone(), rhs.a.clone(), &self.m);
                ReducedInt { a, m: self.m.clone(), r: self.r.clone() }
            }
        }
    };
}
impl_binops!(add, impl Add);
impl_binops!(sub, impl Sub);
impl_binops!(mul, impl Mul);

impl<T: Eq, R: Reducer<T>> Neg for ReducedInt<T, R> {
    type Output = Self;
    #[inline]
    fn neg(self) -> Self::Output {
        let Self { a, m, r } = self;
        let a = r.neg(a, &m);
        Self { a, m, r }
    }
}

impl<T: Eq, R: Reducer<T>> Pow<T> for ReducedInt<T, R> {
    type Output = Self;
    #[inline]
    fn pow(self, rhs: T) -> Self::Output {
        let Self { a, m, r } = self;
        let a = r.pow(a, rhs, &m);
        Self { a, m, r }
    }
}

impl<T: Eq + Clone, R: Reducer<T> + Clone> ModularInteger for ReducedInt<T, R>
{
    type Base = T;

    #[inline]
    fn modulus(&self) -> &T {
        &self.m
    }

    #[inline]
    fn residue(&self) -> T {
        self.r.residue(self.a.clone(), &self.m)
    }

    #[inline]
    fn convert(&self, n: T) -> Self {
        Self {
            a: R::transform(n, &self.m),
            m: self.m.clone(),
            r: self.r.clone(),
        }
    }

    #[inline]
    fn double(self) -> Self {
        let Self { a, m, r } = self;
        let a = r.double(a, &m);
        Self { a, m, r }
    }

    #[inline]
    fn square(self) -> Self {
        let Self { a, m, r } = self;
        let a = r.square(a, &m);
        Self { a, m, r }
    }
}
