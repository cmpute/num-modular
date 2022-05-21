use crate::{Reducer, ModularInteger};
use core::ops::*;
use num_traits::Pow;

/// An integer in a modulo ring
#[derive(Debug, Clone, Copy)]
pub struct ReducedInt<T, R: Reducer<T>> {
    /// The reduced representation of the integer in a modulo ring.
    a: T,

    /// The modulus.
    m: R::Modulus,

    /// The reducer for the integer
    r: R,
}

impl<T, R: Reducer<T>> ReducedInt<T, R> {
    /// Convert n into the modulo ring ℤ/mℤ (i.e. `n % m`)
    #[inline]
    pub fn new(n: T, m: R::Modulus) -> Self {
        let r = R::new(&m);
        let a = R::transform(n, &m);
        ReducedInt { a, m, r }
    }

    #[inline(always)]
    fn check_modulus_eq(&self, rhs: &Self) where R::Modulus: PartialEq {
        if cfg!(debug_assertions) && self.m != rhs.m {
            panic!("The modulus of two operators should be the same!");
        }
    }

    #[inline(always)]
    pub fn repr(&self) -> &T {
        &self.a
    }
}

impl<T: PartialEq, R: Reducer<T>> PartialEq for ReducedInt<T, R> where R::Modulus: PartialEq {
    #[inline]
    fn eq(&self, other: &Self) -> bool {
        self.check_modulus_eq(other);
        self.a == other.a
    }
}

macro_rules! impl_binops {
    ($method:ident, impl $op:ident) => {        
        impl<T: PartialEq, R: Reducer<T>> $op for ReducedInt<T, R> where R::Modulus: PartialEq {
            type Output = Self;
            fn $method(self, rhs: Self) -> Self::Output {
                self.check_modulus_eq(&rhs);
                let Self { a, m, r } = self;
                let a = r.$method(a, rhs.a, &m);
                Self { a, m, r }
            }
        }

        impl<T: PartialEq + Clone, R: Reducer<T>> $op<&Self> for ReducedInt<T, R> where R::Modulus: PartialEq {
            type Output = Self;
            #[inline]
            fn $method(self, rhs: &Self) -> Self::Output {
                self.check_modulus_eq(&rhs);
                let Self { a, m, r } = self;
                let a = r.$method(a, rhs.a.clone(), &m);
                Self { a, m, r }
            }
        }

        impl<T: PartialEq + Clone, R: Reducer<T>> $op<ReducedInt<T, R>> for &ReducedInt<T, R> where R::Modulus: PartialEq {
            type Output = ReducedInt<T, R>;
            #[inline]
            fn $method(self, rhs: ReducedInt<T, R>) -> Self::Output {
                self.check_modulus_eq(&rhs);
                let ReducedInt { a, m, r } = rhs;
                let a = r.$method(self.a.clone(), a, &m);
                ReducedInt { a, m, r }
            }
        }

        impl<T: PartialEq + Clone, R: Reducer<T> + Clone> $op<&ReducedInt<T, R>> for &ReducedInt<T, R> where R::Modulus: PartialEq + Clone {
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

impl<T: PartialEq, R: Reducer<T>> Neg for ReducedInt<T, R> {
    type Output = Self;
    #[inline]
    fn neg(self) -> Self::Output {
        let Self { a, m, r } = self;
        let a = r.neg(a, &m);
        Self { a, m, r }
    }
}

impl<T: PartialEq, R: Reducer<T>> Pow<T> for ReducedInt<T, R> {
    type Output = Self;
    #[inline]
    fn pow(self, rhs: T) -> Self::Output {
        let Self { a, m, r } = self;
        let a = r.pow(a, rhs, &m);
        Self { a, m, r }
    }
}

impl<T: PartialEq + Clone, R: Reducer<T> + Clone> ModularInteger for ReducedInt<T, R> where R::Modulus: PartialEq + Clone
{
    type Base = T;

    #[inline]
    fn modulus(&self) -> T {
        R::modulus(&self.m)
    }

    #[inline(always)]
    fn residue(&self) -> T {
        self.r.residue(self.a.clone(), &self.m)
    }
    
    #[inline(always)]
    fn is_zero(&self) -> bool {
        self.r.is_zero(&self.a, &self.m)
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

// An vanilla reducer is also provided here
/// A plain reducer that just use normal [Rem] operators. It will keep the integer
/// in range [0, modulus) after each operation.
#[derive(Debug, Clone, Copy)]
pub struct Vanilla();

macro_rules! impl_uprim_vanilla_core {
    ($single:ty) => {
        #[inline(always)]
        fn new(_: &Self::Modulus) -> Self {
            Self {}
        }
        #[inline(always)]
        fn transform(target: $single, m: &Self::Modulus) -> $single {
            target % m
        }
        #[inline(always)]
        fn residue(&self, target: $single, _: &Self::Modulus) -> $single {
            target
        }
        #[inline(always)]
        fn modulus(m: &Self::Modulus) -> $single {
            *m
        }
        #[inline(always)]
        fn is_zero(&self, target: &$single, _: &$single) -> bool {
            *target == 0
        }

        #[inline]
        fn add(&self, lhs: $single, rhs: $single, m: &$single) -> $single {
            let (sum, overflow) = lhs.overflowing_add(rhs);
            if overflow {
                sum + m.wrapping_neg()
            } else if &sum >= m {
                sum - m
            } else {
                sum
            }
        }

        #[inline]
        fn double(&self, target: $single, m: &$single) -> $single {
            self.add(target, target, m)
        }

        #[inline]
        fn sub(&self, lhs: $single, rhs: $single, m: &$single) -> $single {
            if lhs >= rhs {
                lhs - rhs
            } else {
                m - (rhs - lhs)
            }
        }

        #[inline]
        fn neg(&self, target: $single, m: &$single) -> $single {
            if target == 0 {
                0
            } else {
                m - target
            }
        }
        
        fn pow(&self, base: $single, exp: $single, m: &$single) -> $single {
            match exp {
                1 => base,
                2 => self.square(base, m),
                e => {
                    let mut multi = base;
                    let mut exp = e;
                    let mut result = 1;
                    while exp > 0 {
                        if exp & 1 != 0 {
                            result = self.mul(result, multi, m);
                        }
                        multi = self.square(multi, m);
                        exp >>= 1;
                    }
                    result
                }
            }
        }
    }
}

macro_rules! impl_uprim_vanilla {
    ($single:ty, $double:ty) => {
        impl Reducer<$single> for Vanilla {
            type Modulus = $single;
            impl_uprim_vanilla_core!($single);

            #[inline]
            fn mul(&self, lhs: $single, rhs: $single, m: &$single) -> $single {
                ((lhs as $double) * (rhs as $double) % (*m as $double)) as $single
            }

            #[inline]
            fn square(&self, target: $single, m: &$single) -> $single {
                let target = target as $double;
                (target * target % (*m as $double)) as $single
            }
        }
    };
}

impl_uprim_vanilla!(u8, u16);

// TODO(v0.5): add test for vanilla integer
