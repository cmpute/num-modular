use crate::{udouble, ModularInteger, ModularUnaryOps, Reducer};
use core::ops::*;
use num_traits::{Inv, Pow};

/// An integer in a modulo ring
#[derive(Debug, Clone, Copy)]
pub struct ReducedInt<T, R: Reducer<T>> {
    /// The reduced representation of the integer in a modulo ring.
    a: T,

    /// The reducer for the integer
    r: R,
}

impl<T, R: Reducer<T>> ReducedInt<T, R> {
    /// Convert n into the modulo ring ℤ/mℤ (i.e. `n % m`)
    #[inline]
    pub fn new(n: T, m: &T) -> Self {
        let r = R::new(m);
        let a = r.transform(n);
        Self { a, r }
    }

    #[inline(always)]
    fn check_modulus_eq(&self, rhs: &Self)
    where
        T: PartialEq,
    {
        // we don't directly compare m because m could be empty in case of Mersenne modular integer
        if cfg!(debug_assertions) && self.r.modulus() != rhs.r.modulus() {
            panic!("The modulus of two operators should be the same!");
        }
    }

    #[inline(always)]
    pub fn repr(&self) -> &T {
        &self.a
    }
}

impl<T: PartialEq, R: Reducer<T>> PartialEq for ReducedInt<T, R> {
    #[inline]
    fn eq(&self, other: &Self) -> bool {
        self.check_modulus_eq(other);
        self.a == other.a
    }
}

macro_rules! impl_binops {
    ($method:ident, impl $op:ident) => {
        impl<T: PartialEq, R: Reducer<T>> $op for ReducedInt<T, R> {
            type Output = Self;
            fn $method(self, rhs: Self) -> Self::Output {
                self.check_modulus_eq(&rhs);
                let Self { a, r } = self;
                let a = r.$method(a, rhs.a);
                Self { a, r }
            }
        }

        impl<T: PartialEq + Clone, R: Reducer<T>> $op<&Self> for ReducedInt<T, R> {
            type Output = Self;
            #[inline]
            fn $method(self, rhs: &Self) -> Self::Output {
                self.check_modulus_eq(&rhs);
                let Self { a, r } = self;
                let a = r.$method(a, rhs.a.clone());
                Self { a, r }
            }
        }

        impl<T: PartialEq + Clone, R: Reducer<T>> $op<ReducedInt<T, R>> for &ReducedInt<T, R> {
            type Output = ReducedInt<T, R>;
            #[inline]
            fn $method(self, rhs: ReducedInt<T, R>) -> Self::Output {
                self.check_modulus_eq(&rhs);
                let ReducedInt { a, r } = rhs;
                let a = r.$method(self.a.clone(), a);
                ReducedInt { a, r }
            }
        }

        impl<T: PartialEq + Clone, R: Reducer<T> + Clone> $op<&ReducedInt<T, R>>
            for &ReducedInt<T, R>
        {
            type Output = ReducedInt<T, R>;
            #[inline]
            fn $method(self, rhs: &ReducedInt<T, R>) -> Self::Output {
                self.check_modulus_eq(&rhs);
                let a = self.r.$method(self.a.clone(), rhs.a.clone());
                ReducedInt {
                    a,
                    r: self.r.clone(),
                }
            }
        }

        impl<T: PartialEq, R: Reducer<T>> $op<T> for ReducedInt<T, R> {
            type Output = Self;
            fn $method(self, rhs: T) -> Self::Output {
                let Self { a, r } = self;
                let rhs = r.transform(rhs);
                let a = r.$method(a, rhs);
                Self { a, r }
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
        let Self { a, r } = self;
        let a = r.neg(a);
        Self { a, r }
    }
}
impl<T: PartialEq + Clone, R: Reducer<T> + Clone> Neg for &ReducedInt<T, R> {
    type Output = ReducedInt<T, R>;
    #[inline]
    fn neg(self) -> Self::Output {
        let a = self.r.neg(self.a.clone());
        ReducedInt { a, r: self.r.clone() }
    }
}

impl<T: PartialEq, R: Reducer<T>> Inv for ReducedInt<T, R> {
    type Output = Self;
    #[inline]
    fn inv(self) -> Self::Output {
        let Self { a, r } = self;
        let a = r.inv(a).expect("the modular inverse doesn't exists.");
        Self { a, r }
    }
}
impl<T: PartialEq + Clone, R: Reducer<T> + Clone> Inv for &ReducedInt<T, R> {
    type Output = ReducedInt<T, R>;
    #[inline]
    fn inv(self) -> Self::Output {
        let a = self
            .r
            .inv(self.a.clone())
            .expect("the modular inverse doesn't exists.");
        ReducedInt { a, r: self.r.clone() }
    }
}

impl<T: PartialEq, R: Reducer<T>> Div for ReducedInt<T, R> {
    type Output = Self;
    #[inline]
    fn div(self, rhs: Self) -> Self::Output {
        self * rhs.inv()
    }
}
impl<T: PartialEq + Clone, R: Reducer<T> + Clone> Div<&ReducedInt<T, R>> for ReducedInt<T, R> {
    type Output = Self;
    #[inline]
    fn div(self, rhs: &Self) -> Self::Output {
        self * rhs.inv()
    }
}
impl<T: PartialEq + Clone, R: Reducer<T> + Clone> Div<ReducedInt<T, R>> for &ReducedInt<T, R> {
    type Output = ReducedInt<T, R>;
    #[inline]
    fn div(self, rhs: ReducedInt<T, R>) -> Self::Output {
        self.clone() * rhs.inv()
    }
}
impl<T: PartialEq + Clone, R: Reducer<T> + Clone> Div<&ReducedInt<T, R>> for &ReducedInt<T, R> {
    type Output = ReducedInt<T, R>;
    #[inline]
    fn div(self, rhs: &ReducedInt<T, R>) -> Self::Output {
        self.clone() * rhs.inv()
    }
}

impl<T: PartialEq, R: Reducer<T>> Pow<T> for ReducedInt<T, R> {
    type Output = Self;
    #[inline]
    fn pow(self, rhs: T) -> Self::Output {
        let Self { a, r } = self;
        let a = r.pow(a, rhs);
        Self { a, r }
    }
}
impl<T: PartialEq + Clone, R: Reducer<T> + Clone> Pow<T> for &ReducedInt<T, R> {
    type Output = ReducedInt<T, R>;
    #[inline]
    fn pow(self, rhs: T) -> Self::Output {
        let a = self.r.pow(self.a.clone(), rhs);
        ReducedInt { a, r: self.r.clone() }
    }
}

impl<T: PartialEq + Clone, R: Reducer<T> + Clone> ModularInteger for ReducedInt<T, R> {
    type Base = T;

    #[inline]
    fn modulus(&self) -> T {
        self.r.modulus()
    }

    #[inline(always)]
    fn residue(&self) -> T {
        self.r.residue(self.a.clone())
    }

    #[inline(always)]
    fn is_zero(&self) -> bool {
        self.r.is_zero(&self.a)
    }

    #[inline]
    fn convert(&self, n: T) -> Self {
        Self {
            a: self.r.transform(n),
            r: self.r.clone(),
        }
    }

    #[inline]
    fn double(self) -> Self {
        let Self { a, r } = self;
        let a = r.double(a);
        Self { a, r }
    }

    #[inline]
    fn square(self) -> Self {
        let Self { a, r } = self;
        let a = r.square(a);
        Self { a, r }
    }
}

// An vanilla reducer is also provided here
/// A plain reducer that just use normal [Rem] operators. It will keep the integer
/// in range [0, modulus) after each operation.
#[derive(Debug, Clone, Copy)]
pub struct Vanilla<T>(T);

macro_rules! impl_reduced_binary_pow {
    ($T:ty, $M:ty) => {
        fn pow(&self, base: $T, exp: $T) -> $T {
            match exp {
                1 => base,
                2 => self.square(base),
                e => {
                    let mut multi = base;
                    let mut exp = e;
                    let mut result = self.transform(1);
                    while exp > 0 {
                        if exp & 1 != 0 {
                            result = self.mul(result, multi);
                        }
                        multi = self.square(multi);
                        exp >>= 1;
                    }
                    result
                }
            }
        }
    };
}

pub(crate) use impl_reduced_binary_pow;

macro_rules! impl_uprim_vanilla_core {
    ($single:ty) => {
        #[inline(always)]
        fn new(m: &$single) -> Self {
            Self(*m)
        }
        #[inline(always)]
        fn transform(&self, target: $single) -> $single {
            target % self.0
        }
        #[inline(always)]
        fn residue(&self, target: $single) -> $single {
            target
        }
        #[inline(always)]
        fn modulus(&self) -> $single {
            self.0
        }
        #[inline(always)]
        fn is_zero(&self, target: &$single) -> bool {
            *target == 0
        }

        #[inline]
        fn add(&self, lhs: $single, rhs: $single) -> $single {
            let (sum, overflow) = lhs.overflowing_add(rhs);
            if overflow {
                sum + self.0.wrapping_neg()
            } else if sum >= self.0 {
                sum - self.0
            } else {
                sum
            }
        }

        #[inline]
        fn double(&self, target: $single) -> $single {
            self.add(target, target)
        }

        #[inline]
        fn sub(&self, lhs: $single, rhs: $single) -> $single {
            if lhs >= rhs {
                lhs - rhs
            } else {
                self.0 - (rhs - lhs)
            }
        }

        #[inline]
        fn neg(&self, target: $single) -> $single {
            if target == 0 {
                0
            } else {
                self.0 - target
            }
        }

        #[inline(always)]
        fn inv(&self, target: $single) -> Option<$single> {
            target.invm(&self.0)
        }

        impl_reduced_binary_pow!($single, $single);
    };
}

macro_rules! impl_uprim_vanilla {
    ($single:ty, $double:ty) => {
        impl Reducer<$single> for Vanilla<$single> {
            impl_uprim_vanilla_core!($single);

            #[inline]
            fn mul(&self, lhs: $single, rhs: $single) -> $single {
                ((lhs as $double) * (rhs as $double) % (self.0 as $double)) as $single
            }

            #[inline]
            fn square(&self, target: $single) -> $single {
                let target = target as $double;
                (target * target % (self.0 as $double)) as $single
            }
        }
    };
}

impl_uprim_vanilla!(u8, u16);
impl_uprim_vanilla!(u16, u32);
impl_uprim_vanilla!(u32, u64);
impl_uprim_vanilla!(u64, u128);
#[cfg(target_pointer_width = "32")]
impl_uprim_vanilla!(usize, u64);
#[cfg(target_pointer_width = "64")]
impl_uprim_vanilla!(usize, u128);

impl Reducer<u128> for Vanilla<u128> {
    impl_uprim_vanilla_core!(u128);

    #[inline]
    fn mul(&self, lhs: u128, rhs: u128) -> u128 {
        udouble::widening_mul(lhs, rhs) % self.0
    }

    #[inline]
    fn square(&self, target: u128) -> u128 {
        udouble::widening_square(target) % self.0
    }
}

/// An integer in modulo ring based on conventional [Rem] operations
pub type VanillaInt<T> = ReducedInt<T, Vanilla<T>>;

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{ModularCoreOps, ModularPow, ModularUnaryOps};
    use rand::random;

    const NRANDOM: u32 = 10;

    #[test]
    fn test_against_prim() {
        macro_rules! tests_for {
            ($($T:ty)*) => ($(
                let m = random::<$T>();
                let e = random::<u8>() as $T;
                let (a, b) = (random::<$T>(), random::<$T>());
                let am = VanillaInt::new(a, &m);
                let bm = VanillaInt::new(b, &m);
                assert_eq!((am + bm).residue(), a.addm(b, &m));
                assert_eq!((am - bm).residue(), a.subm(b, &m));
                assert_eq!((am * bm).residue(), a.mulm(b, &m));
                assert_eq!(am.neg().residue(), a.negm(&m));
                assert_eq!(am.double().residue(), a.dblm(&m));
                assert_eq!(am.square().residue(), a.sqm(&m));
                assert_eq!(am.pow(e).residue(), a.powm(e, &m));
                if let Some(v) = a.invm(&m) {
                    assert_eq!(am.inv().residue(), v);
                }
            )*);
        }

        for _ in 0..NRANDOM {
            tests_for!(u8 u16 u32 u64 u128 usize);
        }
    }
}
