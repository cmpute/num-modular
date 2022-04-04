//! Implementations for modular operations on primitive integers

use crate::{ModularCoreOps, ModularOps, ModularAbs};
use num_integer::Integer;

// TODO: implement the modular functions as const: https://github.com/rust-lang/rust/pull/68847

macro_rules! impl_powm_uprim {
    ($T:ty) => {
        fn powm(self, exp: $T, m: &$T) -> $T {
            match exp {
                1 => self % m,
                2 => self.mulm(self, m),
                _ => {
                    let mut multi = self % m;
                    let mut exp = exp;
                    let mut result = 1;
                    while exp > 0 {
                        if exp & 1 != 0 {
                            result = result.mulm(multi, m);
                        }
                        multi = multi.mulm(multi, m);
                        exp >>= 1;
                    }
                    result
                }
            }
        }
    };
}

macro_rules! impl_jacobi_uprim {
    ($T:ty) => {
        #[inline]
        fn legendre(self, n: &$T) -> i8 {
            match self.powm((n - 1) >> 1, &n) {
                0 => 0,
                1 => 1,
                x if x == n - 1 => -1,
                _ => panic!("n is not prime!"),
            }
        }

        fn jacobi(self, n: &$T) -> i8 {
            if n % 2 == 0 || n < &0 {
                panic!("The Jacobi symbol is only defined for non-negative odd integers!")
            }

            if self == 0 {
                return 0;
            }
            if self == 1 {
                return 1;
            }

            let mut a = self % n;
            let mut n = n.clone();
            let mut t = 1;
            while a > 0 {
                while (a & 1) == 0 {
                    a = a / 2;
                    if n & 7 == 3 || n & 7 == 5 {
                        t *= -1;
                    }
                }
                core::mem::swap(&mut a, &mut n);
                if (a & 3) == 3 && (n & 3) == 3 {
                    t *= -1;
                }
                a = a % n;
            }
            if n == 1 {
                t
            } else {
                0
            }
        }

        #[inline]
        fn kronecker(self, n: &$T) -> i8 {
            match n {
                0 => {
                    if self == 1 {
                        1
                    } else {
                        0
                    }
                }
                1 => 1,
                2 => {
                    if self & 1 == 0 {
                        0
                    } else if self & 7 == 1 || self & 7 == 7 {
                        1
                    } else {
                        -1
                    }
                }
                _ => {
                    let f = n.trailing_zeros();
                    let n = n >> f;
                    ModularOps::<&$T, &$T>::kronecker(self, &2).pow(f)
                        * ModularOps::<&$T, &$T>::jacobi(self, &n)
                }
            }
        }
    };
}

// implement inverse mod using extended euclidean algorithm
macro_rules! impl_invm_uprim {
    ($T:ty) => {
        fn invm(self, m: &$T) -> Option<Self::Output> {
            // TODO: optimize using https://eprint.iacr.org/2020/972.pdf
            let x = if &self >= m { self % m } else { self.clone() };

            let (mut last_r, mut r) = (m.clone(), x);
            let (mut last_t, mut t) = (0, 1);

            while r > 0 {
                let (quo, rem) = last_r.div_rem(&r);
                last_r = r;
                r = rem;

                let new_t = last_t.subm(quo.mulm(t, m), m);
                last_t = t;
                t = new_t;
            }

            // if r = gcd(self, m) > 1, then inverse doesn't exist
            if last_r > 1 {
                None
            } else {
                Some(last_t)
            }
        }
    };
}

macro_rules! impl_mod_arithm_uu {
    ($T:ty, $Tdouble:ty) => {
        impl ModularCoreOps<$T, &$T> for $T {
            type Output = $T;
            #[inline]
            fn addm(self, rhs: $T, m: &$T) -> $T {
                (((self as $Tdouble) + (rhs as $Tdouble)) % (*m as $Tdouble)) as $T
            }
            #[inline]
            fn subm(self, rhs: $T, m: &$T) -> $T {
                let (lhs, rhs) = (self % m, rhs % m);
                if lhs >= rhs {
                    lhs - rhs
                } else {
                    m - (rhs - lhs)
                }
            }
            #[inline]
            fn mulm(self, rhs: $T, m: &$T) -> $T {
                (((self as $Tdouble) * (rhs as $Tdouble)) % (*m as $Tdouble)) as $T
            }
            #[inline]
            fn negm(self, m: &$T) -> $T {
                let x = self % m;
                if x == 0 {
                    0
                } else {
                    m - x
                }
            }
        }

        impl ModularOps<$T, &$T> for $T {
            impl_powm_uprim!($T);
            impl_jacobi_uprim!($T);
            impl_invm_uprim!($T);
        }
    };
}

impl_mod_arithm_uu!(u8, u16);
impl_mod_arithm_uu!(u16, u32);
impl_mod_arithm_uu!(u32, u64);
impl_mod_arithm_uu!(u64, u128);
impl_mod_arithm_uu!(usize, u128);

impl ModularCoreOps<u128, &u128> for u128 {
    type Output = u128;

    // XXX: check if these operations are also faster in u64
    #[inline]
    fn addm(self, rhs: u128, m: &u128) -> u128 {
        if let Some(ab) = self.checked_add(rhs) {
            return ab % m;
        }

        let (lhs, rhs) = (self % m, rhs % m);
        if lhs < m - rhs {
            lhs + rhs
        } else {
            lhs.min(rhs) - (m - lhs.max(rhs))
        }
    }

    #[inline]
    fn subm(self, rhs: u128, m: &u128) -> u128 {
        let (lhs, rhs) = (self % m, rhs % m);
        if lhs >= rhs {
            lhs - rhs
        } else {
            m - (rhs - lhs)
        }
    }

    // XXX: benchmark against http://www.janfeitsma.nl/math/psp2/expmod
    // XXX: benchmark against udouble implementation
    fn mulm(self, rhs: u128, m: &u128) -> u128 {
        if let Some(ab) = self.checked_mul(rhs) {
            return ab % m;
        }

        let mut a = self % m;
        let mut b = rhs % m;

        if let Some(ab) = a.checked_mul(b) {
            return ab % m;
        }

        let mut result: u128 = 0;
        while b > 0 {
            if b & 1 > 0 {
                result = result.addm(a, m);
            }
            a = a.addm(a, m);
            b >>= 1;
        }
        result
    }

    #[inline]
    fn negm(self, m: &u128) -> u128 {
        let x = self % m;
        if x == 0 {
            0
        } else {
            m - x
        }
    }
}
impl ModularOps<u128, &u128> for u128 {
    impl_powm_uprim!(u128);
    impl_jacobi_uprim!(u128);
    impl_invm_uprim!(u128);
}

macro_rules! impl_mod_arithm_by_deref {
    ($($T:ty)*) => {$(
        impl ModularCoreOps<$T, &$T> for &$T {
            type Output = $T;
            #[inline]
            fn addm(self, rhs: $T, m: &$T) -> $T {
                (*self).addm(rhs, &m)
            }
            #[inline]
            fn subm(self, rhs: $T, m: &$T) -> $T {
                (*self).subm(rhs, &m)
            }
            #[inline]
            fn mulm(self, rhs: $T, m: &$T) -> $T {
                (*self).mulm(rhs, &m)
            }
            #[inline]
            fn negm(self, m: &$T) -> $T {
                ModularCoreOps::<$T, &$T>::negm(*self, m)
            }
        }
        impl ModularOps<$T, &$T> for &$T {
            #[inline]
            fn powm(self, exp: $T, m: &$T) -> $T {
                (*self).powm(exp, &m)
            }
            #[inline]
            fn invm(self, m: &$T) -> Option<$T> {
                ModularOps::<$T, &$T>::invm(*self, m)
            }
            #[inline]
            fn legendre(self, n: &$T) -> i8 {
                ModularOps::<$T, &$T>::legendre(*self, n)
            }
            #[inline]
            fn jacobi(self, n: &$T) -> i8 {
                ModularOps::<$T, &$T>::jacobi(*self, n)
            }
            #[inline]
            fn kronecker(self, n: &$T) -> i8 {
                ModularOps::<$T, &$T>::kronecker(*self, n)
            }
        }

        impl ModularCoreOps<&$T, &$T> for $T {
            type Output = $T;
            #[inline]
            fn addm(self, rhs: &$T, m: &$T) -> $T {
                self.addm(*rhs, &m)
            }
            #[inline]
            fn subm(self, rhs: &$T, m: &$T) -> $T {
                self.subm(*rhs, &m)
            }
            #[inline]
            fn mulm(self, rhs: &$T, m: &$T) -> $T {
                self.mulm(*rhs, &m)
            }
            #[inline]
            fn negm(self, m: &$T) -> $T {
                ModularCoreOps::<$T, &$T>::negm(self, m)
            }
        }
        impl ModularOps<&$T, &$T> for $T {
            #[inline]
            fn powm(self, exp: &$T, m: &$T) -> $T {
                self.powm(*exp, &m)
            }
            #[inline]
            fn invm(self, m: &$T) -> Option<$T> {
                ModularOps::<$T, &$T>::invm(self, m)
            }
            #[inline]
            fn legendre(self, n: &$T) -> i8 {
                ModularOps::<$T, &$T>::legendre(self, n)
            }
            #[inline]
            fn jacobi(self, n: &$T) -> i8 {
                ModularOps::<$T, &$T>::jacobi(self, n)
            }
            #[inline]
            fn kronecker(self, n: &$T) -> i8 {
                ModularOps::<$T, &$T>::kronecker(self, n)
            }
        }

        impl ModularCoreOps<&$T, &$T> for &$T {
            type Output = $T;
            #[inline]
            fn addm(self, rhs: &$T, m: &$T) -> $T {
                (*self).addm(*rhs, &m)
            }
            #[inline]
            fn subm(self, rhs: &$T, m: &$T) -> $T {
                (*self).subm(*rhs, &m)
            }
            #[inline]
            fn mulm(self, rhs: &$T, m: &$T) -> $T {
                (*self).mulm(*rhs, &m)
            }
            #[inline]
            fn negm(self, m: &$T) -> $T {
                ModularCoreOps::<$T, &$T>::negm(*self, m)
            }
        }
        impl ModularOps<&$T, &$T> for &$T {
            #[inline]
            fn powm(self, exp: &$T, m: &$T) -> $T {
                (*self).powm(*exp, &m)
            }
            #[inline]
            fn invm(self, m: &$T) -> Option<$T> {
                ModularOps::<$T, &$T>::invm(*self, m)
            }
            #[inline]
            fn legendre(self, n: &$T) -> i8 {
                ModularOps::<$T, &$T>::legendre(*self, n)
            }
            #[inline]
            fn jacobi(self, n: &$T) -> i8 {
                ModularOps::<$T, &$T>::jacobi(*self, n)
            }
            #[inline]
            fn kronecker(self, n: &$T) -> i8 {
                ModularOps::<$T, &$T>::kronecker(*self, n)
            }
        }
    )*};
}

impl_mod_arithm_by_deref!(u8 u16 u32 u64 u128 usize);

macro_rules! impl_absm_for_prim {
    ($($signed:ty => $unsigned:ty;)*) => {$(
        impl ModularAbs<$unsigned> for $signed {
            fn absm(self, m: &$unsigned) -> $unsigned {
                if self >= 0 {
                    (self as $unsigned) % m
                } else {
                    ModularCoreOps::<$unsigned>::negm(&(-self as $unsigned), m)
                }
            }
        }
    )*};
}

impl_absm_for_prim! {
    i8 => u8; i16 => u16; i32 => u32; i64 => u64; i128 => u128; isize => usize;
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::random;
    use core::ops::Neg;

    const ADDM_CASES: [(u8, u8, u8, u8); 10] = [
        // [m, x, y, rem]: x + y = rem (mod m)
        (5, 0, 0, 0),
        (5, 1, 2, 3),
        (5, 2, 1, 3),
        (5, 2, 2, 4),
        (5, 3, 2, 0),
        (5, 2, 3, 0),
        (5, 6, 1, 2),
        (5, 1, 6, 2),
        (5, 11, 7, 3),
        (5, 7, 11, 3),
    ];

    #[test]
    fn addm_test() {
        // fixed cases
        for &(m, x, y, r) in ADDM_CASES.iter() {
            assert_eq!(x.addm(y, &m), r);
            assert_eq!((x as u16).addm(y as u16, &(m as _)), r as _);
            assert_eq!((x as u32).addm(y as u32, &(m as _)), r as _);
            assert_eq!((x as u64).addm(y as u64, &(m as _)), r as _);
            assert_eq!((x as u128).addm(y as u128, &(m as _)), r as _);
        }

        // random cases for u64 and u128
        for _ in 0..10 {
            let a = random::<u32>() as u64;
            let b = random::<u32>() as u64;
            let m = random::<u32>() as u64;
            assert_eq!(a.addm(b, &m), (a + b) % m);
            assert_eq!(a.addm(b, &(1u64 << 32)) as u32, (a as u32).wrapping_add(b as u32));
            
            let a = random::<u64>() as u128;
            let b = random::<u64>() as u128;
            let m = random::<u64>() as u128;
            assert_eq!(a.addm(b, &m), (a + b) % m);
            assert_eq!(a.addm(b, &(1u128 << 64)) as u64, (a as u64).wrapping_add(b as u64));
        }
    }

    const SUBM_CASES: [(u8, u8, u8, u8); 10] = [
        // [m, x, y, rem]: x - y = rem (mod m)
        (7, 0, 0, 0),
        (7, 11, 9, 2),
        (7, 5, 2, 3),
        (7, 2, 5, 4),
        (7, 6, 7, 6),
        (7, 1, 7, 1),
        (7, 7, 1, 6),
        (7, 0, 6, 1),
        (7, 15, 1, 0),
        (7, 1, 15, 0),
    ];

    #[test]
    fn subm_test() {
        // fixed cases
        for &(m, x, y, r) in SUBM_CASES.iter() {
            assert_eq!(x.subm(y, &m), r);
            assert_eq!((x as u16).subm(y as u16, &(m as _)), r as _);
            assert_eq!((x as u32).subm(y as u32, &(m as _)), r as _);
            assert_eq!((x as u64).subm(y as u64, &(m as _)), r as _);
            assert_eq!((x as u128).subm(y as u128, &(m as _)), r as _);
        }
        
        // random cases for u64 and u128
        for _ in 0..10 {
            let a = random::<u32>() as u64;
            let b = random::<u32>() as u64;
            let m = random::<u32>() as u64;
            assert_eq!(a.subm(b, &m), (a as i64 - b as i64).rem_euclid(m as i64) as u64);
            assert_eq!(a.subm(b, &(1u64 << 32)) as u32, (a as u32).wrapping_sub(b as u32));
            
            let a = random::<u64>() as u128;
            let b = random::<u64>() as u128;
            let m = random::<u64>() as u128;
            assert_eq!(a.subm(b, &m), (a as i128 - b as i128).rem_euclid(m as i128) as u128);
            assert_eq!(a.subm(b, &(1u128 << 64)) as u64, (a as u64).wrapping_sub(b as u64));
        }
    }
    
    const NEGM_CASES: [(u8, u8, u8); 5] = [
        // [m, x, rem]: -x = rem (mod m)
        (5, 0, 0),
        (5, 2, 3),
        (5, 1, 4),
        (5, 5, 0),
        (5, 12, 3),
    ];

    #[test]
    fn negm_and_absm_test() {
        // fixed cases
        for &(m, x, r) in NEGM_CASES.iter() {
            assert_eq!(ModularCoreOps::<&u8>::negm(&x, &m), r);
            assert_eq!((x as i8).neg().absm(&m), r);
            assert_eq!(ModularCoreOps::<&u16>::negm(&(x as _), &(m as _)), r as _);
            assert_eq!((x as i16).neg().absm(&(m as u16)), r as _);
            assert_eq!(ModularCoreOps::<&u32>::negm(&(x as _), &(m as _)), r as _);
            assert_eq!((x as i32).neg().absm(&(m as u32)), r as _);
            assert_eq!(ModularCoreOps::<&u64>::negm(&(x as _), &(m as _)), r as _);
            assert_eq!((x as i64).neg().absm(&(m as u64)), r as _);
            assert_eq!(ModularCoreOps::<&u128>::negm(&(x as _), &(m as _)), r as _);
            assert_eq!((x as i128).neg().absm(&(m as u128)), r as _);
        }

        // random cases for u64 and u128
        for _ in 0..10 {
            let a = random::<u32>() as u64;
            let m = random::<u32>() as u64;
            assert_eq!(ModularCoreOps::<&u64>::negm(&a, &m), (a as i64).neg().rem_euclid(m as i64) as u64);
            assert_eq!(ModularCoreOps::<&u64>::negm(&a, &(1u64 << 32)) as u32, (a as u32).wrapping_neg());
            
            let a = random::<u64>() as u128;
            let m = random::<u64>() as u128;
            assert_eq!(ModularCoreOps::<&u128>::negm(&a, &m), (a as i128).neg().rem_euclid(m as i128) as u128);
            assert_eq!(ModularCoreOps::<&u128>::negm(&a, &(1u128 << 64)) as u64, (a as u64).wrapping_neg());
        }
    }

    const MULM_CASES: [(u8, u8, u8, u8); 10] = [
        // [m, x, y, rem]: x*y = rem (mod m)
        (7, 0, 0, 0),
        (7, 11, 9, 1),
        (7, 5, 2, 3),
        (7, 2, 5, 3),
        (7, 6, 7, 0),
        (7, 1, 7, 0),
        (7, 7, 1, 0),
        (7, 0, 6, 0),
        (7, 15, 1, 1),
        (7, 1, 15, 1),
    ];

    #[test]
    fn mulm_test() {
        // fixed cases
        for &(m, x, y, r) in MULM_CASES.iter() {
            assert_eq!(x.mulm(y, &m), r);
            assert_eq!((x as u16).mulm(y as u16, &(m as _)), r as _);
            assert_eq!((x as u32).mulm(y as u32, &(m as _)), r as _);
            assert_eq!((x as u64).mulm(y as u64, &(m as _)), r as _);
            assert_eq!((x as u128).mulm(y as u128, &(m as _)), r as _);
        }

        // random cases for u64 and u128
        for _ in 0..10 {
            let a = random::<u32>() as u64;
            let b = random::<u32>() as u64;
            let m = random::<u32>() as u64;
            assert_eq!(a.mulm(b, &m), (a * b) % m);
            assert_eq!(a.mulm(b, &(1u64 << 32)) as u32, (a as u32).wrapping_mul(b as u32));
            
            let a = random::<u64>() as u128;
            let b = random::<u64>() as u128;
            let m = random::<u64>() as u128;
            assert_eq!(a.mulm(b, &m), (a * b) % m);
            assert_eq!(a.mulm(b, &(1u128 << 32)) as u32, (a as u32).wrapping_mul(b as u32));
        }
    }

    const POWM_CASES: [(u8, u8, u8, u8); 10] = [
        // [m, x, y, rem]: x^y = rem (mod m)
        (7, 0, 0, 1),
        (7, 11, 9, 1),
        (7, 5, 2, 4),
        (7, 2, 5, 4),
        (7, 6, 7, 6),
        (7, 1, 7, 1),
        (7, 7, 1, 0),
        (7, 0, 6, 0),
        (7, 15, 1, 1),
        (7, 1, 15, 1),
    ];

    #[test]
    fn powm_test() {
        // fixed cases
        for &(m, x, y, r) in POWM_CASES.iter() {
            assert_eq!(x.powm(y, &m), r);
            assert_eq!((x as u16).powm(y as u16, &(m as _)), r as _);
            assert_eq!((x as u32).powm(y as u32, &(m as _)), r as _);
            assert_eq!((x as u64).powm(y as u64, &(m as _)), r as _);
            assert_eq!((x as u128).powm(y as u128, &(m as _)), r as _);
        }
    }

    const INVM_CASES: [(u64, u64, u64); 8] = [
        // [a, m, x] s.t. a*x = 1 (mod m) is satisfied
        (5, 11, 9),
        (8, 11, 7),
        (10, 11, 10),
        (3, 5000, 1667),
        (1667, 5000, 3),
        (999, 5000, 3999),
        (999, 9_223_372_036_854_775_807, 3_619_181_019_466_538_655),
        (
            9_223_372_036_854_775_804,
            9_223_372_036_854_775_807,
            3_074_457_345_618_258_602,
        ),
    ];
    
    #[test]
    fn invm_test() {
        // fixed cases
        for &(a, m, x) in INVM_CASES.iter() {
            assert_eq!(ModularOps::<&u64>::invm(&a, &m).unwrap(), x);
        }
    }
}
