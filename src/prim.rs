use crate::{ModularCoreOps, ModularOps};
use num_integer::Integer;

// TODO: implement the modular functions as const: https://github.com/rust-lang/rust/pull/68847
// TODO: provide utility functions to convert signed integers to unsigned during modular operation
// (especially negm and absm)

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

    // TODO: benchmark against http://www.janfeitsma.nl/math/psp2/expmod
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
