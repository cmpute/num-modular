use num_integer::Integer;
use crate::ModularOps;

macro_rules! impl_jacobi_prim {
    ($T:ty) => {
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
                std::mem::swap(&mut a, &mut n);
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
    };
}

// implement inverse mod using extended euclidean algorithm
macro_rules! impl_invm_prim {
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
        impl ModularOps<$T, &$T> for $T {
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
            fn powm(self, exp: $T, m: &$T) -> $T {
                if exp == 1 {
                    return self % m;
                }
                if exp == 2 {
                    return self.mulm(self, m);
                }

                let mut multi = self % m;
                let mut exp = exp;
                let mut result = 1;
                while exp > 0 {
                    if exp & 1 > 0 {
                        result = result.mulm(multi, m);
                    }
                    multi = multi.mulm(multi, m);
                    exp >>= 1;
                }
                result
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
            impl_jacobi_prim!($T);
            impl_invm_prim!($T);
        }
    };
}

impl_mod_arithm_uu!(u8, u16);
impl_mod_arithm_uu!(u16, u32);
impl_mod_arithm_uu!(u32, u64);
impl_mod_arithm_uu!(u64, u128);
impl_mod_arithm_uu!(usize, u128);


impl ModularOps<u128, &u128> for u128 {
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

    fn powm(self, exp: u128, m: &u128) -> u128 {
        if exp == 1 {
            return self % m;
        }

        let mut multi = self % m;
        let mut exp = exp;
        let mut result = 1;
        while exp > 0 {
            if exp & 1 > 0 {
                result = result.mulm(multi, m);
            }
            multi = multi.mulm(multi, m);
            exp >>= 1;
        }
        result
    }

    fn negm(self, m: &u128) -> u128 {
        let x = self % m;
        if x == 0 {
            0
        } else {
            m - x
        }
    }

    impl_jacobi_prim!(u128);
    impl_invm_prim!(u128);
}

macro_rules! impl_mod_arithm_by_deref {
    ($($T:ty)*) => {$(
        impl ModularOps<$T, &$T> for &$T {
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
            fn powm(self, exp: $T, m: &$T) -> $T {
                (*self).powm(exp, &m)
            }
            #[inline]
            fn negm(self, m: &$T) -> $T {
                ModularOps::<$T, &$T>::negm(*self, m)
            }
            #[inline]
            fn invm(self, m: &$T) -> Option<$T> {
                ModularOps::<$T, &$T>::invm(*self, m)
            }
            #[inline]
            fn jacobi(self, n: &$T) -> i8 {
                ModularOps::<$T, &$T>::jacobi(*self, n)
            }
        }

        impl ModularOps<&$T, &$T> for $T {
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
            fn powm(self, exp: &$T, m: &$T) -> $T {
                self.powm(*exp, &m)
            }
            #[inline]
            fn negm(self, m: &$T) -> $T {
                ModularOps::<$T, &$T>::negm(self, m)
            }
            #[inline]
            fn invm(self, m: &$T) -> Option<$T> {
                ModularOps::<$T, &$T>::invm(self, m)
            }
            #[inline]
            fn jacobi(self, n: &$T) -> i8 {
                ModularOps::<$T, &$T>::jacobi(self, n)
            }
        }

        impl ModularOps<&$T, &$T> for &$T {
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
            fn powm(self, exp: &$T, m: &$T) -> $T {
                (*self).powm(*exp, &m)
            }
            #[inline]
            fn negm(self, m: &$T) -> $T {
                ModularOps::<$T, &$T>::negm(*self, m)
            }
            #[inline]
            fn invm(self, m: &$T) -> Option<$T> {
                ModularOps::<$T, &$T>::invm(*self, m)
            }
            #[inline]
            fn jacobi(self, n: &$T) -> i8 {
                ModularOps::<$T, &$T>::jacobi(*self, n)
            }
        }
    )*};
}

impl_mod_arithm_by_deref!(u8 u16 u32 u64 u128 usize);
