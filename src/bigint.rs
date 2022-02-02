use crate::ModularOps;

#[cfg(feature="num-bigint")]
use num_traits::{Zero, One, ToPrimitive};
#[cfg(feature="num-bigint")]
use num_bigint::BigUint;

#[cfg(feature="num-bigint")]
impl ModularOps<&BigUint, &BigUint> for &BigUint {
    type Output = BigUint;

    #[inline]
    fn addm(self, rhs: &BigUint, m: &BigUint) -> BigUint {
        (self + rhs) % m
    }
    fn subm(self, rhs: &BigUint, m: &BigUint) -> BigUint {
        let (lhs, rhs) = (self % m, rhs % m);
        if lhs >= rhs {
            lhs - rhs
        } else {
            m - (rhs - lhs)
        }
    }

    fn mulm(self, rhs: &BigUint, m: &BigUint) -> BigUint {
        let a = self % m;
        let b = rhs % m;

        if let Some(sm) = m.to_u64() {
            let sself = a.to_u64().unwrap();
            let srhs = b.to_u64().unwrap();
            return BigUint::from(sself.mulm(srhs, &sm));
        }

        (a * b) % m
    }

    #[inline]
    fn powm(self, exp: &BigUint, m: &BigUint) -> BigUint {
        self.modpow(&exp, m)
    }

    #[inline]
    fn negm(self, m: &BigUint) -> BigUint {
        let x = self % m;
        if x.is_zero() {
            BigUint::zero()
        } else {
            m - x
        }
    }

    fn jacobi(self, n: &BigUint) -> i8 {
        debug_assert!(n.is_odd());

        if self.is_zero() {
            return 0;
        }
        if self.is_one() {
            return 1;
        }

        let three = BigUint::from(3u8);
        let five = BigUint::from(5u8);
        let seven = BigUint::from(7u8);

        let mut a = self % n;
        let mut n = n.clone();
        let mut t = 1;
        while a > BigUint::zero() {
            while a.is_even() {
                a >>= 1;
                if &n & &seven == three || &n & &seven == five {
                    t *= -1;
                }
            }
            std::mem::swap(&mut a, &mut n);
            if (&a & &three) == three && (&n & &three) == three {
                t *= -1;
            }
            a %= &n;
        }
        if n.is_one() {
            t
        } else {
            0
        }
    }

    fn invm(self, m: &BigUint) -> Option<Self::Output> {
        let x = if self >= m { self % m } else { self.clone() };

        let (mut last_r, mut r) = (m.clone(), x);
        let (mut last_t, mut t) = (BigUint::zero(), BigUint::one());

        while r > BigUint::zero() {
            let (quo, rem) = last_r.div_rem(&r);
            last_r = r;
            r = rem;

            let new_t = last_t.subm(&quo.mulm(&t, m), m);
            last_t = t;
            t = new_t;
        }

        // if r = gcd(self, m) > 1, then inverse doesn't exist
        if last_r > BigUint::one() {
            None
        } else {
            Some(last_t)
        }
    }
}

#[cfg(feature="num-bigint")]
macro_rules! impl_mod_arithm_by_ref {
    ($T:ty) => {
        impl ModularOps<$T, &$T> for &$T {
            type Output = $T;
            #[inline]
            fn addm(self, rhs: $T, m: &$T) -> $T {
                self.addm(&rhs, &m)
            }
            #[inline]
            fn subm(self, rhs: $T, m: &$T) -> $T {
                self.subm(&rhs, &m)
            }
            #[inline]
            fn mulm(self, rhs: $T, m: &$T) -> $T {
                self.mulm(&rhs, &m)
            }
            #[inline]
            fn powm(self, exp: $T, m: &$T) -> $T {
                self.powm(&exp, &m)
            }
            #[inline]
            fn negm(self, m: &$T) -> $T {
                ModularOps::<&$T, &$T>::negm(self, m)
            }
            #[inline]
            fn invm(self, m: &$T) -> Option<$T> {
                ModularOps::<&$T, &$T>::invm(self, m)
            }
            #[inline]
            fn jacobi(self, n: &$T) -> i8 {
                ModularOps::<&$T, &$T>::jacobi(self, n)
            }
        }

        impl ModularOps<&$T, &$T> for $T {
            type Output = $T;
            #[inline]
            fn addm(self, rhs: &$T, m: &$T) -> $T {
                (&self).addm(rhs, &m)
            }
            #[inline]
            fn subm(self, rhs: &$T, m: &$T) -> $T {
                (&self).subm(rhs, &m)
            }
            #[inline]
            fn mulm(self, rhs: &$T, m: &$T) -> $T {
                (&self).mulm(rhs, &m)
            }
            #[inline]
            fn powm(self, exp: &$T, m: &$T) -> $T {
                (&self).powm(exp, &m)
            }
            #[inline]
            fn negm(self, m: &$T) -> $T {
                ModularOps::<&$T, &$T>::negm(&self, m)
            }
            #[inline]
            fn invm(self, m: &$T) -> Option<$T> {
                ModularOps::<&$T, &$T>::invm(&self, m)
            }
            #[inline]
            fn jacobi(self, n: &$T) -> i8 {
                ModularOps::<&$T, &$T>::jacobi(&self, n)
            }
        }

        impl ModularOps<$T, &$T> for $T {
            type Output = $T;
            #[inline]
            fn addm(self, rhs: $T, m: &$T) -> $T {
                (&self).addm(&rhs, &m)
            }
            #[inline]
            fn subm(self, rhs: $T, m: &$T) -> $T {
                (&self).subm(&rhs, &m)
            }
            #[inline]
            fn mulm(self, rhs: $T, m: &$T) -> $T {
                (&self).mulm(&rhs, &m)
            }
            #[inline]
            fn powm(self, exp: $T, m: &$T) -> $T {
                (&self).powm(&exp, &m)
            }
            #[inline]
            fn negm(self, m: &$T) -> $T {
                ModularOps::<&$T, &$T>::negm(&self, m)
            }
            #[inline]
            fn invm(self, m: &$T) -> Option<$T> {
                ModularOps::<&$T, &$T>::invm(&self, m)
            }
            #[inline]
            fn jacobi(self, n: &$T) -> i8 {
                ModularOps::<&$T, &$T>::jacobi(&self, n)
            }
        }
    };
}

#[cfg(feature="num-bigint")]
impl_mod_arithm_by_ref!(BigUint);
