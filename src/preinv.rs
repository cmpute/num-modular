use crate::{DivExact, ModularUnaryOps};

/// Pre-computation for fast divisibility check.
///
/// This struct stores the modular inverse of a divisor, and a limit for divisibility check.
/// See <https://math.stackexchange.com/a/1251328> for the explanation of the trick
#[derive(Debug, Clone, Copy)]
pub struct PreInv<T> {
    d_inv: T, // modular inverse of divisor
    q_lim: T, // limit of residue
}

macro_rules! impl_preinv_for_prim_int {
    ($t:ty, $tdouble:ty) => {
        impl PreInv<$t> {
            /// Construct the preinv instance with raw values.
            ///
            /// This function can be used to initialize preinv in a constant context, the divisor d
            /// is required only for verification of d_inv and q_lim.
            #[inline]
            pub const fn new(d_inv: $t, q_lim: $t) -> Self {
                Self { d_inv, q_lim }
            }

            // check if the divisor is consistent in debug mode
            #[inline]
            fn debug_check(&self, d: $t) {
                debug_assert!(d % 2 != 0, "only odd divisors are supported");
                debug_assert!(d.wrapping_mul(self.d_inv) == 1);
                debug_assert!(self.q_lim * d > (<$t>::MAX - d));
            }
        }

        impl From<$t> for PreInv<$t> {
            #[inline]
            fn from(v: $t) -> Self {
                debug_assert!(v % 2 != 0, "only odd divisors are supported");
                let d_inv = (v as $tdouble)
                    .invm(&((1 as $tdouble) << <$t>::BITS))
                    .unwrap() as $t;
                let q_lim = <$t>::MAX / v;
                Self { d_inv, q_lim }
            }
        }

        impl DivExact<$t, PreInv<$t>> for $t {
            type Output = $t;
            #[inline]
            fn div_exact(self, d: $t, pre: &PreInv<$t>) -> Option<Self> {
                pre.debug_check(d);
                let q = self.wrapping_mul(pre.d_inv);
                if q <= pre.q_lim {
                    Some(q)
                } else {
                    None
                }
            }
        }

        impl DivExact<$t, PreInv<$t>> for $tdouble {
            type Output = $tdouble;
            #[inline]
            fn div_exact(self, d: $t, pre: &PreInv<$t>) -> Option<$tdouble> {
                pre.debug_check(d);

                // this implementation comes from GNU factor,
                // see https://math.stackexchange.com/q/4436380/815652 for explanation

                let (n1, n0) = ((self >> <$t>::BITS) as $t, self as $t);
                let q0 = n0.wrapping_mul(pre.d_inv);
                let nr0 = (q0 as $tdouble) * (d as $tdouble);
                let nr0 = (nr0 >> <$t>::BITS) as $t;
                if nr0 > n1 {
                    return None;
                }
                let nr1 = n1 - nr0;
                let q1 = nr1.wrapping_mul(pre.d_inv);
                if q1 > pre.q_lim {
                    return None;
                }
                Some(((q1 as $tdouble) << <$t>::BITS) + q0 as $tdouble)
            }
        }
    };
}
impl_preinv_for_prim_int!(u8, u16);
impl_preinv_for_prim_int!(u16, u32);
impl_preinv_for_prim_int!(u32, u64);
impl_preinv_for_prim_int!(u64, u128);

// XXX: we could implement the div_exact for big integers, in a similar way to support double width
// REF: https://gmplib.org/manual/Exact-Division

// XXX: unchecked div_exact can be introduced by not checking the q_lim,
//      investigate this after `exact_div` is introduced or removed from core lib
//      https://github.com/rust-lang/rust/issues/85122
// REF: GMP `mpz_divexact`, `mpn_divexact`

#[cfg(test)]
mod tests {
    use super::*;
    use rand::random;

    #[test]
    fn div_exact_test() {
        const N: u8 = 100;
        for _ in 0..N {
            // u8 test
            let d = random::<u8>() | 1;
            let pre: PreInv<_> = d.into();

            let n: u8 = random();
            let expect = if n % d == 0 { Some(n / d) } else { None };
            assert_eq!(n.div_exact(d, &pre), expect, "{} / {}", n, d);
            let n: u16 = random();
            let expect = if n % (d as u16) == 0 {
                Some(n / (d as u16))
            } else {
                None
            };
            assert_eq!(n.div_exact(d, &pre), expect, "{} / {}", n, d);

            // u16 test
            let d = random::<u16>() | 1;
            let pre: PreInv<_> = d.into();

            let n: u16 = random();
            let expect = if n % d == 0 { Some(n / d) } else { None };
            assert_eq!(n.div_exact(d, &pre), expect, "{} / {}", n, d);
            let n: u32 = random();
            let expect = if n % (d as u32) == 0 {
                Some(n / (d as u32))
            } else {
                None
            };
            assert_eq!(n.div_exact(d, &pre), expect, "{} / {}", n, d);

            // u32 test
            let d = random::<u32>() | 1;
            let pre: PreInv<_> = d.into();

            let n: u32 = random();
            let expect = if n % d == 0 { Some(n / d) } else { None };
            assert_eq!(n.div_exact(d, &pre), expect, "{} / {}", n, d);
            let n: u64 = random();
            let expect = if n % (d as u64) == 0 {
                Some(n / (d as u64))
            } else {
                None
            };
            assert_eq!(n.div_exact(d, &pre), expect, "{} / {}", n, d);

            // u64 test
            let d = random::<u64>() | 1;
            let pre: PreInv<_> = d.into();

            let n: u64 = random();
            let expect = if n % d == 0 { Some(n / d) } else { None };
            assert_eq!(n.div_exact(d, &pre), expect, "{} / {}", n, d);
            let n: u128 = random();
            let expect = if n % (d as u128) == 0 {
                Some(n / (d as u128))
            } else {
                None
            };
            assert_eq!(n.div_exact(d, &pre), expect, "{} / {}", n, d);
        }
    }
}
