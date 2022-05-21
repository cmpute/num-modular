use crate::{udouble, umax, ModularInteger, ModularUnaryOps};
use core::ops::*;
use num_traits::{Inv, Pow};

// FIXME: use unchecked operators to speed up calculation (after https://github.com/rust-lang/rust/issues/85122)
/// An unsigned integer modulo (pseudo) Mersenne primes `2^P - K`, it supports `P` up to 127 and `K < 2^(P-1)`
///
/// IMPORTANT NOTE: this class assumes that `2^P-K` is a prime. During compliation, we don't do full check
/// of the primality of `2^P-K`. If it's not a prime, then the modular division and inverse will panic.
#[derive(Clone, Copy, PartialEq, Eq)]
pub struct MersenneInt<const P: u8, const K: umax>(umax); // the underlying integer is in the half-open range [0, 2^P-K)

// XXX: support other primes as modulo, such as solinas prime, proth prime

impl<const P: u8, const K: umax> MersenneInt<P, K> {
    const BITMASK: umax = (1 << P) - 1;
    const MODULUS: umax = (1 << P) - K;

    // Calculate v % Self::MODULUS, where v is a umax integer
    const fn reduce_single(v: umax) -> umax {
        let mut lo = v & Self::BITMASK;
        let mut hi = v >> P;
        while hi > 0 {
            let sum = if K == 1 { hi + lo } else { hi * K + lo };
            lo = sum & Self::BITMASK;
            hi = sum >> P;
        }
    
        if K == 1 {
            lo
        } else {
            if lo >= Self::MODULUS {
                lo - Self::MODULUS
            } else {
                lo
            }
        }
    }

    // Calculate v % Self::MODULUS, where v is a udouble integer
    fn reduce_double(v: udouble) -> umax {
        // reduce modulo
        let mut lo = v.lo & Self::BITMASK;
        let mut hi = v >> P;
        while hi.hi > 0 {
            // first reduce until high bits fit in umax
            let sum = if K == 1 { hi + lo } else { hi * K + lo };
            lo = sum.lo & Self::BITMASK;
            hi = sum >> P;
        }

        let mut hi = hi.lo;
        while hi > 0 {
            // then reduce the smaller high bits
            let sum = if K == 1 { hi + lo } else { hi * K + lo };
            lo = sum & Self::BITMASK;
            hi = sum >> P;
        }

        if K == 1 {
            lo
        } else {
            if lo >= Self::MODULUS {
                lo - Self::MODULUS
            } else {
                lo
            }
        }
    }

    /// Create a new MersenneInt instance from a normal integer (by modulo `2^P-K`)
    #[inline]
    pub const fn new(n: umax) -> Self {
        // FIXME: use compile time checks, maybe after https://github.com/rust-lang/rust/issues/76560
        assert!(P <= 127);
        assert!(K > 0 && K < 2u128.pow(P as u32 - 1) && K % 2 == 1);
        assert!(
            Self::MODULUS % 3 != 0
                && Self::MODULUS % 5 != 0
                && Self::MODULUS % 7 != 0
                && Self::MODULUS % 11 != 0
                && Self::MODULUS % 13 != 0
        ); // error on easy composites
        Self(Self::reduce_single(n))
    }
}

impl<const P: u8, const K: umax> From<umax> for MersenneInt<P, K> {
    fn from(v: umax) -> Self {
        Self(v)
    }
}

impl<const P: u8, const K: umax> From<MersenneInt<P, K>> for umax {
    fn from(v: MersenneInt<P, K>) -> Self {
        v.0
    }
}

impl<const P: u8, const K: umax> Add for MersenneInt<P, K> {
    type Output = Self;
    #[inline]
    fn add(self, rhs: Self) -> Self {
        let sum = self.0 + rhs.0;
        Self(if sum >= Self::MODULUS {
            sum - Self::MODULUS
        } else {
            sum
        })
    }
}

impl<const P: u8, const K: umax> Sub for MersenneInt<P, K> {
    type Output = Self;
    #[inline]
    fn sub(self, rhs: Self) -> Self {
        Self(if self.0 >= rhs.0 {
            self.0 - rhs.0
        } else {
            Self::MODULUS - (rhs.0 - self.0)
        })
    }
}

impl<const P: u8, const K: umax> Mul for MersenneInt<P, K> {
    type Output = Self;
    #[inline]
    fn mul(self, rhs: Self) -> Self::Output {
        if (P as u32) < (umax::BITS / 2) {
            Self(Self::reduce_single(self.0 * rhs.0))
        } else {
            Self(Self::reduce_double(udouble::widening_mul(self.0, rhs.0)))
        }
    }
}

impl<const P: u8, const K: umax> Pow<umax> for MersenneInt<P, K> {
    type Output = Self;

    fn pow(self, rhs: umax) -> Self::Output {
        match rhs {
            1 => self,
            2 => self * self,
            _ => {
                let mut multi = self;
                let mut exp = rhs;
                let mut result = Self(1);
                while exp > 0 {
                    if exp & 1 != 0 {
                        result = result * multi;
                    }
                    multi = multi.square();
                    exp >>= 1;
                }
                result
            }
        }
    }
}

impl<const P: u8, const K: umax> Neg for MersenneInt<P, K> {
    type Output = Self;
    fn neg(self) -> Self::Output {
        Self(if self.0 == 0 {
            0
        } else {
            Self::MODULUS - self.0
        })
    }
}

impl<const P: u8, const K: umax> Inv for MersenneInt<P, K> {
    type Output = Self;
    fn inv(self) -> Self::Output {
        // It seems that extended gcd is faster than using fermat's theorem a^-1 = a^(p-2) mod p
        // For faster inverse using fermat theorem, refer to https://eprint.iacr.org/2018/1038.pdf (haven't benchmarked with this)
        Self(if (P as u32) < usize::BITS {
            (self.0 as usize)
                .invm(&(Self::MODULUS as usize))
                .expect("the modulus shoud be a prime") as umax
        } else {
            self.0
                .invm(&Self::MODULUS)
                .expect("the modulus shoud be a prime")
        })
    }
}

impl<const P: u8, const K: umax> Div for MersenneInt<P, K> {
    type Output = Self;
    fn div(self, rhs: Self) -> Self::Output {
        self * rhs.inv()
    }
}

impl<const P: u8, const K: umax> ModularInteger for MersenneInt<P, K> {
    type Base = umax;
    #[inline(always)]
    fn modulus(&self) -> Self::Base {
        Self::MODULUS
    }

    #[inline(always)]
    fn residue(&self) -> Self::Base {
        self.0
    }

    #[inline(always)]
    fn is_zero(&self) -> bool {
        self.0 == 0
    }

    #[inline(always)]
    fn convert(&self, n: Self::Base) -> Self {
        Self::new(n)
    }

    #[inline]
    fn double(self) -> Self {
        let sum = self.0 << 1;
        Self(if sum > Self::MODULUS {
            sum - Self::MODULUS
        } else {
            sum
        })
    }

    #[inline]
    fn square(self) -> Self {
        if (P as u32) < (umax::BITS / 2) {
            Self(Self::reduce_single(self.0 * self.0))
        } else {
            Self(Self::reduce_double(udouble::widening_square(self.0)))
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{ModularCoreOps, ModularPow};
    use rand::random;

    const P1: u128 = (1 << 31) - 1;
    const P2: u128 = (1 << 61) - 1;
    const P3: u128 = (1 << 127) - 1;
    const P4: u128 = (1 << 32) - 5;
    const P5: u128 = (1 << 56) - 5;
    const P6: u128 = (1 << 122) - 3;

    const NRANDOM: u32 = 10;

    #[test]
    fn creation_test() {
        // random creation test
        for _ in 0..NRANDOM {
            let a = random::<u128>();

            assert_eq!(MersenneInt::<31, 1>::new(a).residue(), a % P1);
            assert_eq!(MersenneInt::<61, 1>::new(a).residue(), a % P2);
            assert_eq!(MersenneInt::<127, 1>::new(a).residue(), a % P3);
            assert_eq!(MersenneInt::<32, 5>::new(a).residue(), a % P4);
            assert_eq!(MersenneInt::<56, 5>::new(a).residue(), a % P5);
            assert_eq!(MersenneInt::<122, 3>::new(a).residue(), a % P6);
        }
    }

    #[test]
    fn test_against_prim() {
        for _ in 0..NRANDOM {
            let (a, b) = (random::<u128>(), random::<u128>());
            let e = random::<u8>();

            // mod 2^31-1
            let am = MersenneInt::<31, 1>::new(a);
            let bm = MersenneInt::<31, 1>::new(b);
            assert_eq!((am + bm).residue(), a.addm(b, &P1));
            assert_eq!((am - bm).residue(), a.subm(b, &P1));
            assert_eq!((am * bm).residue(), a.mulm(b, &P1));
            assert_eq!((am / bm).residue(), a.mulm(b.invm(&P1).unwrap(), &P1));
            assert_eq!(am.neg().residue(), a.negm(&P1));
            assert_eq!(am.double().residue(), a.dblm(&P1));
            assert_eq!(am.square().residue(), a.sqm(&P1));
            assert_eq!(am.pow(e as u128).residue(), a.powm(e as u128, &P1));

            // mod 2^61-1
            let am = MersenneInt::<61, 1>::new(a);
            let bm = MersenneInt::<61, 1>::new(b);
            assert_eq!((am + bm).residue(), a.addm(b, &P2));
            assert_eq!((am - bm).residue(), a.subm(b, &P2));
            assert_eq!((am * bm).residue(), a.mulm(b, &P2));
            assert_eq!((am / bm).residue(), a.mulm(b.invm(&P2).unwrap(), &P2));
            assert_eq!(am.neg().residue(), a.negm(&P2));
            assert_eq!(am.double().residue(), a.dblm(&P2));
            assert_eq!(am.square().residue(), a.sqm(&P2));
            assert_eq!(am.pow(e as u128).residue(), a.powm(e as u128, &P2));

            // mod 2^127-1
            let am = MersenneInt::<127, 1>::new(a);
            let bm = MersenneInt::<127, 1>::new(b);
            assert_eq!((am + bm).residue(), a.addm(b, &P3));
            assert_eq!((am - bm).residue(), a.subm(b, &P3));
            assert_eq!((am * bm).residue(), a.mulm(b, &P3));
            assert_eq!((am / bm).residue(), a.mulm(b.invm(&P3).unwrap(), &P3));
            assert_eq!(am.neg().residue(), a.negm(&P3));
            assert_eq!(am.double().residue(), a.dblm(&P3));
            assert_eq!(am.square().residue(), a.sqm(&P3));
            assert_eq!(am.pow(e as u128).residue(), a.powm(e as u128, &P3));

            // mod 2^32-5
            let am = MersenneInt::<32, 5>::new(a);
            let bm = MersenneInt::<32, 5>::new(b);
            assert_eq!((am + bm).residue(), a.addm(b, &P4));
            assert_eq!((am - bm).residue(), a.subm(b, &P4));
            assert_eq!((am * bm).residue(), a.mulm(b, &P4));
            assert_eq!((am / bm).residue(), a.mulm(b.invm(&P4).unwrap(), &P4));
            assert_eq!(am.neg().residue(), a.negm(&P4));
            assert_eq!(am.double().residue(), a.dblm(&P4));
            assert_eq!(am.square().residue(), a.sqm(&P4));
            assert_eq!(am.pow(e as u128).residue(), a.powm(e as u128, &P4));

            // mod 2^56-5
            let am = MersenneInt::<56, 5>::new(a);
            let bm = MersenneInt::<56, 5>::new(b);
            assert_eq!((am + bm).residue(), a.addm(b, &P5));
            assert_eq!((am - bm).residue(), a.subm(b, &P5));
            assert_eq!((am * bm).residue(), a.mulm(b, &P5));
            assert_eq!((am / bm).residue(), a.mulm(b.invm(&P5).unwrap(), &P5));
            assert_eq!(am.neg().residue(), a.negm(&P5));
            assert_eq!(am.double().residue(), a.dblm(&P5));
            assert_eq!(am.square().residue(), a.sqm(&P5));
            assert_eq!(am.pow(e as u128).residue(), a.powm(e as u128, &P5));

            // mod 2^122-3
            let am = MersenneInt::<122, 3>::new(a);
            let bm = MersenneInt::<122, 3>::new(b);
            assert_eq!((am + bm).residue(), a.addm(b, &P6));
            assert_eq!((am - bm).residue(), a.subm(b, &P6));
            assert_eq!((am * bm).residue(), a.mulm(b, &P6));
            assert_eq!((am / bm).residue(), a.mulm(b.invm(&P6).unwrap(), &P6));
            assert_eq!(am.neg().residue(), a.negm(&P6));
            assert_eq!(am.double().residue(), a.dblm(&P6));
            assert_eq!(am.square().residue(), a.sqm(&P6));
            assert_eq!(am.pow(e as u128).residue(), a.powm(e as u128, &P6));
        }
    }
}
