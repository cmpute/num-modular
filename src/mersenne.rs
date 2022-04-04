use crate::{ModularInteger, udouble, umax};
use core::ops::*;
use num_traits::Pow;

// TODO(v0.2.2): static_assert check P <= 127, K < 2^(P-1)
// TODO: use unchecked operators to speed up calculation
/// An unsigned integer modulo (pseudo) Mersenne primes `2^P-K`, it supports P up to 127 and K < 2^(P-1)
/// 
/// IMPORTANT NOTE: this class assumes that `2^P-K` is a prime. We don't do compile time check
/// on the primality of `2^P-K`. If it's not a prime, then the modular division and inverse
/// will panic.
#[derive(Clone, Copy, PartialEq, Eq)]
pub struct MersenneInt<const P: u8, const K: umax>(umax); // the underlying integer is in the inclusive range [0, 2^P-K]

impl<const P: u8, const K: umax> MersenneInt<P, K> {
    const BITMASK: umax = (1 << P) - 1;
    const MODULUS: umax = (1 << P) - K;

    #[inline]
    pub const fn new(n: umax) -> Self {
        let mut lo = n & Self::BITMASK;
        let mut hi = n >> P;
        while hi > 0 {
            let sum = if K == 1 { lo + hi } else { lo + hi * K };
            lo = sum & Self::BITMASK;
            hi = sum >> P;
        }

        let v = if K == 1 {
            lo
        } else {
            if lo > Self::MODULUS { lo - Self::MODULUS } else { lo }
        };
        Self(v)
    }
}

impl<const P: u8, const K: umax> From<umax> for MersenneInt<P, K> {
    fn from(v: umax) -> Self {
        Self(v)
    }
}

impl<const P: u8, const K: umax> From<i128> for MersenneInt<P, K> {
    fn from(v: i128) -> Self {
        unimplemented!() // TODO: implement
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
        Self(if sum > Self::MODULUS {
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
    fn mul(self, rhs: Self) -> Self::Output {
        let prod = udouble::widening_mul(self.0, rhs.0);
        
        // reduce modulo
        let mut lo = prod.lo & Self::BITMASK;
        let mut hi = prod >> P.into();
        while hi.hi > 0 { // first reduce until high bits fit in umax
            let sum = if K == 1 { hi + lo } else { hi * K + lo };
            lo = sum.lo & Self::BITMASK;
            hi = sum >> P.into();
        }

        let mut hi = hi.lo;
        while hi > 0 { // then reduce the smaller high bits
            let sum = if K == 1 { hi + lo } else { hi * K + lo };
            lo = sum & Self::BITMASK;
            hi = sum >> P;
        }

        let v = if K == 1 {
            lo
        } else {
            if lo > Self::MODULUS { lo - Self::MODULUS } else { lo }
        };
        Self(v)
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
                    multi = multi * multi;
                    exp >>= 1;
                }
                result
            }
        }
    }
}

// TODO: implement inverse and division using a^-1 = a^(p-2) mod p, assuming p is prime
// We can explicitly check cases for K < 3 using static assertion
// It seems that extended gcd is faster
// To improve inverse speed, reference to https://eprint.iacr.org/2018/1038.pdf

impl<const P: u8, const K: umax> Neg for MersenneInt<P, K> {
    type Output = Self;
    fn neg(self) -> Self::Output {
        Self(Self::MODULUS - self.0)
    }
}

impl<const P: u8, const K: umax> ModularInteger for MersenneInt<P, K> {
    type Base = umax;
    #[inline]
    fn modulus(&self) -> &Self::Base {
        &Self::MODULUS
    }

    #[inline]
    fn residue(&self) -> Self::Base {
        if self.0 == Self::MODULUS {
            0
        } else {
            self.0
        }
    }

    #[inline]
    fn new(&self, n: Self::Base) -> Self {
        Self::new(n)
    }
}

#[cfg(test)]
mod tests {
    // TODO: add tests
}
