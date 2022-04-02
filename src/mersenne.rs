use crate::{ModularInteger, udouble};
use core::ops::*;

// TODO: static_assert check P <= 127, K < 2^(P-1)
/// An unsigned integer modulo `2^P - K`, it supports P up to 127 and K < 2^(P-1)
#[derive(Clone, Copy, PartialEq, Eq)]
pub struct MersenneInt<const P: u8, const K: u128>(u128);

impl<const P: u8, const K: u128> MersenneInt<P, K> {
    const BITMASK: u128 = (1 << P) - 1;
    const MODULUS: u128 = (1 << P) - K;

    #[inline]
    pub const fn new(n: u128) -> Self {
        let mut lo = n & Self::BITMASK;
        let mut hi = n >> P;
        while hi > 0 {
            let sum = if K == 1 { lo + hi } else { lo + hi * K };
            lo = sum & Self::BITMASK;
            hi = sum >> P;
        }

        let v = if K == 1 {
            if lo == Self::BITMASK { 0 } else { lo }
        } else {
            if lo >= Self::MODULUS { lo - Self::MODULUS } else { lo }
        };
        Self(v)
    }
}

impl<const P: u8, const K: u128> From<u128> for MersenneInt<P, K> {
    fn from(v: u128) -> Self {
        Self(v)
    }
}

impl<const P: u8, const K: u128> From<MersenneInt<P, K>> for u128 {
    fn from(v: MersenneInt<P, K>) -> Self {
        v.0
    }
}

impl<const P: u8, const K: u128> Add for MersenneInt<P, K> {
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

impl<const P: u8, const K: u128> Sub for MersenneInt<P, K> {
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

impl<const P: u8, const K: u128> Mul for MersenneInt<P, K> {
    type Output = Self;
    fn mul(self, rhs: Self) -> Self::Output {
        let prod = udouble::widening_mul(self.0, rhs.0);
        
        // reduce modulo
        let mut lo = prod.lo & Self::BITMASK;
        let mut hi = prod >> P.into();
        while hi.hi > 0 { // first reduce until high bits fit in u128
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
            if lo == Self::BITMASK { 0 } else { lo }
        } else {
            if lo >= Self::MODULUS { lo - Self::MODULUS } else { lo }
        };
        Self(v)
    }
}

// TODO: implement inverse and division using a^-1 = a^(p-2) mod p, assuming p is prime
// We can explicitly check cases for K < 3 using static assertion

impl<const P: u8, const K: u128> Neg for MersenneInt<P, K> {
    type Output = Self;
    fn neg(self) -> Self::Output {
        Self(if self.0 == 0 {
            0
        } else {
            Self::MODULUS - self.0
        })
    }
}

impl<const P: u8, const K: u128> ModularInteger for MersenneInt<P, K> {
    type Base = u128;
    #[inline]
    fn modulus(&self) -> &Self::Base {
        &Self::MODULUS
    }

    #[inline]
    fn residue(&self) -> Self::Base {
        self.0
    }

    #[inline]
    fn new(&self, n: Self::Base) -> Self {
        Self::new(n)
    }
}
