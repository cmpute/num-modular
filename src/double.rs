//! This module implements a double width integer type based on the largest built-in integer (u128)

use core::ops::{Add, Rem, Shl, Shr, ShlAssign, ShrAssign, AddAssign, Sub, SubAssign};
use num_traits::Zero;
use num_traits::ops::overflowing::OverflowingAdd;

/// Alias of the builtin integer type with max width
#[allow(non_camel_case_types)]
pub type umax = u128;

#[allow(non_camel_case_types)]
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
/// A double width integer type based on the largest built-in integer (u128)
/// 
/// It's used to support double-width operations on u128. Although it can be considered as u256,
/// it's not as feature-rich as other crates since it's only designed to support this crate.
pub struct udouble {
    pub hi: umax,
    pub lo: umax
}

impl udouble {
    pub fn widening_add(lhs: umax, rhs: umax) -> Self {
        let (sum, carry) = lhs.overflowing_add(rhs);
        udouble { hi: carry as umax, lo: sum }
    }

    /// Calculate multiplication of two [umax] integers with result represented in double width integer
    // equivalent to umul_ppmm
    pub fn widening_mul(lhs: umax, rhs: umax) -> Self {
        const HALF_BITS: u32 = umax::BITS >> 1;
        let halves = |x| (x >> HALF_BITS, x & !(0 as umax) >> HALF_BITS);
        let ((x1, x0), (y1, y0)) = (halves(lhs), halves(rhs));
        let (z2, z0) = (x1 * y1, x0 * y0);
        let (z1, c1) = (x1 * y0).overflowing_add(x0 * y1);
        let (lo, c0) = umax::overflowing_add(z0, z1 << HALF_BITS);
        Self { hi: z2 + (z1 >> HALF_BITS) + c0 as umax + ((c1 as umax) << HALF_BITS), lo }
    }
}

impl From<umax> for udouble {
    fn from(v: umax) -> Self {
        Self { lo: v, hi: 0 }
    }
}

impl Add for udouble {
    type Output = udouble;

    // equivalent to add_ssaaaa
    fn add(self, rhs: Self) -> Self::Output {
        let (lo, carry) = self.lo.overflowing_add(rhs.lo);
        let hi = self.hi + rhs.hi + carry as umax;
        Self { lo, hi }
    }
}

impl Add<umax> for udouble {
    type Output = udouble;
    fn add(self, rhs: umax) -> Self::Output {
        let (lo, carry) = self.lo.overflowing_add(rhs);
        let hi = if carry { self.hi + 1 } else { self.hi };
        Self { lo, hi }
    }
}

impl AddAssign for udouble {
    fn add_assign(&mut self, rhs: Self) {
        let (lo, carry) = self.lo.overflowing_add(rhs.lo);
        self.lo = lo;
        self.hi += rhs.hi + carry as umax;
    }
}

impl AddAssign<umax> for udouble {
    fn add_assign(&mut self, rhs: umax) {
        let (lo, carry) = self.lo.overflowing_add(rhs);
        self.lo = lo;
        if carry { self.hi += 1 }
    }
}

impl OverflowingAdd for udouble {
    fn overflowing_add(&self, rhs: &Self) -> (Self, bool) {
        let (lo, carry) = self.lo.overflowing_add(rhs.lo);
        let (hi, of1) = self.hi.overflowing_add(rhs.hi);
        let (hi, of2) = hi.overflowing_add(carry as umax);
        (Self { lo, hi }, of1 || of2)
    }
}

impl Sub for udouble {
    type Output = Self;
    fn sub(self, rhs: Self) -> Self::Output {
        let carry = self.lo < rhs.lo;
        let lo = self.lo.wrapping_sub(rhs.lo);
        let hi = self.hi - rhs.hi - carry as umax;
        Self { lo, hi }
    }
}

impl SubAssign for udouble {
    fn sub_assign(&mut self, rhs: Self) {
        let carry = self.lo < rhs.lo;
        self.lo = self.lo.wrapping_sub(rhs.lo);
        self.hi -= rhs.hi + carry as umax;
    }
}

impl Zero for udouble {
    fn zero() -> Self {
        Self { lo: 0, hi: 0}
    }

    fn is_zero(&self) -> bool {
        self.lo == 0 && self.hi == 0
    }
}

impl Shl<u32> for udouble {
    type Output = Self;
    fn shl(self, rhs: u32) -> Self::Output {
        match rhs {
            0 => self,
            s if s >= umax::BITS => Self {
                hi: self.lo << (s - umax::BITS),
                lo: 0,
            },
            s => Self {
                lo: self.lo << s,
                hi: (self.hi << s) | (self.lo >> (umax::BITS - s)),
            }
        }
    }
}

impl ShlAssign<u32> for udouble {
    fn shl_assign(&mut self, rhs: u32) {
        match rhs {
            0 => {},
            s if s >= umax::BITS => {
                self.hi = self.lo << (s - umax::BITS);
                self.lo = 0;
            },
            s => {
                self.hi <<= s;
                self.hi |= self.lo >> (umax::BITS - s);
                self.lo <<= s;
            }
        }
    }
}

impl Shr<u32> for udouble {
    type Output = Self;
    fn shr(self, rhs: u32) -> Self::Output {
        match rhs {
            0 => self,
            s if s >= umax::BITS => Self {
                lo: self.hi >> (rhs - umax::BITS),
                hi: 0,
            },
            s => Self {
                hi: self.hi >> s,
                lo: (self.lo >> s) | (self.hi << (umax::BITS - s)),
            }
        }
    }
}

impl ShrAssign<u32> for udouble {
    fn shr_assign(&mut self, rhs: u32) {
        match rhs {
            0 => {},
            s if s >= umax::BITS => {
                self.lo = self.hi >> (rhs - umax::BITS);
                self.hi = 0;
            },
            s => {
                self.lo >>= s;
                self.lo |= self.hi << (umax::BITS - s);
                self.hi >>= s;
            }
        }
    }
}

impl udouble {
    pub fn leading_zeros(self) -> u32 {
        if self.hi == 0 {
            self.lo.leading_zeros() + umax::BITS
        } else {
            self.hi.leading_zeros()
        }
    }

    // REF:
    // https://docs.rs/u256/0.1.0/src/u256/lib.rs.html#165-199
    // https://github.com/coreutils/coreutils/blob/master/src/factor.c#L284-L304
    // similar to udiv_qrnnd
    pub fn div_rem(self, other: Self) -> (Self, Self) {
        let mut n = self; // numerator
        let mut d = other; // denominator
        let mut q = Self::zero(); // quotient

        let nbits = 2 * umax::BITS - n.leading_zeros();
        let dbits = 2 * umax::BITS - d.leading_zeros();

        // Check for division by 0
        assert!(dbits != 0);

        // Early return in case we are dividing by a larger number than us
        if nbits < dbits {
            return (q, n);
        }

        // Bitwise long division
        let mut shift = nbits - dbits;
        d <<= shift;
        while shift > 0 {
            if n >= d {
                q += 1;
                n -= d;
            }
            d >>= 1;
            q <<= 1;
            shift -= 1;
        }
        (q, n)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::random;

    #[test]
    fn test_construction() {
        // construct zero
        assert!(udouble { lo: 0, hi: 0 }.is_zero());

        // from widening operators
        assert_eq!(udouble { hi: 0, lo: 2 }, udouble::widening_add(1, 1));
        assert_eq!(udouble { hi: 1, lo: umax::MAX - 1 }, udouble::widening_add(umax::MAX, umax::MAX));

        assert_eq!(udouble { hi: 0, lo: 1 }, udouble::widening_mul(1, 1));
        assert_eq!(udouble { hi: 1 << 32, lo: 0 }, udouble::widening_mul(1 << 80, 1 << 80));
        assert_eq!(udouble { hi: 1 << 32, lo: 2 << 120 | 1 << 80 }, udouble::widening_mul(1 << 80 | 1 << 40, 1 << 80 | 1 << 40));
        assert_eq!(udouble { hi: umax::MAX - 1, lo: 1 }, udouble::widening_mul(umax::MAX, umax::MAX));
    }

    #[test]
    fn test_ops() {
        const ONE: udouble = udouble { hi: 0, lo: 1 };
        const MAX: udouble = udouble { hi: 0, lo: umax::MAX };
        const ONEZERO: udouble = udouble { hi: 1, lo: 0 };
        const ONEMAX: udouble = udouble { hi: 1, lo: umax::MAX };
        const TWOZERO: udouble = udouble { hi: 2, lo: 0 };

        assert_eq!(ONE + MAX, ONEZERO);
        assert_eq!(ONE + ONEMAX, TWOZERO);
        assert_eq!(ONEZERO - ONE, MAX);
        assert_eq!(ONEZERO - MAX, ONE);
        assert_eq!(TWOZERO - ONE, ONEMAX);
        assert_eq!(TWOZERO - ONEMAX, ONE);

        assert_eq!(ONE << umax::BITS, ONEZERO);
        assert_eq!((MAX << 1) + 1, ONEMAX);
        assert_eq!(ONEZERO >> umax::BITS, ONE);
        assert_eq!(ONEMAX >> 1, MAX);
    }

    #[test]
    fn test_assign_ops() {
        for _ in 0..10 {
            let x = udouble { hi: random::<u32>() as umax, lo: random() };
            let y = udouble { hi: random::<u32>() as umax, lo: random() };
            let mut z = x;

            z += y; assert_eq!(z, x + y);
            z -= y; assert_eq!(z, x);
        }
    }
}
