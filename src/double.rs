//! This module implements a double width integer type based on the largest built-in integer (u128)

use core::ops::*;
use num_traits::{One, Zero};

/// Alias of the builtin integer type with max width
#[allow(non_camel_case_types)]
pub type umax = u128;

#[allow(non_camel_case_types)]
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
/// A double width integer type based on the largest built-in integer (currently u128) and
/// to support double-width operations on it is the only goal for this type.
///
/// Although it can be regarded as u256, it's not as feature-rich as in other crates
/// since it's only designed to support this crate and few other crates (will be noted in comments).
pub struct udouble {
    /// Most significant part
    pub hi: umax,
    /// Least significant part
    pub lo: umax,
}

// TODO: port optimizations from and benchmark against
//       https://github.com/nlordell/ethnum-rs
//       https://github.com/0xProject/OpenZKP/tree/master/algebra/u256
//       https://github.com/paritytech/parity-common/blob/master/uint
// TODO(v0.2.2): further benchmark against
//       https://docs.rs/bigint/latest/bigint/uint/struct.U256.html
//       https://trussed-dev.github.io/trussed/crypto_bigint/type.U256.html
//       https://crates.io/crates/crypto-bigint
// TODO(v0.3): remove excessive APIs that are not (likely) used (to reduce optimization burden)

impl udouble {
    pub fn widening_add(lhs: umax, rhs: umax) -> Self {
        let (sum, carry) = lhs.overflowing_add(rhs);
        udouble {
            hi: carry as umax,
            lo: sum,
        }
    }

    /// Calculate multiplication of two [umax] integers with result represented in double width integer
    // equivalent to umul_ppmm, can be implemented efficiently with carrying_mul and widening_mul implemented (rust#85532)
    pub fn widening_mul(lhs: umax, rhs: umax) -> Self {
        const HALF_BITS: u32 = umax::BITS >> 1;
        let halves = |x| (x >> HALF_BITS, x & !(0 as umax) >> HALF_BITS);
        let ((x1, x0), (y1, y0)) = (halves(lhs), halves(rhs));
        let (z2, z0) = (x1 * y1, x0 * y0);
        // it's possible to use Karatsuba multiplication, but overflow checking is much easier here
        let (z1, c1) = (x1 * y0).overflowing_add(x0 * y1);
        let (lo, c0) = umax::overflowing_add(z0, z1 << HALF_BITS);
        Self {
            hi: z2 + (z1 >> HALF_BITS) + c0 as umax + ((c1 as umax) << HALF_BITS),
            lo,
        }
    }

    pub fn overflowing_add(&self, rhs: Self) -> (Self, bool) {
        let (lo, carry) = self.lo.overflowing_add(rhs.lo);
        let (hi, of1) = self.hi.overflowing_add(rhs.hi);
        let (hi, of2) = hi.overflowing_add(carry as umax);
        (Self { lo, hi }, of1 || of2)
    }

    pub fn overflowing_mul(&self, rhs: Self) -> (Self, bool) {
        let c2 = self.hi != 0 && rhs.hi != 0;
        let Self { lo: z0, hi: c0 } = Self::widening_mul(self.lo, rhs.lo);
        let (z1x, c1x) = u128::overflowing_mul(self.lo, rhs.hi);
        let (z1y, c1y) = u128::overflowing_mul(self.hi, rhs.lo);
        let (z1z, c1z) = u128::overflowing_add(z1x, z1y);
        let (z1, c1) = z1z.overflowing_add(c0);
        (Self { hi: z1, lo: z0 }, c1x | c1y | c1z | c1 | c2)
    }

    /// Multiplication of double width and single width
    pub fn overflowing_muls(&self, rhs: umax) -> (Self, bool) {
        let Self { lo: z0, hi: c0 } = Self::widening_mul(self.lo, rhs);
        let (z1, c1) = u128::overflowing_mul(self.hi, rhs);
        let (z1, cs1) = z1.overflowing_add(c0);
        (Self { hi: z1, lo: z0 }, c1 | c1 | cs1)
    }

    pub fn checked_shl(self, rhs: u32) -> Option<Self> {
        if rhs < umax::BITS * 2 {
            Some(self << rhs)
        } else {
            None
        }
    }

    pub fn checked_shr(self, rhs: u32) -> Option<Self> {
        if rhs < umax::BITS * 2 {
            Some(self >> rhs)
        } else {
            None
        }
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
        if carry {
            self.hi += 1
        }
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
        Self { lo: 0, hi: 0 }
    }

    fn is_zero(&self) -> bool {
        self.lo == 0 && self.hi == 0
    }
}

impl One for udouble {
    fn one() -> Self {
        Self { lo: 1, hi: 0 }
    }
}

// TODO: support only Shr/Shl<u8> and <u16>

impl Shl<u8> for udouble {
    type Output = Self;
    fn shl(self, rhs: u8) -> Self::Output {
        match rhs {
            0 => self,
            s if s >= umax::BITS as u8 => Self {
                hi: self.lo << (s - umax::BITS as u8),
                lo: 0,
            },
            s => Self {
                lo: self.lo << s,
                hi: (self.hi << s) | (self.lo >> (umax::BITS as u8 - s)),
            },
        }
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
            },
        }
    }
}

impl ShlAssign<u8> for udouble {
    fn shl_assign(&mut self, rhs: u8) {
        match rhs {
            0 => {}
            s if s >= umax::BITS as u8 => {
                self.hi = self.lo << (s - umax::BITS as u8);
                self.lo = 0;
            }
            s => {
                self.hi <<= s;
                self.hi |= self.lo >> (umax::BITS as u8 - s);
                self.lo <<= s;
            }
        }
    }
}

impl ShlAssign<u32> for udouble {
    fn shl_assign(&mut self, rhs: u32) {
        match rhs {
            0 => {}
            s if s >= umax::BITS => {
                self.hi = self.lo << (s - umax::BITS);
                self.lo = 0;
            }
            s => {
                self.hi <<= s;
                self.hi |= self.lo >> (umax::BITS - s);
                self.lo <<= s;
            }
        }
    }
}

impl Shr<u8> for udouble {
    type Output = Self;
    fn shr(self, rhs: u8) -> Self::Output {
        match rhs {
            0 => self,
            s if s >= umax::BITS as u8 => Self {
                lo: self.hi >> (rhs - umax::BITS as u8),
                hi: 0,
            },
            s => Self {
                hi: self.hi >> s,
                lo: (self.lo >> s) | (self.hi << (umax::BITS as u8 - s)),
            },
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
            },
        }
    }
}

impl ShrAssign<u8> for udouble {
    fn shr_assign(&mut self, rhs: u8) {
        match rhs {
            0 => {}
            s if s >= umax::BITS as u8 => {
                self.lo = self.hi >> (rhs - umax::BITS as u8);
                self.hi = 0;
            }
            s => {
                self.lo >>= s;
                self.lo |= self.hi << (umax::BITS as u8 - s);
                self.hi >>= s;
            }
        }
    }
}

impl ShrAssign<u32> for udouble {
    fn shr_assign(&mut self, rhs: u32) {
        match rhs {
            0 => {}
            s if s >= umax::BITS => {
                self.lo = self.hi >> (rhs - umax::BITS);
                self.hi = 0;
            }
            s => {
                self.lo >>= s;
                self.lo |= self.hi << (umax::BITS - s);
                self.hi >>= s;
            }
        }
    }
}

impl BitAnd for udouble {
    type Output = Self;
    fn bitand(self, rhs: Self) -> Self::Output {
        Self { lo: self.lo & rhs.lo, hi: self.hi & rhs.hi }
    }
}

impl BitAndAssign for udouble {
    fn bitand_assign(&mut self, rhs: Self) {
        self.lo &= rhs.lo;
        self.hi &= rhs.hi;
    }
}

impl BitOr for udouble {
    type Output = Self;
    fn bitor(self, rhs: Self) -> Self::Output {
        Self { lo: self.lo | rhs.lo, hi: self.hi | rhs.hi }
    }
}

impl BitOrAssign for udouble {
    fn bitor_assign(&mut self, rhs: Self) {
        self.lo |= rhs.lo;
        self.hi |= rhs.hi;
    }
}

impl BitXor for udouble {
    type Output = Self;
    fn bitxor(self, rhs: Self) -> Self::Output {
        Self { lo: self.lo ^ rhs.lo, hi: self.hi ^ rhs.hi }
    }
}

impl BitXorAssign for udouble {
    fn bitxor_assign(&mut self, rhs: Self) {
        self.lo ^= rhs.lo;
        self.hi ^= rhs.hi;
    }
}

impl Not for udouble {
    type Output = Self;
    fn not(self) -> Self::Output {
        Self { lo: !self.lo, hi: !self.hi }
    }
}

impl udouble {
    pub const fn leading_zeros(self) -> u32 {
        if self.hi == 0 {
            self.lo.leading_zeros() + umax::BITS
        } else {
            self.hi.leading_zeros()
        }
    }

    // similar to udiv_qrnnd
    pub fn div_rem(self, other: Self) -> (Self, Self) {
        let mut n = self; // numerator
        let mut d = other; // denominator
        let mut q = Self::zero(); // quotient

        let nbits = (2 * umax::BITS - n.leading_zeros()) as u8; // assuming umax = u128
        let dbits = (2 * umax::BITS - d.leading_zeros()) as u8;
        assert!(dbits != 0, "division by zero");

        // Early return in case we are dividing by a larger number than us
        if nbits < dbits {
            return (q, n);
        }

        // Bitwise long division
        let mut shift = nbits - dbits;
        d <<= shift;
        loop {
            if n >= d {
                q += 1;
                n -= d;
            }
            if shift == 0 {
                break;
            }

            d >>= 1u8;
            q <<= 1u8;
            shift -= 1;
        }
        (q, n)
    }
}

impl Mul for udouble {
    type Output = Self;
    fn mul(self, rhs: Self) -> Self::Output {
        // TODO: use checked_mul instead of overflowing_mul
        let (m, overflow) = self.overflowing_mul(rhs);
        assert!(!overflow, "multiplication overflow!");
        m
    }
}

impl Div for udouble {
    type Output = Self;
    fn div(self, rhs: Self) -> Self::Output {
        self.div_rem(rhs).0
    }
}

impl Rem for udouble {
    type Output = Self;
    fn rem(self, rhs: Self) -> Self::Output {
        self.div_rem(rhs).1
    }
}

impl Mul<umax> for udouble {
    type Output = Self;
    fn mul(self, rhs: umax) -> Self::Output {
        let (m, overflow) = self.overflowing_muls(rhs);
        assert!(!overflow, "multiplication overflow!");
        m
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
        assert!(udouble { lo: 1, hi: 0 }.is_one());

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
        const TWO: udouble = udouble { hi: 0, lo: 2 };
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
        assert_eq!((MAX << 1u8) + 1, ONEMAX);
        assert_eq!(ONE << 200u8, udouble { lo: 0, hi: 1 << (200 - umax::BITS) });
        assert_eq!(ONEZERO >> umax::BITS, ONE);
        assert_eq!(ONEMAX >> 1u8, MAX);

        assert_eq!(ONE * MAX, MAX);
        assert_eq!(ONE * ONEMAX, ONEMAX);
        assert_eq!(TWO * ONEMAX, ONEMAX + ONEMAX);
        assert_eq!(MAX / ONE, MAX);
        assert_eq!(MAX / MAX, ONE);
        assert_eq!(ONE / MAX, udouble::zero());
        assert_eq!(ONEMAX / ONE, ONEMAX);
        assert_eq!(ONEMAX / ONEMAX, ONE);
        assert_eq!(ONEMAX / MAX, TWO);
        assert_eq!(ONEMAX / TWO, MAX);
        assert_eq!(ONE % MAX, ONE);
        assert_eq!(TWO % MAX, TWO);
        assert_eq!(ONEMAX % MAX, ONE);
        assert_eq!(ONEMAX % TWO, ONE);

        let (m, overflow) = ONEMAX.overflowing_mul(MAX);
        assert_eq!(m, udouble { lo: 1, hi: umax::MAX - 2});
        assert!(overflow);
        let (m, overflow) = TWOZERO.overflowing_mul(MAX);
        assert_eq!(m, udouble { lo: 0, hi: umax::MAX - 1});
        assert!(overflow);
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
