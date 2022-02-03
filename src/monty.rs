use crate::ModularOps;
use num_integer::Integer;
use num_traits::WrappingNeg;
use std::borrow::Borrow;
use std::convert::TryInto;
use std::rc::Rc;

/// Operations of a integer represented in Montgomery form. This data type can
/// be used in place of a normal integer with regard to modular arithmetics.
///
/// The generic type T represents the underlying integer representation, and
/// R=2^B will be used as the auxiliary modulus, where B is automatically selected
/// based on the size of T.
pub trait Montgomery: Sized {
    /// The type for inversion of the modulus.
    ///
    /// This type is usually the same as Self, but it can be smaller when using
    /// Montgomery form on multi-precision integer representations.
    type Inv;

    /// The type of integer with doubled width
    type Double;

    /// Calculate -(m^-1) mod R
    fn neginv(m: &Self) -> Self::Inv;

    /// Transform a normal integer into Montgomery form (compute `target*R mod m`)
    fn transform(target: Self, m: &Self) -> Self;

    /// Transform a montgomery form back to normal integer (compute `monty/R mod m`)
    fn reduce(monty: Self::Double, m: &Self, minv: &Self::Inv) -> Self;

    /// Calculate (lhs + rhs) mod m in Montgomery form
    fn mul(lhs: Self, rhs: Self, m: &Self, minv: &Self::Inv) -> Self;

    /// Calculate base ^ exp mod m in Montgomery form
    fn pow(base: Self, exp: Self, m: &Self, minv: &Self::Inv) -> Self;
}

// TODO: implement Montgomery for u32, u64, biguint
// REF: https://github.com/uutils/coreutils/blob/main/src/uu/factor/src/numeric/montgomery.rs#L68
//      https://crates.io/crates/modulo-n-tools
//      https://docs.rs/ibig/latest/ibig/modular/index.html
//      https://docs.rs/ring-algorithm/latest/ring_algorithm/
//      https://github.com/vks/discrete-log/blob/master/src/main.rs

/* Entry i contains (2i+1)^(-1) mod 2^8.  */
const BINVERT_TABLE: [u8; 128] = [
  0x01, 0xAB, 0xCD, 0xB7, 0x39, 0xA3, 0xC5, 0xEF,
  0xF1, 0x1B, 0x3D, 0xA7, 0x29, 0x13, 0x35, 0xDF,
  0xE1, 0x8B, 0xAD, 0x97, 0x19, 0x83, 0xA5, 0xCF,
  0xD1, 0xFB, 0x1D, 0x87, 0x09, 0xF3, 0x15, 0xBF,
  0xC1, 0x6B, 0x8D, 0x77, 0xF9, 0x63, 0x85, 0xAF,
  0xB1, 0xDB, 0xFD, 0x67, 0xE9, 0xD3, 0xF5, 0x9F,
  0xA1, 0x4B, 0x6D, 0x57, 0xD9, 0x43, 0x65, 0x8F,
  0x91, 0xBB, 0xDD, 0x47, 0xC9, 0xB3, 0xD5, 0x7F,
  0x81, 0x2B, 0x4D, 0x37, 0xB9, 0x23, 0x45, 0x6F,
  0x71, 0x9B, 0xBD, 0x27, 0xA9, 0x93, 0xB5, 0x5F,
  0x61, 0x0B, 0x2D, 0x17, 0x99, 0x03, 0x25, 0x4F,
  0x51, 0x7B, 0x9D, 0x07, 0x89, 0x73, 0x95, 0x3F,
  0x41, 0xEB, 0x0D, 0xF7, 0x79, 0xE3, 0x05, 0x2F,
  0x31, 0x5B, 0x7D, 0xE7, 0x69, 0x53, 0x75, 0x1F,
  0x21, 0xCB, 0xED, 0xD7, 0x59, 0xC3, 0xE5, 0x0F,
  0x11, 0x3B, 0x5D, 0xC7, 0x49, 0x33, 0x55, 0xFF
];

impl Montgomery for u8 {
    type Inv = u8;
    type Double = u16;

    fn neginv(m: &Self) -> Self {
        BINVERT_TABLE[((m >> 1) & 0x7F) as usize].wrapping_neg()
    }

    fn transform(target: Self, m: &Self) -> Self {
        (((target as u16) << 8) % (*m as u16)) as _
    }

    fn reduce(monty: Self::Double, m: &Self, minv: &Self::Inv) -> Self {
        // REDC algorithm
        debug_assert!(monty < ((*m as u16) << u8::BITS));

        let tm = (monty as u8).wrapping_mul(*minv);
        let (t, overflow) = monty.overflowing_add((tm as u16) * (*m as u16));
        let t = (t >> u8::BITS) as u8;
        
        // in case of overflow, we need to add another `R mod m` = `R - m`
        let t = if overflow {
            t + m.wrapping_neg()
        } else { t };

        if &t >= m { return t-m } else { return t }
    }

    fn mul(lhs: Self, rhs: Self, m: &Self, minv: &Self::Inv) -> Self {
        Montgomery::reduce((lhs as u16) * (rhs as u16), m, minv)
    }

    fn pow(base: Self, exp: Self, m: &Self, minv: &Self::Inv) -> Self {
        match exp {
            1 => base,
            2 => Montgomery::mul(base, base, m, minv),
            _ => {
                let mut multi = base;
                let mut exp = exp;
                let mut result = 1;
                while exp > 0 {
                    if exp & 1 > 0 {
                        result = Montgomery::mul(result, multi, m, minv);
                    }
                    multi = Montgomery::mul(multi, multi, m, minv);
                    exp >>= 1;
                }
                result
            }
        }
    }
}

/// A integer represented in Montgomery form, which can be used in place of normal integers.
pub struct MontgomeryInt<T: Integer + Montgomery> {
    /// The Montgomery representation of the integer.
    a: T,

    /// The modulus and its negated modular inverse.
    /// 
    /// It's stored as a pointer to prevent frequent copying. It also allows
    /// quick checking of the equity of two moduli.
    minv: Rc<(T, T::Inv)>
}

impl<T: Integer + Montgomery> MontgomeryInt<T> where T::Double : From<T> {
    pub fn modulus(&self) -> &T {
        &Borrow::<(T, T::Inv)>::borrow(&self.minv).0
    }

    pub fn residue(self) -> T {
        let minv = Borrow::<(T, T::Inv)>::borrow(&self.minv);
        Montgomery::reduce(T::Double::from(self.a), &minv.0, &minv.1)
    }

    /// Create a new instance
    pub fn new(n: T, m: T) -> Self {
        let inv = Montgomery::neginv(&m);
        let a = Montgomery::transform(n, &m);
        MontgomeryInt { a, minv: Rc::new((m, inv)) }
    }
}

// TODO: implement ModularInteger for MontgomeryInt
