use crate::ModularOps;
use num_integer::Integer;
use std::borrow::Borrow;
use std::rc::Rc;

/// Operations of a integer represented in Montgomery form. This data type can
/// be used in place of a normal integer with regard to modular arithmetics.
///
/// The generic type T represents the underlying integer representation, and
/// R=2^B will be used as the auxiliary modulus, where B is automatically selected
/// based on the size of T.
pub trait Montgomery: ModularOps + Sized {
    /// The type for inversion of the modulus.
    ///
    /// This type is usually the same as Self, but it can be smaller when using
    /// Montgomery form on multi-precision integer representations.
    type Inv;

    /// Calculate m^-1 mod R
    fn inv(m: &Self) -> Self::Inv;

    /// Transform a normal integer into Montgomery form (monty = target*R mod m)
    fn transform(target: Self, m: &Self, minv: &Self::Inv) -> Self;

    /// Transform a montgomery form back to normal integer. (target = monty/R mod m)
    fn reduce(monty: Self, m: &Self, minv: &Self::Inv) -> Self;

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

/// A integer represented in Montgomery form, which can be used in place of normal integers.
pub struct MontgomeryInt<T: Integer + Montgomery> {
    /// The Montgomery representation of the integer if the modulus is present. Otherwise
    /// the original number will be stored in a.
    a: T,

    /// The modulus. It's stored as a pointer to prevent frequent copying
    m: Option<Rc<T>>,

    /// The modular inverse of the modulus. It's calculated only when necessary
    minv: Option<T::Inv>,
}

impl<T: Integer + Montgomery> MontgomeryInt<T> {
    /// Get the modulus value
    pub fn modulus(&self) -> Option<&T> {
        self.m.as_ref().map(Borrow::borrow)
    }

    /// Recover the original integer
    pub fn reduce(self) -> T {
        match self.m {
            Some(rc) => {
                let inv = match self.minv {
                    Some(v) => v,
                    None => <T as Montgomery>::inv(rc.borrow()),
                };
                Montgomery::reduce(self.a, rc.borrow(), &inv)
            }
            None => self.a,
        }
    }
}

// TODO: implement trait Integer and ModularOps for MontgomeryInt
