//! This crate provides efficient Modular arithmetic operations for various integer types,
//! including primitive integers and `num-bigint`. The latter option is enabled optionally.
//!
//! To achieve fast modular arithmetics, convert integers to any [ModularInteger] implementation
//! use static `new()` or associated [ModularInteger::new()]. [MontgomeryInt] is a builtin implementation
//! of [ModularInteger] based on the Montgomery form.
//!
//! Example code:
//! ```rust
//! use num_modular::{ModularCoreOps, ModularInteger, MontgomeryInt};
//!
//! // directly using methods in ModularCoreOps
//! let (x, y, m) = (12u8, 13u8, 5u8);
//! assert_eq!(x.mulm(y, &m), x * y % m);
//!
//! // convert integers into ModularInteger
//! let mx = MontgomeryInt::new(x, m);
//! let my = mx.new(y); // faster than static new()
//! assert_eq!((mx * my).residue(), x * y % m);
//! ```
//!

// XXX: consider implementing lookup table based modulo?
// REF: https://eprint.iacr.org/2014/040.pdf

#![no_std]
#[cfg(any(feature = "std", test))]
extern crate std;

use core::ops::{Add, Mul, Neg, Sub};

/// This trait describes core modular arithmetic operations
pub trait ModularCoreOps<Rhs = Self, Modulus = Self> {
    type Output;

    /// Return (self + rhs) % m
    fn addm(self, rhs: Rhs, m: Modulus) -> Self::Output;

    /// Return (self + rhs) % m
    fn subm(self, rhs: Rhs, m: Modulus) -> Self::Output;

    /// Return (self * rhs) % m
    fn mulm(self, rhs: Rhs, m: Modulus) -> Self::Output;

    /// Return (-self) % m and make sure the result is normalized in range [0,m)
    fn negm(self, m: Modulus) -> Self::Output;
}

// TODO (v0.3): checked_invm, checked_jacobi, checked_kronecker, checked_legendre
/// This trait describes modular arithmetic operations
pub trait ModularOps<Rhs = Self, Modulus = Self>: ModularCoreOps<Rhs, Modulus> {
    /// Return (self ^ exp) % m
    fn powm(self, exp: Rhs, m: Modulus) -> Self::Output;

    /// Calculate modular inverse (x such that self*x = 1 mod m).
    ///
    /// This operation is only available for integer that is coprime to `m`
    fn invm(self, m: Modulus) -> Option<Self::Output>
    where
        Self: Sized;

    /// Calculate Jacobi Symbol (a|n), where a is self
    ///
    /// # Panics
    /// if n is negative or even
    fn jacobi(self, n: Modulus) -> i8;

    /// Calculate Kronecker Symbol (a|n), where a is self
    fn kronecker(self, n: Modulus) -> i8;

    /// Calculate Legendre Symbol (a|n), where a is self.
    ///
    /// Note that this function doesn't perform primality check, since
    /// is costly. So if n is not a prime, the result is not reasonable.
    ///
    /// # Panics
    /// if n is not prime
    fn legendre(self, n: Modulus) -> i8;

    // TODO: Modular sqrt aka Quadratic residue, follow the behavior of FLINT `n_sqrtmod`
    // fn sqrtm(self, m: Modulus);

    // TODO: Discrete log aka index, follow the behavior of FLINT `n_discrete_log_bsgs`
    // fn logm(self, base: Modulus, m: Modulus);
}

/// Provides a utility function to convert signed integers into unsigned modular form
pub trait ModularAbs<Modulus> {
    /// Return |self| % m
    fn absm(self, m: &Modulus) -> Modulus;
}

/// Represents an number defined in a modulo ring ℤ/nℤ
///
/// The operators should panic if the modulus of two number
/// are not the same.
pub trait ModularInteger:
    Sized
    + PartialEq
    + Add<Self, Output = Self>
    + Sub<Self, Output = Self>
    + Neg<Output = Self>
    + Mul<Self, Output = Self>
{
    /// The underlying representation type of the integer
    type Base;

    /// Return the modulus of the ring
    fn modulus(&self) -> &Self::Base;

    /// Return the normalized residue of this integer in the ring
    fn residue(&self) -> Self::Base;

    /// Convert an normal integer into the same ring.
    ///
    /// This method should be perferred over the static
    /// constructor to prevent unnecessary overhead of pre-computation.
    fn new(&self, n: Self::Base) -> Self;
}

mod barret;
mod double;
mod monty;
mod prim;
mod mersenne;

pub use double::{udouble, umax};
#[cfg(std)]
pub use monty::MontgomeryBigint;
pub use monty::{Montgomery, MontgomeryInt};

#[cfg(feature = "num-bigint")]
mod bigint;

// tests for ModularOps goes here
#[cfg(test)]
mod tests;
