//! This crate provides efficient Modular arithmetic operations for various integer types,
//! including primitive integers and `num-bigint`. The latter option is enabled optionally.
//!
//! To achieve fast modular arithmetics, convert integers to any [ModularInteger] implementation
//! use static `new()` or associated [ModularInteger::convert()] functions. Some builtin implementations
//! of [ModularInteger] includes [MontgomeryInt] and [MersenneInt].
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
//! let my = mx.convert(y); // faster than static MontgomeryInt::new(y, m)
//! assert_eq!((mx * my).residue(), x * y % m);
//! ```
//!

// XXX: Other fast modular arithmetic tricks
// REF: https://github.com/lemire/fastmod & https://arxiv.org/pdf/1902.01961.pdf
// REF: https://eprint.iacr.org/2014/040.pdf
// REF: https://github.com/ridiculousfish/libdivide/
// REF: Faster Interleaved Modular Multiplication Based on Barrett and Montgomery Reduction Methods (work for modulus in certain form)

#![no_std]
#[cfg(any(feature = "std", test))]
extern crate std;

use core::ops::{Add, Mul, Neg, Sub};

/// Core modular arithmetic operations.
///
/// Note that all functions will panic if the modulus is zero.
pub trait ModularCoreOps<Rhs = Self, Modulus = Self> {
    type Output;

    /// Return (self + rhs) % m
    fn addm(self, rhs: Rhs, m: Modulus) -> Self::Output;

    /// Return (self - rhs) % m
    fn subm(self, rhs: Rhs, m: Modulus) -> Self::Output;

    /// Return (self * rhs) % m
    fn mulm(self, rhs: Rhs, m: Modulus) -> Self::Output;
}

/// Core unary modular arithmetics
///
/// Note that all functions will panic if the modulus is zero.
pub trait ModularUnaryOps<Modulus = Self> {
    type Output;

    /// Return (-self) % m and make sure the result is normalized in range [0,m)
    fn negm(self, m: Modulus) -> Self::Output;

    /// Calculate modular inverse (x such that self*x = 1 mod m).
    ///
    /// This operation is only available for integer that is coprime to `m`. If not,
    /// the result will be [None].
    fn invm(self, m: Modulus) -> Option<Self::Output>;

    /// Calculate modular double ( x+x mod m)
    fn dblm(self, m: Modulus) -> Self::Output;

    /// Calculate modular square ( x*x mod m )
    fn sqm(self, m: Modulus) -> Self::Output;

    // TODO: Modular sqrt aka Quadratic residue, follow the behavior of FLINT `n_sqrtmod`
    // fn sqrtm(self, m: Modulus) -> Option<Self::Output>;
    // REF: https://stackoverflow.com/questions/6752374/cube-root-modulo-p-how-do-i-do-this
}

/// Modular power functions
pub trait ModularPow<Exp = Self, Modulus = Self> {
    type Output;

    /// Return (self ^ exp) % m
    fn powm(self, exp: Exp, m: Modulus) -> Self::Output;
}

/// Math symbols related to modular arithmetics
pub trait ModularSymbols<Modulus = Self> {
    /// Calculate Legendre Symbol (a|n), where a is self.
    ///
    /// Note that this function doesn't perform primality check, since
    /// is costly. So if n is not a prime, the result is not reasonable.
    ///
    /// # Panics
    /// if n is not prime
    #[inline]
    fn legendre(&self, n: Modulus) -> i8 {
        self.checked_legendre(n).expect("n shoud be a prime")
    }

    /// Checked version of [legendre()][ModularSymbols::legendre], return [None] if n is not prime
    fn checked_legendre(&self, n: Modulus) -> Option<i8>;

    /// Calculate Jacobi Symbol (a|n), where a is self
    ///
    /// # Panics
    /// if n is negative or even
    #[inline]
    fn jacobi(&self, n: Modulus) -> i8 {
        self.checked_jacobi(n)
            .expect("the Jacobi symbol is only defined for non-negative odd integers")
    }

    /// Checked version of [jacobi()][ModularSymbols::jacobi], return [None] if n is negative or even
    fn checked_jacobi(&self, n: Modulus) -> Option<i8>;

    /// Calculate Kronecker Symbol (a|n), where a is self
    fn kronecker(&self, n: Modulus) -> i8;
}

// TODO: Discrete log aka index, follow the behavior of FLINT `n_discrete_log_bsgs`
// REF: https://github.com/vks/discrete-log
// fn logm(self, base: Modulus, m: Modulus);

/// Collection of common modular arithmetic operations
pub trait ModularOps<Rhs = Self, Modulus = Self, Output = Self>:
    ModularCoreOps<Rhs, Modulus, Output = Output>
    + ModularUnaryOps<Modulus, Output = Output>
    + ModularPow<Rhs, Modulus, Output = Output>
    + ModularSymbols<Modulus>
{
}
impl<T, Rhs, Modulus> ModularOps<Rhs, Modulus> for T where
    T: ModularCoreOps<Rhs, Modulus, Output = T>
        + ModularUnaryOps<Modulus, Output = T>
        + ModularPow<Rhs, Modulus, Output = T>
        + ModularSymbols<Modulus>
{
}

/// Collection of operations similar to [ModularOps], but takes operands with references
pub trait ModularRefOps: for<'r> ModularOps<&'r Self, &'r Self> + Sized {}
impl<T> ModularRefOps for T where T: for<'r> ModularOps<&'r T, &'r T> {}

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
    fn convert(&self, n: Self::Base) -> Self;

    // Calculate the value of self + self
    fn double(self) -> Self;

    // Calculate the value of self * self
    fn square(self) -> Self;
}

// XXX: implement this trait for ff::PrimeField?
// TODO: implement invm_range (Modular inverse in certain range) and crt (Chinese Remainder Theorem), REF: bubblemath crate

mod barret;
mod double;
mod mersenne;
mod monty;
mod prim;

pub use double::{udouble, umax};
pub use mersenne::MersenneInt;
pub use monty::{Montgomery, MontgomeryInt};

#[cfg(feature = "num-bigint")]
mod bigint;
