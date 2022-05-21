use crate::{udouble, umax, ModularUnaryOps, Reducer};
use crate::reduced::impl_reduced_binary_pow;

// FIXME: use unchecked operators to speed up calculation (after https://github.com/rust-lang/rust/issues/85122)
/// An unsigned integer modulo (pseudo) Mersenne primes `2^P - K`, it supports `P` up to 127 and `K < 2^(P-1)`
///
/// IMPORTANT NOTE: this class assumes that `2^P-K` is a prime. During compliation, we don't do full check
/// of the primality of `2^P-K`. If it's not a prime, then the modular division and inverse will panic.
#[derive(Debug, Clone, Copy)]
pub struct Mersenne<const P: u8, const K: umax>();

// TODO(v0.5): implement Mersenne as a Reducer
// XXX: support other primes as modulo, such as solinas prime, proth prime

impl<const P: u8, const K: umax> Mersenne<P, K> {
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
}

impl<const P: u8, const K: umax> Reducer<umax> for Mersenne<P, K> {
    type Modulus = umax;

    #[inline]
    fn new(m: &Self::Modulus) -> Self {
        debug_assert!(m == &Self::MODULUS);
        assert!(P <= 127);
        assert!(K > 0 && K < 2u128.pow(P as u32 - 1) && K % 2 == 1);
        assert!(
            Self::MODULUS % 3 != 0
                && Self::MODULUS % 5 != 0
                && Self::MODULUS % 7 != 0
                && Self::MODULUS % 11 != 0
                && Self::MODULUS % 13 != 0
        ); // error on easy composites
        Self {}
    }
    #[inline]
    fn transform(target: umax, m: &Self::Modulus) -> umax {
        debug_assert!(m == &Self::MODULUS);
        Self::reduce_single(target)
    }
    #[inline]
    fn residue(&self, target: umax, m: &Self::Modulus) -> umax {
        debug_assert!(m == &Self::MODULUS);
        target
    }
    #[inline]
    fn modulus(m: &Self::Modulus) -> umax {
        debug_assert!(m == &Self::MODULUS);
        Self::MODULUS
    }
    #[inline]
    fn is_zero(&self, target: &umax, m: &Self::Modulus) -> bool {
        debug_assert!(m == &Self::MODULUS);
        target == &0
    }

    #[inline]
    fn add(&self, lhs: umax, rhs: umax, m: &Self::Modulus) -> umax {
        debug_assert!(m == &Self::MODULUS);
        let mut sum = lhs + rhs;
        if sum >= Self::MODULUS {
            sum -= Self::MODULUS
        }
        sum
    }
    #[inline]
    fn sub(&self, lhs: umax, rhs: umax, m: &Self::Modulus) -> umax {
        debug_assert!(m == &Self::MODULUS);
        if lhs >= rhs {
            lhs - rhs
        } else {
            Self::MODULUS - (rhs - lhs)
        }
    }
    #[inline]
    fn double(&self, target: umax, m: &Self::Modulus) -> umax {
        debug_assert!(m == &Self::MODULUS);
        self.add(target, target, &Self::MODULUS)
    }
    #[inline]
    fn neg(&self, target: umax, m: &Self::Modulus) -> umax {
        debug_assert!(m == &Self::MODULUS);
        if target == 0 {
            0
        } else {
            Self::MODULUS - target
        }
    }
    #[inline]
    fn mul(&self, lhs: umax, rhs: umax, m: &Self::Modulus) -> umax {
        debug_assert!(m == &Self::MODULUS);
        if (P as u32) < (umax::BITS / 2) {
            Self::reduce_single(lhs * rhs)
        } else {
            Self::reduce_double(udouble::widening_mul(lhs, rhs))
        }
    }
    #[inline]
    fn inv(&self, target: umax, m: &Self::Modulus) -> Option<umax> {
        debug_assert!(m == &Self::MODULUS);
        if (P as u32) < usize::BITS {
            (target as usize)
                .invm(&(Self::MODULUS as usize))
                .map(|v| v as umax)
        } else {
            target.invm(&Self::MODULUS)
        }
    }
    #[inline]
    fn square(&self, target: umax, m: &Self::Modulus) -> umax {
        debug_assert!(m == &Self::MODULUS);
        if (P as u32) < (umax::BITS / 2) {
            Self::reduce_single(target * target)
        } else {
            Self::reduce_double(udouble::widening_square(target))
        }
    }

    impl_reduced_binary_pow!(u128);
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{ModularCoreOps, ModularPow};
    use rand::random;

    type M1 = Mersenne::<31, 1>;
    const P1: u128 = (1 << 31) - 1;
    type M2 = Mersenne::<61, 1>;
    const P2: u128 = (1 << 61) - 1;
    type M3 = Mersenne::<127, 1>;
    const P3: u128 = (1 << 127) - 1;
    type M4 = Mersenne::<32, 5>;
    const P4: u128 = (1 << 32) - 5;
    type M5 = Mersenne::<56, 5>;
    const P5: u128 = (1 << 56) - 5;
    type M6 = Mersenne::<122, 3>;
    const P6: u128 = (1 << 122) - 3;

    const NRANDOM: u32 = 10;

    #[test]
    fn creation_test() {
        // random creation test
        for _ in 0..NRANDOM {
            let a = random::<u128>();

            assert_eq!(M1::new(&P1).residue(M1::transform(a, &P1), &P1), a % P1);
            assert_eq!(M2::new(&P2).residue(M2::transform(a, &P2), &P2), a % P2);
            assert_eq!(M3::new(&P3).residue(M3::transform(a, &P3), &P3), a % P3);
            assert_eq!(M4::new(&P4).residue(M4::transform(a, &P4), &P4), a % P4);
            assert_eq!(M5::new(&P5).residue(M5::transform(a, &P5), &P5), a % P5);
            assert_eq!(M6::new(&P6).residue(M6::transform(a, &P6), &P6), a % P6);
        }
    }

    #[test]
    fn test_against_prim() {
        for _ in 0..NRANDOM {
            let (a, b) = (random::<u128>(), random::<u128>());
            let e = random::<u8>();

            // mod 2^31-1
            let am = M1::transform(a, &P1);
            let bm = M1::transform(b, &P1);
            let r = M1::new(&P1);
            assert_eq!(r.add(am, bm, &P1), a.addm(b, &P1));
            assert_eq!(r.sub(am, bm, &P1), a.subm(b, &P1));
            assert_eq!(r.mul(am, bm, &P1), a.mulm(b, &P1));
            assert_eq!(r.neg(am, &P1), a.negm(&P1));
            assert_eq!(r.inv(am, &P1), a.invm(&P1));
            assert_eq!(r.double(am, &P1), a.dblm(&P1));
            assert_eq!(r.square(am, &P1), a.sqm(&P1));
            assert_eq!(r.pow(am, e as u128, &P1), a.powm(e as u128, &P1));

            // mod 2^61-1
            let am = M2::transform(a, &P2);
            let bm = M2::transform(b, &P2);
            let r = M2::new(&P2);
            assert_eq!(r.add(am, bm, &P2), a.addm(b, &P2));
            assert_eq!(r.sub(am, bm, &P2), a.subm(b, &P2));
            assert_eq!(r.mul(am, bm, &P2), a.mulm(b, &P2));
            assert_eq!(r.neg(am, &P2), a.negm(&P2));
            assert_eq!(r.inv(am, &P2), a.invm(&P2));
            assert_eq!(r.double(am, &P2), a.dblm(&P2));
            assert_eq!(r.square(am, &P2), a.sqm(&P2));
            assert_eq!(r.pow(am, e as u128, &P2), a.powm(e as u128, &P2));

            // mod 2^127-1
            let am = M3::transform(a, &P3);
            let bm = M3::transform(b, &P3);
            let r = M3::new(&P3);
            assert_eq!(r.add(am, bm, &P3), a.addm(b, &P3));
            assert_eq!(r.sub(am, bm, &P3), a.subm(b, &P3));
            assert_eq!(r.mul(am, bm, &P3), a.mulm(b, &P3));
            assert_eq!(r.neg(am, &P3), a.negm(&P3));
            assert_eq!(r.inv(am, &P3), a.invm(&P3));
            assert_eq!(r.double(am, &P3), a.dblm(&P3));
            assert_eq!(r.square(am, &P3), a.sqm(&P3));
            assert_eq!(r.pow(am, e as u128, &P3), a.powm(e as u128, &P3));

            // mod 2^32-5
            let am = M4::transform(a, &P4);
            let bm = M4::transform(b, &P4);
            let r = M4::new(&P4);
            assert_eq!(r.add(am, bm, &P4), a.addm(b, &P4));
            assert_eq!(r.sub(am, bm, &P4), a.subm(b, &P4));
            assert_eq!(r.mul(am, bm, &P4), a.mulm(b, &P4));
            assert_eq!(r.neg(am, &P4), a.negm(&P4));
            assert_eq!(r.inv(am, &P4), a.invm(&P4));
            assert_eq!(r.double(am, &P4), a.dblm(&P4));
            assert_eq!(r.square(am, &P4), a.sqm(&P4));
            assert_eq!(r.pow(am, e as u128, &P4), a.powm(e as u128, &P4));

            // mod 2^56-5
            let am = M5::transform(a, &P5);
            let bm = M5::transform(b, &P5);
            let r = M5::new(&P5);
            assert_eq!(r.add(am, bm, &P5), a.addm(b, &P5));
            assert_eq!(r.sub(am, bm, &P5), a.subm(b, &P5));
            assert_eq!(r.mul(am, bm, &P5), a.mulm(b, &P5));
            assert_eq!(r.neg(am, &P5), a.negm(&P5));
            assert_eq!(r.inv(am, &P5), a.invm(&P5));
            assert_eq!(r.double(am, &P5), a.dblm(&P5));
            assert_eq!(r.square(am, &P5), a.sqm(&P5));
            assert_eq!(r.pow(am, e as u128, &P5), a.powm(e as u128, &P5));

            // mod 2^122-3
            let am = M6::transform(a, &P6);
            let bm = M6::transform(b, &P6);
            let r = M6::new(&P6);
            assert_eq!(r.add(am, bm, &P6), a.addm(b, &P6));
            assert_eq!(r.sub(am, bm, &P6), a.subm(b, &P6));
            assert_eq!(r.mul(am, bm, &P6), a.mulm(b, &P6));
            assert_eq!(r.neg(am, &P6), a.negm(&P6));
            assert_eq!(r.inv(am, &P6), a.invm(&P6));
            assert_eq!(r.double(am, &P6), a.dblm(&P6));
            assert_eq!(r.square(am, &P6), a.sqm(&P6));
            assert_eq!(r.pow(am, e as u128, &P6), a.powm(e as u128, &P6));
        }
    }
}
