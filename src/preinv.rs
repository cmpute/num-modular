
// implement modular inverse based quick divisibility check here (as used in num-prime)
// it differs from barret reduction that, barret finds a rational approximation of 1/d, while this method finds d^-1 mod B

use crate::ModularUnaryOps;

/// See https://math.stackexchange.com/a/1251328 for the explanation of the trick
#[derive(Debug, Clone, Copy)]
pub struct PreInv<T> {
    d_inv: T, // modular inverse of divisor
    q_lim: T, // limit of residue
}

impl PreInv<u8> {
    /// Construct the preinv instance with raw values.
    /// 
    /// This function can be used to initialize preinv in a constant context, the divisor d
    /// is required only for verification of d_inv and q_lim.
    #[inline]
    pub const fn new(d: u8, d_inv: u8, q_lim: u8) -> Self {
        debug_assert!(d % 2 != 0, "only odd divisors are supported");
        debug_assert!(d.wrapping_mul(d_inv) == 1);
        debug_assert!(q_lim * d > (u8::MAX - d));

        Self { d_inv, q_lim }
    }
}

impl From<u8> for PreInv<u8> {
    #[inline]
    fn from(v: u8) -> Self {
        debug_assert!(v % 2 != 0, "only odd divisors are supported");
        let d_inv = (v as u16).invm(&(1u16 << u8::BITS)).unwrap() as u8;
        let q_lim = u8::MAX / v;
        Self { d_inv, q_lim }
    }
}

// TODO(v0.4.2): move to crate root
trait DivExact<Rhs, Precompute> : Sized {
    type Output;
    fn div_exact(self, d: Rhs, pre: Precompute) -> Option<Self::Output>;
}

impl DivExact<u8, PreInv<u8>> for u8 {
    type Output = u8;
    #[inline]
    fn div_exact(self, _: u8, pre: PreInv<u8>) -> Option<Self> {
        let q = self.wrapping_mul(pre.d_inv);
        if q <= pre.q_lim {
            Some(q)
        } else {
            None
        }
    }
}

impl DivExact<u8, PreInv<u8>> for u16 {
    type Output = u16;
    #[inline]
    fn div_exact(self, d: u8, pre: PreInv<u8>) -> Option<u16> {
        // this code comes from GNU factor
        // TODO: someone explain this? can this be extended to multi-precision integers?
        //       https://math.stackexchange.com/q/4436380/815652

        let (n1, n0) = ((self >> u8::BITS) as u8, self as u8);
        let q0 = n0.wrapping_mul(pre.d_inv);
        let nr0 = (q0 as u16) * (d as u16);
        let nr0 = (nr0 >> u8::BITS) as u8;
        if nr0 > n1 {
            return None
        }
        let nr1 = n1 - nr0;
        let q1 = nr1.wrapping_mul(pre.d_inv);
        if q1 > pre.q_lim {
            return None
        }
        Some(((q1 as u16) << u8::BITS) + q0 as u16)
    }
}

// XXX: we could implement the div_exact for big integers, in a similar way to support double width

#[cfg(test)]
mod tests {
    use rand::random;
    use super::*;

    #[test]
    fn div_exact_test() {
        const N: u8 = 100;
        for _ in 0..N {
            // u8 test
            let d = random::<u8>() | 1;
            let pre: PreInv<u8> = d.into();
    
            let n: u8 = random();
            let expect = if n % d == 0 { Some(n/d) } else { None };
            assert_eq!(n.div_exact(d, pre), expect, "{} / {}", n, d);
            let n: u16 = random();
            let expect = if n % (d as u16) == 0 { Some(n/(d as u16)) } else { None };
            assert_eq!(n.div_exact(d, pre), expect, "{} / {}", n, d);

            // u16 test
        }
    }
}
