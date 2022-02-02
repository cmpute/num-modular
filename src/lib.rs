//! This crate provides efficient Modular arithmetic operations
//! for various integer types, including primitive integers and
//! `num-bigint`. The latter option is enabled optionally.

/// This trait describes modular arithmetic operations
pub trait ModularOps<Rhs = Self, Modulus = Self> {
    type Output;

    /// Return (self + rhs) % m
    fn addm(self, rhs: Rhs, m: Modulus) -> Self::Output;

    /// Return (self + rhs) % m
    fn subm(self, rhs: Rhs, m: Modulus) -> Self::Output;

    /// Return (self * rhs) % m
    fn mulm(self, rhs: Rhs, m: Modulus) -> Self::Output;

    /// Return (self ^ exp) % m
    fn powm(self, exp: Rhs, m: Modulus) -> Self::Output;

    /// Return (-self) % m and make sure the result is normalized in range [0,m)
    fn negm(self, m: Modulus) -> Self::Output;

    /// Calculate modular inverse (x such that self*x = 1 mod m).
    fn invm(self, m: Modulus) -> Option<Self::Output>
    where
        Self: Sized;

    /// Calculate Jacobi Symbol (a|n), where a is self
    fn jacobi(self, n: Modulus) -> i8;

    // TODO: Calculate Kronecker Symbol (a|n), where a is self
    // fn kronecker(self, n: Modulus) -> i8;

    // TODO: ModularOps sqrt aka Quadratic residue
    // fn sqrtm(self, m: Modulus);
}

mod bigint;
mod prim;
mod monty;

// tests for ModularOps goes here
#[cfg(test)]
mod tests {
    use super::*;
    use rand;

    #[test]
    fn u64_basic_mod_test() {
        let a = rand::random::<u64>() % 100000;
        let m = rand::random::<u64>() % 100000;
        assert_eq!(a.addm(a, &m), (a + a) % m);
        assert_eq!(a.mulm(a, &m), (a * a) % m);
        assert_eq!(a.powm(3, &m), a.pow(3) % m);
    }

    #[test]
    #[cfg(feature="num-bigint")]
    fn biguint_basic_mod_test() {
        let mut rng = rand::thread_rng();
        let a = rand::random::<u128>();
        let ra = &BigUint::from(a);
        let m = rand::random::<u128>();
        let rm = &BigUint::from(m);
        assert_eq!(ra.addm(ra, rm), (ra + ra) % rm);
        assert_eq!(ra.mulm(ra, rm), (ra * ra) % rm);
        assert_eq!(ra.powm(BigUint::from(3u8), rm), ra.pow(3) % rm);
    }

    #[test]
    fn addm_test() {
        let m = 5;

        let test_cases: [(u8, u8, u8); 10] = [
            // [x, y, rem]: x + y = rem (mod m)
            (0, 0, 0),
            (1, 2, 3),
            (2, 1, 3),
            (2, 2, 4),
            (3, 2, 0),
            (2, 3, 0),
            (6, 1, 2),
            (1, 6, 2),
            (11, 7, 3),
            (7, 11, 3),
        ];

        for (x, y, r) in test_cases.iter() {
            assert_eq!(x.addm(y, &m), *r, "u8 x: {}, y: {}", x, y);
            assert_eq!(
                (*x as u16).addm(*y as u16, &(m as u16)),
                *r as u16,
                "u16 x: {}, y: {}",
                x,
                y
            );
            assert_eq!(
                (*x as u32).addm(*y as u32, &(m as u32)),
                *r as u32,
                "u32 x: {}, y: {}",
                x,
                y
            );
            assert_eq!(
                (*x as u64).addm(*y as u64, &(m as u64)),
                *r as u64,
                "u64 x: {}, y: {}",
                x,
                y
            );
            assert_eq!(
                (*x as u128).addm(*y as u128, &(m as u128)),
                *r as u128,
                "u128 x: {}, y: {}",
                x,
                y
            );

            #[cfg(feature="num-bigint")]
            {
                assert_eq!(
                    BigUint::from(*x).addm(BigUint::from(*y), &BigUint::from(m)),
                    BigUint::from(*r),
                    "biguint x: {}, y: {}",
                    x,
                    y
                );
            }
        }
    }

    #[test]
    fn subm_test() {
        let m = 7;

        let test_cases: [(u8, u8, u8); 10] = [
            // [x, y, rem]: x - y = rem (mod modu)
            (0, 0, 0),
            (11, 9, 2),
            (5, 2, 3),
            (2, 5, 4),
            (6, 7, 6),
            (1, 7, 1),
            (7, 1, 6),
            (0, 6, 1),
            (15, 1, 0),
            (1, 15, 0),
        ];

        for (x, y, r) in test_cases.iter() {
            assert_eq!(x.subm(y, &m), *r, "u8 x: {}, y: {}", x, y);
            assert_eq!(
                (*x as u16).subm(*y as u16, &(m as u16)),
                *r as u16,
                "u16 x: {}, y: {}",
                x,
                y
            );
            assert_eq!(
                (*x as u32).subm(*y as u32, &(m as u32)),
                *r as u32,
                "u32 x: {}, y: {}",
                x,
                y
            );
            assert_eq!(
                (*x as u64).subm(*y as u64, &(m as u64)),
                *r as u64,
                "u64 x: {}, y: {}",
                x,
                y
            );
            assert_eq!(
                (*x as u128).subm(*y as u128, &(m as u128)),
                *r as u128,
                "u128 x: {}, y: {}",
                x,
                y
            );
            
            #[cfg(feature="num-bigint")]
            {
                assert_eq!(
                    BigUint::from(*x).subm(BigUint::from(*y), &BigUint::from(m)),
                    BigUint::from(*r),
                    "biguint x: {}, y: {}",
                    x,
                    y
                );
            }
        }
    }

    #[test]
    fn invm_test() {
        let test_cases: [(u64, u64, u64); 8] = [
            // [a, m, x] s.t. a*x = 1 (mod m) is satisfied
            (5, 11, 9),
            (8, 11, 7),
            (10, 11, 10),
            (3, 5000, 1667),
            (1667, 5000, 3),
            (999, 5000, 3999),
            (999, 9_223_372_036_854_775_807, 3_619_181_019_466_538_655),
            (
                9_223_372_036_854_775_804,
                9_223_372_036_854_775_807,
                3_074_457_345_618_258_602,
            ),
        ];

        for (a, m, x) in test_cases.iter() {
            assert_eq!(
                ModularOps::<&u64>::invm(a, m).unwrap(),
                *x,
                "a: {}, m: {}",
                a,
                m
            );
            
            #[cfg(feature="num-bigint")]
            {
                assert_eq!(
                    ModularOps::<&BigUint>::invm(&BigUint::from(*a), &BigUint::from(*m)).unwrap(),
                    BigUint::from(*x),
                    "a: {}, m: {}",
                    a,
                    m
                );
            }
        }
    }

    #[test]
    fn jacobi_test() {
        let test_cases: [(u8, u8, i8); 15] = [
            (1, 1, 1),
            (15, 1, 1),
            (2, 3, -1),
            (29, 9, 1),
            (4, 11, 1),
            (17, 11, -1),
            (19, 29, -1),
            (10, 33, -1),
            (11, 33, 0),
            (12, 33, 0),
            (14, 33, -1),
            (15, 33, 0),
            (15, 37, -1),
            (29, 59, 1),
            (30, 59, -1),
        ];

        for (a, n, res) in test_cases.iter() {
            assert_eq!(ModularOps::<&u8>::jacobi(a, n), *res, "u8 a: {}, n: {}", a, n);
            assert_eq!(
                ModularOps::<&u16>::jacobi(&(*a as u16), &(*n as u16)),
                *res,
                "u16 a: {}, n: {}",
                a,
                n
            );
            assert_eq!(
                ModularOps::<&u32>::jacobi(&(*a as u32), &(*n as u32)),
                *res,
                "u32 a: {}, n: {}",
                a,
                n
            );
            assert_eq!(
                ModularOps::<&u64>::jacobi(&(*a as u64), &(*n as u64)),
                *res,
                "u64 a: {}, n: {}",
                a,
                n
            );
            assert_eq!(
                ModularOps::<&u128>::jacobi(&(*a as u128), &(*n as u128)),
                *res,
                "u128 a: {}, n: {}",
                a,
                n
            );

            #[cfg(feature="num-bigint")]
            {
                assert_eq!(
                    ModularOps::<&BigUint>::jacobi(&(BigUint::from(*a)), &(BigUint::from(*n))),
                    *res,
                    "biguint a: {}, n: {}",
                    a,
                    n
                );
            }
        }
    }
}
