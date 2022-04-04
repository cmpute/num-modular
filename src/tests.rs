use super::*;
use num_traits::Pow;
use rand;

#[cfg(feature = "num-bigint")]
use num_bigint::BigUint;

#[test]
fn u64_basic_mod_test() {
    let a = rand::random::<u64>() % 100000;
    let m = rand::random::<u64>() % 100000;
    assert_eq!(a.addm(a, &m), (a + a) % m);
    assert_eq!(a.mulm(a, &m), (a * a) % m);
    assert_eq!(a.powm(3, &m), a.pow(3) % m);
}

#[test]
#[cfg(feature = "num-bigint")]
fn biguint_basic_mod_test() {
    let a = rand::random::<u128>();
    let ra = &BigUint::from(a);
    let m = rand::random::<u128>();
    let rm = &BigUint::from(m);
    assert_eq!(ra.addm(ra, rm), (ra + ra) % rm);
    assert_eq!(ra.mulm(ra, rm), (ra * ra) % rm);
    assert_eq!(ra.powm(BigUint::from(3u8), rm), ra.pow(3) % rm);
}

#[test]
fn legendre_test() {
    let test_cases: [(u8, u8, i8); 18] = [
        (0, 11, 0),
        (1, 11, 1),
        (2, 11, -1),
        (4, 11, 1),
        (7, 11, -1),
        (10, 11, -1),
        (0, 17, 0),
        (1, 17, 1),
        (2, 17, 1),
        (4, 17, 1),
        (9, 17, 1),
        (10, 17, -1),
        (0, 101, 0),
        (1, 101, 1),
        (2, 101, -1),
        (4, 101, 1),
        (9, 101, 1),
        (10, 101, -1),
    ];

    for (a, n, res) in test_cases.iter() {
        assert_eq!(ModularOps::<&u8>::legendre(a, n), *res);
        assert_eq!(
            ModularOps::<&u16>::legendre(&(*a as u16), &(*n as u16)),
            *res
        );
        assert_eq!(
            ModularOps::<&u32>::legendre(&(*a as u32), &(*n as u32)),
            *res
        );
        assert_eq!(
            ModularOps::<&u64>::legendre(&(*a as u64), &(*n as u64)),
            *res
        );
        assert_eq!(
            ModularOps::<&u128>::legendre(&(*a as u128), &(*n as u128)),
            *res
        );

        #[cfg(feature = "num-bigint")]
        {
            assert_eq!(
                ModularOps::<&BigUint>::legendre(&(BigUint::from(*a)), &(BigUint::from(*n))),
                *res
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
        assert_eq!(ModularOps::<&u8>::jacobi(a, n), *res);
        assert_eq!(ModularOps::<&u16>::jacobi(&(*a as u16), &(*n as u16)), *res);
        assert_eq!(ModularOps::<&u32>::jacobi(&(*a as u32), &(*n as u32)), *res);
        assert_eq!(ModularOps::<&u64>::jacobi(&(*a as u64), &(*n as u64)), *res);
        assert_eq!(
            ModularOps::<&u128>::jacobi(&(*a as u128), &(*n as u128)),
            *res
        );

        #[cfg(feature = "num-bigint")]
        {
            assert_eq!(
                ModularOps::<&BigUint>::jacobi(&(BigUint::from(*a)), &(BigUint::from(*n))),
                *res
            );
        }
    }
}

#[test]
fn kronecker_test() {
    let test_cases: [(u8, u8, i8); 18] = [
        (0, 15, 0),
        (1, 15, 1),
        (2, 15, 1),
        (4, 15, 1),
        (7, 15, -1),
        (10, 15, 0),
        (0, 14, 0),
        (1, 14, 1),
        (2, 14, 0),
        (4, 14, 0),
        (9, 14, 1),
        (10, 14, 0),
        (0, 11, 0),
        (1, 11, 1),
        (2, 11, -1),
        (4, 11, 1),
        (9, 11, 1),
        (10, 11, -1),
    ];

    for (a, n, res) in test_cases.iter() {
        assert_eq!(ModularOps::<&u8>::kronecker(a, n), *res);
        assert_eq!(
            ModularOps::<&u16>::kronecker(&(*a as u16), &(*n as u16)),
            *res
        );
        assert_eq!(
            ModularOps::<&u32>::kronecker(&(*a as u32), &(*n as u32)),
            *res
        );
        assert_eq!(
            ModularOps::<&u64>::kronecker(&(*a as u64), &(*n as u64)),
            *res
        );
        assert_eq!(
            ModularOps::<&u128>::kronecker(&(*a as u128), &(*n as u128)),
            *res
        );

        #[cfg(feature = "num-bigint")]
        {
            assert_eq!(
                ModularOps::<&BigUint>::kronecker(&(BigUint::from(*a)), &(BigUint::from(*n))),
                *res
            );
        }
    }
}
