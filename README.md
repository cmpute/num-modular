# num-modular

A generic implementation of modular arithmetics in Rust. It provide basic operators and an type to represent integers in a modulo ring (using Montgomery form). Specifically the following operations are supported:

- Common modular arithmetics: `Add`, `Sub`, `Mul`, `Div`, `Neg`, `Inv`, `Pow`
- Efficient modular arithmetics in Montgomery form, or using pseudo Mersenne primes
- Legendre, Jacobi and Kronecker symbols

It also support various integer type backends, including primitive integers and `num-bigint`. Note that this crate also supports `[no_std]`. To enable `std` related functionalities, enable the `std` feature of the crate.

<!-- TODO for v1: support ibig-rs -->
<!-- TODO for v1: const functions -->
