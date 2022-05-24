//! All methods that using pre-computed inverse of the modulus will be contained in this module,
//! as it shares the idea of barret reduction.
//!
//! Version 1: Vanilla barret reduction (for x mod n, x < n^2)
//! - Choose k = ceil(log2(n))
//! - Precompute r = floor(2^(k+1)/n)
//! - t = x - floor(x*r/2^(k+1)) * n
//! - if t > n, t -= n
//! - return t
//!
//! Version 2: Full width barret reduction
//! - Similar to version 1 but support n up to full width
//! - Ref (u128): <https://math.stackexchange.com/a/3455956/815652>
//!
//! Version 3: Floating point barret reduction
//! - Using floating point to store r
//! - Ref: <http://flintlib.org/doc/ulong_extras.html#c.n_mulmod_precomp>
//!
//! Version 4: "Improved division by invariant integers" by Granlund
//! - Ref: <https://gmplib.org/~tege/division-paper.pdf>
//!        <https://gmplib.org/~tege/divcnst-pldi94.pdf>
//!
//! Comparison between vanilla Barret reduction and Montgomery reduction:
//! - Barret reduction requires one 2k-by-k bits and one k-by-k bits multiplication while Montgomery only involves two k-by-k multiplications
//! - Extra conversion step is required for Montgomery form to get a normal integer
//! (Referece: <https://www.nayuki.io/page/barrett-reduction-algorithm>)
//!
//! The latter two versions are efficient and practical for use.
// TODO: implement version 3 and version 4
