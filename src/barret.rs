//! Comparison between Barret reduction and Montgomery reduction:
//! - Barret reduction requires one 2k-by-k bits and one k-by-k bits multiplication while Montgomery only involves two k-by-k multiplications
//! - Extra conversion step is required for Montgomery form to get a normal integer
//! (Referece: https://www.nayuki.io/page/barrett-reduction-algorithm)

/// Operations related to Barret reduction
pub trait Barret {}

// TODO: implement methods with Barret reduction (like) methods
// TODO(v0.5.x) Version 1: Original barret reduction (for x mod n)
// - Choose k = ceil(log2(n))
// - Precompute r = floor(2^(k+1)/n)
// - t = x - floor(x*r/2^(k+1)) * n
// - if t > n, t -= n
// - return t
//
// Version 2: Fixed point barret reduction
// - Similar to version 1
// - Ref (u128): https://math.stackexchange.com/questions/3455277/barrett-reduction-possible-without-overflow-and-floating-point-arithmetic (special case for 128bit)
// - Ref: https://stackoverflow.com/a/58470455
//
// Version 3: Floating point barret reduction (this method is very limited, don't consider)
// - Using floating point to store r
// - Ref: http://flintlib.org/doc/ulong_extras.html#c.n_mulmod_precomp
//
// Version 4: "Improved division by invariant integers"
// - Ref: https://gmplib.org/~tege/division-paper.pdf
//        https://gmplib.org/~tege/divcnst-pldi94.pdf
//