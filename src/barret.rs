/// Operations related to Barret reduction
pub trait Barret {}

// TODO: implement methods with Barret reduction (like) methods
// Version 1: Original barret reduction (for x mod n)
// - Choose k = ceil(log2(n))
// - Precompute r = floor(2^(k+1)/n)
// - t = x - floor(x*r/2^(k+1)) * n
// - if t > n, t -= n
// - return t
//
// Version 2: Fixed point barret reduction
// - Similar to version 1
// - Ref (u128): https://math.stackexchange.com/questions/3455277/barrett-reduction-possible-without-overflow-and-floating-point-arithmetic
// - Ref: https://stackoverflow.com/a/58470455
//
// Version 3: Floating point barret reduction
// - Using floating point to store r
// - Ref: http://flintlib.org/doc/ulong_extras.html#c.n_mulmod_precomp
//
// Version 4: "Improved division by invariant integers"
// - Ref: https://gmplib.org/~tege/division-paper.pdf
//