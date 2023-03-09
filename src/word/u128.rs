type Word = u128;
type DoubleWord = crate::double::udouble;

#[inline]
pub const fn ones(n: u32) -> Word {
    if n == 0 {
        0
    } else {
        Word::MAX >> (Word::BITS - n)
    }
}

#[inline]
pub const fn extend(word: Word) -> DoubleWord {
    DoubleWord { lo: word, hi: 0 }
}

#[inline]
pub const fn split(dw: DoubleWord) -> (Word, Word) {
    (dw.lo, dw.hi)
}

#[inline]
pub const fn merge(low: Word, high: Word) -> DoubleWord {
    DoubleWord { lo: low, hi: high }
}

#[inline]
pub const fn wmul(a: Word, b: Word) -> DoubleWord {
    DoubleWord::widening_mul(a, b)
}
