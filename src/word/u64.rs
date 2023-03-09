type Word = u64;
type DoubleWord = u128;

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
    word as DoubleWord
}

#[inline]
pub const fn split(dw: DoubleWord) -> (Word, Word) {
    (dw as Word, (dw >> Word::BITS) as Word)
}

#[inline]
pub const fn merge(low: Word, high: Word) -> DoubleWord {
    extend(low) | extend(high) << Word::BITS
}

#[inline]
pub const fn wmul(a: Word, b: Word) -> DoubleWord {
    extend(a) * extend(b)
}
