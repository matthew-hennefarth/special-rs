//**********************************************************************
// This file is part of sci-rs                                         *
// Copyright 2023 Matthew R. Hennefarth                                *
//**********************************************************************

use crate::special::Factorial;
use num_traits::{FromPrimitive, PrimInt};

/// Tangent and Secant numbers.
pub trait ZigZag: Sized {
    /// Returns the first $n$ Tangent numbers ($T_n$).
    /// $$
    /// T_n = \frac{2^{2n}\left(2^{2n}-1\right)|B_{2n}|}{2n}
    /// $$
    /// where $B_n$ is a Bernoulli number. Tangent numbers are also called 'zag' numbers. They are also related to combinatorics where they represent the number of alternative permutations on $n=1, 3, 5, 7, \ldots$ items (where permutations are equivalent if they are reverses of each other)[Wolframalpha]. For more details on Tangent number, see the [dlmf] page.
    ///
    /// # Notes:
    /// The implementation here is based on the TangentNumbers algorithm described in section 6.1 of "Fast Computation of Bernoulli, Tangent and Secant Numbers" ([arxiv:1109.0286]). As such, the numbers do not suffer any numerical instability since they are directly computed via integer arithmetic.
    ///
    /// [Wolframalpha]: https://mathworld.wolfram.com/TangentNumber.html
    /// [dlmf]: https://dlmf.nist.gov/search/search?q=tangent%20numbers
    /// [arxiv:1109.0286]: https://doi.org/10.48550/arXiv.1108.0286
    fn zag(self) -> Vec<Self>;
}

macro_rules! zigzag_primint_impl {
    ($($T: ty)*) => ($(
        impl ZigZag for $T {
            #[inline(always)]
            fn zag(self) -> Vec<Self> {
                zag_nums(self)
            }
        }
    )*)
}

zigzag_primint_impl! {u8 u16 u32 u64 u128 usize i8 i16 i32 i64 i128 isize}

fn zag_nums<T>(n: T) -> Vec<T>
where
    T: PrimInt + Factorial + FromPrimitive,
{
    let mut zags: Vec<T> = (0..n.to_usize().unwrap())
        .map(|k| T::factorial(T::from_usize(k).unwrap()))
        .collect::<Vec<T>>();

    for k in 1..zags.len() {
        for j in k..zags.len() {
            zags[j] = T::from_usize(j - k).unwrap() * zags[j - 1]
                + T::from_usize(j - k + 2).unwrap() * zags[j];
        }
    }

    return zags;
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_tangent_nums() {
        // Reference values taken from OEIS: A000182
        const REFERENCE: [u128; 17] = [
            1,
            2,
            16,
            272,
            7936,
            353792,
            22368256,
            1903757312,
            209865342976,
            29088885112832,
            4951498053124096,
            1015423886506852352,
            246921480190207983616,
            70251601603943959887872,
            23119184187809597841473536,
            8713962757125169296170811392,
            3729407703720529571097509625856,
        ];
        assert_eq!(17.zag(), REFERENCE);
    }
}
