//**********************************************************************
// This file is part of sci-rs                                         *
// Copyright 2023 Matthew R. Hennefarth                                *
//**********************************************************************

use crate::special::ZigZag;
use num_traits::{Float, FromPrimitive, PrimInt};

/// The Bernoulli numbers.
pub trait Bernoulli: Sized {
    /// Computes the first $2n$ Bernoulli numbers, omitting the odd ones.
    /// $$
    /// \sum_{n\geq 0} B_n \frac{z^n}{n!} = \frac{z}{e^z - 1}
    /// $$
    /// The Bernoulli numbers ($B_n$) are rational numbers defined by the above generating function. The first few Bernoulli numbers are
    /// $$
    /// B_0 = 1, B_1 = \pm \frac{1}{2}, B_2 = \frac{1}{6}, B_3 = 0, B_4 = -\frac{1}{30}
    /// $$
    /// Note that for all odd $n > 1$, $B_n = 0$. The reason for the $\pm$ for $B_1$ comes from the two conventions in the literature. Because of these two facts, we never compute the odd Bernoulli numbers. That is given an input $n$, we compute $B_0,\ldots, B_{2n}$ and return as a vector. For more details on the Bernoulli numbers, see the [wiki] and [dlmf] pages.
    ///
    /// # Examples
    /// ```
    /// use sci_rs::special::Bernoulli;
    /// assert_eq!(0.bernoulli_b2n::<f32>(), vec![1.0]); // [B_0]
    /// assert_eq!(2.bernoulli_b2n::<f32>(), vec![1.0, 1.0/6.0, -1.0/30.0]); // [B_0, B_2, B_4]
    /// ```
    /// Of course, we can control the input and output size as
    /// ```
    /// use sci_rs::special::Bernoulli;
    /// assert_eq!(1_i8.bernoulli_b2n::<f64>(), vec![1.0, 1.0/6.0]); // [B_0, B_2]
    /// assert_eq!(2_u64.bernoulli_b2n::<f64>(), vec![1.0, 1.0/6.0, -1.0/30.0]); // [B_0, B_2, B_4]
    /// ```
    /// # Notes
    /// The implementation here is based on the TangentNumber algorithm described in section 6.2 of "Fast Computation of Bernoulli, Tangent and Secant Numbers" ([arxiv:1109.0286]). These are computed using only integer arithmetic and are numerically stable. It then uses the following relationship to build the Bernoulli numbers from the tangent numbers:
    /// $$
    /// T_n = (-1)^{n-1}2^{2n}\left(2^{2n} - 1\right)\frac{B_{2n}}{2n}
    /// $$
    /// This leads to highly accurate, floating-point Bernoulli numbers.
    ///
    /// [wiki]: https://en.wikipedia.org/wiki/Bernoulli_number
    /// [dlmf]: https://dlmf.nist.gov/24.7#i
    /// [arxiv:1109.0286]: https://doi.org/10.48550/arXiv.1108.0286
    fn bernoulli_b2n<Output>(self) -> Vec<Output>
    where
        Output: Float + FromPrimitive;
}

macro_rules! bernoulli_impl {
    ($($T: ty)*) => ($(
        impl Bernoulli for $T {
            #[inline(always)]
            fn bernoulli_b2n<Output>(self) -> Vec<Output>
            where
                Output: Float + FromPrimitive,
            {
                bernoulli_b2n_zag::<Self, Output>(self)
            }
        }
    )*)
}

bernoulli_impl! {u8 u16 u32 u64 u128 usize i8 i16 i32 i64 i128 isize}

/// Expected output is the first $2n$ Bernoulli numbers such that:
///
/// 0 -> [1.0]
/// 1 -> [1.0, 1/6]
/// 2 -> [1.0, 1/6, -1.0 / 30.0]
///
/// Uses algorithm from [here](https://arxiv.org/abs/1108.0286), specifically eq. 14. Essentially using the relationship that
/// $$
/// T_n = (-1)^{n-1}2^{2n}\left(2^{2n} - 1\right)\frac{B_{2n}}{2n}
/// $$
fn bernoulli_b2n_zag<T, Output>(n: T) -> Vec<Output>
where
    T: PrimInt + FromPrimitive + ZigZag,
    Output: Float + FromPrimitive,
{
    if n.is_zero() {
        return vec![Output::one()];
    }

    let zags = n.zag();

    let modified_zags = zags
        .iter()
        .enumerate()
        .map(|(n, &t)| {
            let j = 2 * (n + 1) as u32;
            let abs_bn = Output::from(t * T::from_u32(j).unwrap()).unwrap()
                / Output::from_usize(4_usize.pow(j) - 2_usize.pow(j)).unwrap();
            if n & 1 == 0 {
                abs_bn
            } else {
                -abs_bn
            }
        })
        .collect::<Vec<Output>>();

    [vec![Output::one()], modified_zags].concat()
    //return zags.iter().map(|&z| Output::from(z).unwrap()).collect();
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_bernoulli_zag() {
        const ABSOLUTE_KNOWN_VALUES: [f64; 6] = [
            1.0,
            1.0 / 6.0,
            -1.0 / 30.0,
            1.0 / 42.0,
            -1.0 / 30.0,
            5.0 / 66.0,
        ];
        for i in 0..ABSOLUTE_KNOWN_VALUES.len() {
            assert_eq!(
                bernoulli_b2n_zag::<_, f64>(i),
                ABSOLUTE_KNOWN_VALUES[..i + 1]
            );
        }
    }
}
