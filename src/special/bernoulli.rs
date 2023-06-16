//**********************************************************************
// This file is part of sci-rs                                         *
// Copyright 2023 Matthew R. Hennefarth                                *
//**********************************************************************

use crate::special::Comb;
use num_traits::{cast, Float, FromPrimitive, PrimInt};
use std::ops::AddAssign;

pub trait Bernoulli {
    fn bernoulli<Output>(self) -> Vec<Output>
    where
        Output: Float + AddAssign + FromPrimitive;
}

impl<Int> Bernoulli for Int
where
    Int: PrimInt + FromPrimitive,
{
    fn bernoulli<Output>(self) -> Vec<Output>
    where
        Output: Float + AddAssign + FromPrimitive,
    {
        bernoulli_rec::<Self, Output>(self)
    }
}

/// Returns first $n$ Bernoulli numbers using recurssion.
///
/// bernoulli_rec(3) -> [0.0, -0.5, 1/6, 0.0] (list of 4)
pub(crate) fn bernoulli_rec<T, Output>(n: T) -> Vec<Output>
where
    T: PrimInt + FromPrimitive,
    Output: Float + AddAssign + FromPrimitive,
{
    assert!(n >= T::zero());

    let result_len = (n + T::one()).to_usize().unwrap();
    const CACHE_SIZE: usize = 4;
    let mut result = vec![
        Output::one(),
        cast(-0.5).unwrap(),
        cast(1.0 / 6.0).unwrap(),
        Output::zero(),
    ];
    debug_assert_eq!(result.len(), CACHE_SIZE);

    if result_len <= CACHE_SIZE {
        return result[..result_len].to_vec();
    }

    result.reserve_exact(result_len - CACHE_SIZE);

    for i in CACHE_SIZE..result_len {
        if i % 2 == 1 {
            result.push(Output::zero());
        } else {
            let np1 = i + 1; // n+1
            let mut bn = Output::zero(); // B_n
            for (k, bk) in result.iter().enumerate() {
                bn += Output::from(np1.choose(k)).unwrap() * *bk;
            }
            bn = -bn / Output::from_usize(np1).unwrap();
            result.push(bn);
        }
    }
    result
}

#[cfg(test)]
mod tests {
    use super::*;

    const PRECISION: f64 = 1E-14;

    #[test]
    fn test_bernoulli() {
        for i in 0..4 {
            assert_eq!(i.bernoulli::<f32>().len(), i + 1);
        }

        assert_almost_eq!(
            *4.bernoulli::<f64>().last().unwrap(),
            -1.0 / 30.0,
            PRECISION
        );
        assert_almost_eq!(*6.bernoulli::<f64>().last().unwrap(), 1.0 / 42.0, PRECISION);
        assert_almost_eq!(
            *8.bernoulli::<f64>().last().unwrap(),
            -1.0 / 30.0,
            PRECISION
        );
        assert_almost_eq!(
            *10.bernoulli::<f64>().last().unwrap(),
            5.0 / 66.0,
            PRECISION
        );

        for i in (5..=15).step_by(2) {
            assert_eq!(*i.bernoulli::<f64>().last().unwrap(), 0.0);
        }
    }

    #[test]
    fn test_bernoulli_rec() {
        for i in 0..4 {
            assert_eq!(bernoulli_rec::<usize, f32>(i).len(), i + 1);
        }

        assert_almost_eq!(
            *bernoulli_rec::<usize, f64>(4).last().unwrap(),
            -1.0 / 30.0,
            PRECISION
        );
        assert_almost_eq!(
            *bernoulli_rec::<usize, f64>(6).last().unwrap(),
            1.0 / 42.0,
            PRECISION
        );
        assert_almost_eq!(
            *bernoulli_rec::<usize, f64>(8).last().unwrap(),
            -1.0 / 30.0,
            PRECISION
        );
        assert_almost_eq!(
            *bernoulli_rec::<usize, f64>(10).last().unwrap(),
            5.0 / 66.0,
            PRECISION
        );

        for i in (5..=15).step_by(2) {
            assert_eq!(*bernoulli_rec::<usize, f64>(i).last().unwrap(), 0.0);
        }
    }
}
