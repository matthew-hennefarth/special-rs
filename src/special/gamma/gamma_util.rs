//**********************************************************************
// This file is part of Sci-rs                                         *
//                                                                     *
// Sci-rs is licensed under the Apache License, Version 2.0 (the       *
// "License"); you may not use this file except in compliance with the *
// License. You may obtain a copy of the License at                    *
//                                                                     *
//     http://www.apache.org/licenses/LICENSE-2.0                      *
//                                                                     *
// Unless required by applicable law or agreed to in writing, software *
// distributed under the License is distributed on an "AS IS" BASIS,   *
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or     *
// implied. See the License for the specific language governing        *
// permissions and limitations under the License.                      *
//                                                                     *
// Copyright 2023 Matthew R. Hennefarth                                *
//**********************************************************************

use crate::traits::FloatSciConst;
use num_traits::{Float, Zero};
use std::ops::{Add, Mul};

/// Evaluate an $n$-degree polynomial at a specific value $x$.
///
/// Evaluates an $n$-degree polynomial where the coefficients are in
/// reversed order. That is if $\text{coeffs}\[i\] = c_i$, then evaluate
/// $$
/// c_0x^n + c_1x^{n-1} + \ldots + c_n
/// $$
pub(crate) fn eval_poly<T>(x: T, coeffs: &[T]) -> T
where
    T: Copy + Zero + Mul<Output = T>,
    <T as Mul>::Output: Add<T>,
{
    match coeffs.len() {
        0 => T::zero(),
        1 => coeffs[0],
        _ => coeffs[1..]
            .iter()
            .fold(coeffs[0], |result, &c| (result * x) + c),
    }
}

#[inline]
pub(crate) fn euler_reflection_prefactor<T>(x_abs: T, x_floor: T) -> T
where
    T: Float + FloatSciConst,
{
    let z = if (x_abs - x_floor) > (T::one() + T::one()).recip() {
        (x_floor + T::one()) - x_abs
    } else {
        x_abs - x_floor
    };
    x_abs * (T::PI() * z).sin()
}

#[cfg(test)]
mod tests {
    use super::*;

    const PRECISION: f64 = 1E-14;

    #[test]
    fn test_eval_poly() {
        assert_eq!(eval_poly(1.0, &[1.0, 1.0]), 2.0);
        assert_eq!(eval_poly(0.0, &[1.0, 1.0]), 1.0);
        assert_eq!(eval_poly(2.0, &[1.0, 1.0]), 3.0);

        for i in 0..10 {
            let i = i as f64;
            for j in 0..10 {
                let j = j as f64;
                assert_eq!(eval_poly(i, &[j]), j);
                assert_eq!(eval_poly(i, &[0.0, j]), j);
                assert_eq!(eval_poly(i, &[0.0, 0.0, j]), j);
            }
        }
        assert_almost_eq!(
            eval_poly(72.2, &[-6.42, 5.111219, 0.12]),
            -33097.2827882,
            PRECISION
        );

        assert_almost_eq!(
            eval_poly(-6.124, &[0.615, -2.801, 0.837, -4.701, 7.357]),
            1575.84328434321037,
            PRECISION
        );
    }
}
