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
use num_traits::{Float, One, Zero};
use std::ops::{Add, Mul, Sub};

/// Evaluate an $n$-degree polynomial at a specific value $x$.
///
/// Evaluates an $n$-degree polynomial where the coefficients are in reversed order. That is if $\text{coeffs}\[i\] = c_i$, then evaluate
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

/// Evaluate the Chebyshev series at $x$.
///
/// The Chebyshev series is defined as
/// $$
/// y = \sum_{i=0}^{N-1} c_{N-1-i} T_i (x)
/// $$
/// where $T_i(x)$ are Chebyshev polynomials of the first kind ([wiki]). We evaluate this series using the recursive [Clenshaw algorithm]. Note that the coefficients are stored in the reverse order such that the constant term is stored in the last element of the slice. This is similar to [eval_poly].
///
/// The Chebyshev polynomials are defined only on the interval $(-1, 1\$. Therefore, if the coefficients are for the interval $(a, b)$, then $x$ must be transformed as
/// $$
/// x \rightarrow 2 \frac{2x-b-a}{b-a}
/// $$
/// before entering the function. This affine mapping sends the interval $(a,b)\rightarrow (-1,1)$.
///
/// If the coefficients are for the inverted interval, in which $(a,b)$ is mapped to $(1/a, 1/b)$, then the affine transformation required is
/// $$
/// x\rightarrow 2\frac{2\frac{ab}{x} - b - a}{b-a}
/// $$
///
/// ## Note
/// This differs from the `chbevl` in the [cephes library] since it does not compute the series at $x/2$. Additionally the constant term in that series is always divided by two.
///
/// [wiki]: https://en.wikipedia.org/wiki/Chebyshev_polynomials
/// [Clenshaw algorithm]: https://en.wikipedia.org/wiki/Clenshaw_algorithm#Special_case_for_Chebyshev_series
/// [eval_poly]: crate::special::gamma_util::eval_poly
/// [cephes library]: https://github.com/scipy/scipy/blob/46081a85c3a6ca4c45610f4207abf791985e17e0/scipy/special/cephes/chbevl.c#L63
pub(crate) fn eval_cheby<T>(x: T, coeffs: &[T]) -> T
where
    T: One + Zero + Copy + Sub<Output = T>,
{
    let two = T::one() + T::one();
    match coeffs.len() {
        0 => T::zero(),
        1 => coeffs[0],
        2 => coeffs[1] + coeffs[0] * x,
        _ => {
            // This uses the Clenshaw's recurrece formula
            let mut bk = coeffs[0];
            let mut bk1 = T::zero();
            let mut _bk2 = T::zero(); // Compiler warns about unused, but it is used.

            let last_element = coeffs.len() - 1;

            for p in &coeffs[1..last_element] {
                _bk2 = bk1;
                bk1 = bk;
                bk = *p + two * x * bk1 - _bk2;
            }
            coeffs[last_element] + x * bk - bk1
        }
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

/// Coefficients are
/// $$
/// \frac{B_{2n}}{2n(2n-1)}
/// $$
/// where $B_{2n}$ is the $2n$th Bernoulli number.
pub(crate) trait LnGammaStirlingConsts: Sized {
    const LNGAMMA_STIRLING_COEFFS: [Self; 8];
}

macro_rules! impl_lngammastirlingconsts {
    ($($T: ty)*) => ($(
        impl LnGammaStirlingConsts for $T {
            const LNGAMMA_STIRLING_COEFFS: [Self; 8] = [
                -2.955065359477124183e-2, 6.4102564102564102564e-3,
                -1.9175269175269175269e-3, 8.4175084175084175084e-4,
                -5.952380952380952381e-4, 7.9365079365079365079e-4,
                -2.7777777777777777778e-3, 8.3333333333333333333e-2
            ];
        }
)*)
}

impl_lngammastirlingconsts! {f32 f64}

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

    #[test]
    fn test_eval_cheby() {
        for i in 0..5 {
            let i = i as f32;
            assert_eq!(eval_cheby(i, &[1.0, 1.0]), i + 1.0);
            assert_eq!(eval_cheby(i, &[2.0, 1.0]), 2.0 * i + 1.0);
            assert_eq!(eval_cheby(i, &[2.0, 3.0]), 2.0 * i + 3.0);
        }
        for i in 0..10 {
            let i = i as f64;
            for j in 0..10 {
                let j = j as f64;
                assert_eq!(eval_cheby(i, &[j]), j);
                assert_eq!(eval_cheby(i, &[0.0, j]), j);
                assert_eq!(eval_cheby(i, &[0.0, 0.0, j]), j);
                assert_eq!(eval_cheby(i, &[1.0, 0.0, j]), 2.0 * i * i - 1.0 + j);
                assert_eq!(
                    eval_cheby(i, &[1.0, 0.0, 0.0, j]),
                    4.0 * i * i * i - 3.0 * i + j
                );
            }
        }

        assert_almost_eq!(
            eval_cheby(0.25, &[-6.42, 5.111219, 0.12, 0.0]),
            -0.0285666250000001,
            PRECISION
        );

        assert_almost_eq!(
            eval_cheby(0.25, &[-6.42, 5.111219, 0.12, 0.5]),
            0.4714333749999999,
            PRECISION
        );

        assert_almost_eq!(
            eval_cheby(-0.325, &[-6.42, 48.2, -0.2212, 1.1, -0.2]),
            38.4254039375000019,
            PRECISION
        );
    }
}
