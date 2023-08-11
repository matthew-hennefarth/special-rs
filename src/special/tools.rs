//**********************************************************************
// This file is part of sci-rs                                         *
// Copyright 2023 Matthew R. Hennefarth                                *
//**********************************************************************

use num_traits::{Float, One, PrimInt, Zero};
use std::ops::{Add, Mul, Sub};

/// Evaluate an $n$-degree polynomial at a specific value $x$.
///
/// Evaluates an $n$-degree polynomial where the coefficients are in reversed order. That is if $\text{coeffs}\[i\] = c_i$, then evaluate
/// $$
/// c_0x^n + c_1x^{n-1} + \ldots + c_n
/// $$
pub(crate) fn eval_poly<T, Scalar>(x: T, coeffs: &[Scalar]) -> T
where
    T: Copy + One + Zero + Mul<T> + Mul<Scalar, Output = T> + Add<Scalar, Output = T>,
    Scalar: Copy,
{
    match coeffs.len() {
        0 => T::zero(),
        1 => T::one() * coeffs[0],
        _ => coeffs[1..]
            .iter()
            .fold(T::one() * coeffs[0], |result, &c| (result * x) + c),
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
pub(crate) fn eval_cheby<T, Scalars>(x: T, coeffs: &[Scalars]) -> T
where
    T: One + Zero + Copy + Sub<Output = T> + Mul<Scalars, Output = T> + Add<Scalars, Output = T>,
    Scalars: Copy,
{
    let two = T::one() + T::one();
    match coeffs.len() {
        0 => T::zero(),
        1 => T::one() * coeffs[0],
        2 => x * coeffs[0] + coeffs[1],
        _ => {
            // This uses the Clenshaw's recurrece formula
            let mut bk = T::one() * coeffs[0];
            let mut bk1 = T::zero();
            let mut _bk2 = T::zero(); // Compiler warns about unused, but it is used.

            let last_element = coeffs.len() - 1;

            for p in &coeffs[1..last_element] {
                _bk2 = bk1;
                bk1 = bk;
                bk = two * x * bk1 - _bk2 + *p;
            }
            x * bk - bk1 + coeffs[last_element]
        }
    }
}

pub(crate) fn frexp<T, I>(x: T) -> (T, I)
where
    T: Float,
    I: PrimInt,
{
    if x.is_zero() {
        return (x, I::zero());
    } else {
        let lg = x.abs().log2();
        let s = (lg - lg.floor() - T::one()).exp2();
        let exp = lg.floor() + T::one();
        (x.signum() * s, I::from(exp).unwrap())
    }
}

pub(crate) fn ldexp<T>(num: T, exp: i64) -> T
where
    T: Float,
{
    if exp >= 0 {
        num * T::from(2u64.pow(exp as u32)).unwrap()
    } else {
        num / T::from(2u64.pow(-exp as u32)).unwrap()
    }
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
