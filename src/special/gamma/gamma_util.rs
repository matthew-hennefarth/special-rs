//**********************************************************************
// This file is part of sci-rs                                         *
// Copyright 2023 Matthew R. Hennefarth                                *
//**********************************************************************

use crate::traits::FloatSciConst;

use num_complex::ComplexFloat;
use num_traits::{Float, FloatConst, One, Zero};
use std::ops::{Add, Div, Mul, Sub};

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

/// Determines if $z$ is at a pole of the Gamma function (0, -1, -2, etc).
#[inline]
pub(crate) fn is_gamma_pole<T>(z: T) -> bool
where
    T: ComplexFloat,
{
    z.re() <= T::Real::zero() && z.re() == z.re().floor() && z.im().is_zero()
}

/// Computes the pre-factor need to flip the sign of the arg in the Gamma function.
/// $$
/// \Gamma(z)\Gamma(-z) = - \frac{\pi}{z\sin(\pi z)}
/// $$
#[inline]
pub(crate) fn euler_reflection_prefactor<T>(z: T) -> T
where
    T: ComplexFloat + Mul<<T as ComplexFloat>::Real, Output = T>,
    <T as ComplexFloat>::Real: FloatSciConst,
{
    -(z * (z * T::Real::PI()).sin()).recip() * T::Real::PI()
}

/// Trait to tag types which have stirling coefficients expansions
/// Will just be f32 and f64, but I don't want to copy and paste.
pub(crate) trait StirlingSeriesCoefficients: Sized {
    const STIR_COEFFICIENTS: [Self; 5];
}

macro_rules! impl_stirseries_coefficients {
    ($($T: ty)*) => ($(
        impl StirlingSeriesCoefficients for $T {
            // Taken from OEIS: A001164
            // Values pre-computed in rust
            const STIR_COEFFICIENTS: [Self; 5] = [
                7.84039221720066615423E-4,  // 163879/209018880
                -2.29472093621399167830E-4, // -571/2488320
                -2.68132716049382727186E-3, // -139/51840
                3.47222222222222202948E-3,  // 1/288
                8.33333333333333287074E-2,  // 1/12
            ];
        }
)*)
}

impl_stirseries_coefficients! {f32 f64}

/// Stirlings Formula
///
/// Compute the Stirling series for a given value $z$.
/// $$
/// \sqrt{\frac{2\pi}{z}} \left(\frac{z}{e}\right)^n \left(1 + \frac{1}{12 z} + \frac{1}{288 z^2} - \ldots \right)
/// $$
/// See [here](https://dlmf.nist.gov/5.11) for a detailed explanation of
/// the Stirling series and its relationship to the Gamma function.
///
/// # Notes
/// The implementation expands to 6th order and the coefficients are taken from OEIS: [A001164] and [A001163]
///
/// [A001164]: https://oeis.org/A001164
/// [A001163]: https://oeis.org/A001163
pub(crate) fn gamma_stirling_series<T>(z: T) -> T
where
    T: ComplexFloat
        + Add<<T as ComplexFloat>::Real, Output = T>
        + Mul<<T as ComplexFloat>::Real, Output = T>
        + Div<<T as ComplexFloat>::Real, Output = T>,
    <T as ComplexFloat>::Real: FloatSciConst + StirlingSeriesCoefficients,
{
    let rz = z.recip();
    let series = T::one() + rz * eval_poly(rz, &T::Real::STIR_COEFFICIENTS);
    let prefactor = (z * z.ln() - z).exp();
    prefactor * series / z.sqrt() * T::Real::SQRT_TAU()
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

/// Stirling approximation for lngamma(z)
/// $$
/// \ln \Gamma(z) \approx \ln \left(\sqrt{2\pi} z^{z-0.5}e^{-z}\right) + \sum^{\infty}\_{n=1} \frac{B_{2n}}{2n(2n-1)z^{2n-1}}
/// $$
/// where the $B_{2n}$ are Bernoulli numbers
#[inline]
pub(crate) fn lngamma_stirling<T>(z: T) -> T
where
    T: ComplexFloat
        + Add<<T as ComplexFloat>::Real, Output = T>
        + Mul<<T as ComplexFloat>::Real, Output = T>,
    <T as ComplexFloat>::Real: FloatSciConst + LnGammaStirlingConsts,
{
    let rz: T = z.recip();
    let rzz: T = rz / z;

    // (z - 0.5)ln(z) - z + ln(sqrt(2pi)) + (1/z)(Stirling)
    let q = (z - (T::one() + T::one()).recip()) * z.ln() - z + T::Real::LOG_SQRT_2_PI();
    q + rz * eval_poly(rzz, &T::Real::LNGAMMA_STIRLING_COEFFS)
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

    #[test]
    fn test_stirlings_series() {
        const REFERENCE_VALUES: [f64; 10] = [
            1.0002224601164145, 1.0000024896493827, 2.0000002868007112, 6.0000001000594825,
            24.0000000672158009, 120.00000006770064, 720.0000000819621, 5040.0000000665805,
            40319.99999966742, 362879.9999961068,
        ];
        const REFERENCE_140: f64 = 961572319694071.10E224;
        const REFERENCE_MAXSTIR: f64 = 2919114949633263.6E230;
        const REFERENCE_MAXSTIR_P_EPSILON: f64 = 29191294268938814.0E229;

        const REFERENCE_150: f64 = 3808922637630237.0E245;
        for i in 0..10 {
            assert_almost_eq!(
                gamma_stirling_series((i + 1) as f64),
                REFERENCE_VALUES[i],
                PRECISION
            );
        }
        assert_almost_eq!(gamma_stirling_series(140.0), REFERENCE_140, PRECISION);
        assert_almost_eq!(
            gamma_stirling_series(143.01608),
            REFERENCE_MAXSTIR,
            PRECISION
        );
        assert_almost_eq!(
            gamma_stirling_series(143.016081),
            REFERENCE_MAXSTIR_P_EPSILON,
            PRECISION
        );
        assert_almost_eq!(gamma_stirling_series(150.0), REFERENCE_150, PRECISION);
    }
}
