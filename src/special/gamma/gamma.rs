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

use crate::special::gamma::{euler_reflection_prefactor, eval_poly};
use crate::traits::FloatSciConst;
use num_traits::{Float, Zero};

/// Implementation of the Gamma and related functions for both real and complex-valued inputs.
pub trait Gamma {
    /// The Gamma function.
    ///
    /// The Gamma function is defined as
    /// $$
    /// \Gamma(z) = \int^{\infty}_0 t^{z-1}e^{-t}dt
    /// $$
    /// where $\Re (z) > 0$. It is defined for the entire complex plane through analytic continuation. It is a generalization of the
    /// [factorial] function to integer values.
    /// $$
    /// \Gamma(n+1) = n!
    /// $$
    /// For a more thorough explanation of the Gamma function and all of its properties, it is recommended to read through the [DLMF] or [wiki] page.
    ///
    /// # Examples
    /// ```
    /// use sci_rs::special::gamma;
    /// assert_eq!(gamma(4.0_f32), 6.0); // Gamma(4) = 3!
    /// assert!(gamma(0.0_f64).is_nan()); // Gamma(0) is undefined
    /// assert!((gamma(4.5_f64) - 11.6317283).abs() <  1e-5);
    /// ```
    /// ## Notes
    /// The implementation uses a few different methods. Firstly, for real-valued arguments ($z=x$) and $|x| > 33$, then we utilize the Stirling series which is given by
    /// $$
    /// \sqrt{\frac{2\pi}{x}} \left(\frac{x}{e}\right)^x \left(1 + \frac{1}{12 x} + \frac{1}{288 x^2} - \frac{139}{51840 x^3} - \frac{571}{2488320 x^4} + \ldots \right)
    /// $$
    /// We additionally utilize the Euler's reflection formula for the Gamma function to relate negative values to positive values.
    /// $$
    /// \Gamma(-x)\Gamma(x) = -\frac{\pi}{x\sin\pi x}
    /// $$
    /// Otherwise we recursively put the value into the range of $(2,3)$ using
    /// $$
    /// \Gamma(x+1) =x\Gamma(x)
    /// $$
    /// Then we use 2 rational functions of degree 6 and 7 to approximate the Gamma function in this interval. This implementation is based off of the implementation of Scipy which comes from the
    /// [cephes implementation].
    ///
    /// # References
    /// - [DLMF]
    /// - [wiki]
    /// - [cephes implementation]
    ///
    /// [comment]: <> (Reference hyperlinks)
    /// [DLMF]: https://dlmf.nist.gov/5.2
    /// [wiki]: https://en.wikipedia.org/wiki/Gamma_function
    /// [cephes implementation]: https://github.com/scipy/scipy/blob/main/scipy/special/cephes/gamma.c
    /// [factorial]: crate::special::Factorial::factorial
    fn gamma(self) -> Self;
}

/// Gamma function evaluated at $z$.
///
/// Has the same semantics as [gamma] in the [Gamma trait].
///
/// # Examples
/// ```
/// use sci_rs::special::{Gamma, gamma};
/// assert_eq!(gamma(1.0_f32), 1.0_f32.gamma());
/// assert_eq!(gamma(2.24_f64), 2.24_f64.gamma());
/// ```
///
/// [Gamma trait]: crate::special::Gamma
/// [gamma]: crate::special::Gamma::gamma
#[inline(always)]
pub fn gamma<T>(z: T) -> T
where
    T: Gamma,
{
    z.gamma()
}

macro_rules! float_gamma_impl {
    ($($T: ty)*) => ($(
        impl Gamma for $T {
            fn gamma(self) -> Self {
                fn gamma_singularity() -> $T {
                    return <$T>::NAN;
                }

                if self.is_zero() {
                    return gamma_singularity();
                }
                if !self.is_finite() {
                    return self;
                }

                const MIN_TO_USE_STIRLING: $T = 33.0;

                let x_abs = self.abs();
                if x_abs > MIN_TO_USE_STIRLING {
                    if self.is_sign_positive() {
                        return stirling_series(self);
                    }

                    let x_abs_floor = x_abs.floor();
                    // Gamma function has poles at the negative integers
                    if x_abs_floor == x_abs {
                        return gamma_singularity();
                    }

                    let is_positive_sign = (x_abs_floor as usize) % 2 == 1;

                    // Utilize the Euler's reflection formula for the gamma function
                    // Gamma(-z)Gamma(z) = -\frac{\pi}{z\sin\pi z}
                    let z = euler_reflection_prefactor(x_abs, x_abs_floor);
                    if z.is_zero() {
                        return gamma_singularity();
                    }
                    return Self::PI() / (z.abs() * stirling_series(x_abs)) * if is_positive_sign { 1.0 } else { -1.0 };
                }

                fn small(x: $T, z: $T) -> $T {
                    if x.is_zero() {
                        return gamma_singularity();
                    } else {
                        z / ((1.0 + 0.5772156649015329 * x) * x)
                    }
                }

                let mut x = self;
                let mut z = 1.0;

                while x >= 3.0 {
                    x -= 1.0;
                    z *= x;
                }

                while x.is_sign_negative() {
                    if x > -1.0E-9 {
                        return small(x, z);
                    }
                    z /= x;
                    x += 1.0;
                }
                while x < 2.0 {
                    if x < 1.0E-9 {
                        return small(x, z);
                    }
                    z /= x;
                    x += 1.0;
                }

                if x == 2.0 {
                    return z;
                }

                const P: [$T; 7] = [
                    1.60119522476751861407E-4,
                    1.19135147006586384913E-3,
                    1.04213797561761569935E-2,
                    4.76367800457137231464E-2,
                    2.07448227648435975150E-1,
                    4.94214826801497100753E-1,
                    9.99999999999999996796E-1,
                ];
                const Q: [$T; 8] = [
                    -2.31581873324120129819E-5,
                    5.39605580493303397842E-4,
                    -4.45641913851797240494E-3,
                    1.18139785222060435552E-2,
                    3.58236398605498653373E-2,
                    -2.34591795718243348568E-1,
                    7.14304917030273074085E-2,
                    1.00000000000000000320E0,
                ];

                let p = eval_poly(x - 2.0, &P);
                let q = eval_poly(x - 2.0, &Q);

                z * p / q
            }
        }
    )*)
}

float_gamma_impl! {f32 f64}

/// Trait to tag types which have stirling coefficients expansions
/// Will just be f32 and f64, but I don't want to copy and paste.
trait StirlingSeriesCoefficients: Sized {
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
/// Compute the Stirling series for a given real-valued $x$.
/// $$
/// \sqrt{\frac{2\pi}{x}} \left(\frac{x}{e}\right)^n \left(1 + \frac{1}{12 x} + \frac{1}{288 x^2} - \ldots \right)
/// $$
/// See [here](https://dlmf.nist.gov/5.11) for a detailed explanation of
/// the Stirling series and its relationship to the Gamma function.
///
/// ## Notes
/// The implementation expands to 6th order and the coefficients are taken from OEIS: [A001164] and [A001163]
///
/// [A001164]: https://oeis.org/A001164
/// [A001163]: https://oeis.org/A001163
fn stirling_series<T>(x: T) -> T
where
    T: Float + FloatSciConst + StirlingSeriesCoefficients,
{
    let series = x.recip();
    let series = T::one() + series * eval_poly(series, &T::STIR_COEFFICIENTS);
    let prefactor = (x / T::E()).powf(x);
    T::SQRT_TAU() / x.sqrt() * prefactor * series
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::constants::f64::SQRT_PI;
    use crate::special::Factorial;

    const PRECISION: f64 = 1E-14;

    #[test]
    fn test_stirlings_series() {
        const REFERENCE_VALUES: [f64; 10] = [
            1.0002224601164145,
            1.0000024896493827,
            2.0000002868007112,
            6.0000001000594825,
            24.0000000672158009,
            120.0000000677005829,
            720.0000000819628667,
            5040.0000000665886546,
            40319.9999996674159775,
            362879.9999961055582389,
        ];
        const REFERENCE_140: f64 = 961572319694109.0E224;
        const REFERENCE_MAXSTIR: f64 = 2919114949633048.0E230;
        const REFERENCE_MAXSTIR_P_EPSILON: f64 = 29191294268940784.0E229;

        const REFERENCE_150: f64 = 3808922637630618.0E245;
        for i in 0..10 {
            assert_almost_eq!(
                stirling_series((i + 1) as f64),
                REFERENCE_VALUES[i],
                PRECISION
            );
        }
        assert_almost_eq!(stirling_series(140.0), REFERENCE_140, PRECISION);
        assert_almost_eq!(stirling_series(143.01608), REFERENCE_MAXSTIR, PRECISION);
        assert_almost_eq!(
            stirling_series(143.016081),
            REFERENCE_MAXSTIR_P_EPSILON,
            PRECISION
        );
        assert_almost_eq!(stirling_series(150.0), REFERENCE_150, PRECISION);
    }

    #[test]
    fn test_gamma() {
        for i in 1..10 {
            assert_eq!(gamma(i as f64), (i - 1).factorial() as f64);
            assert!(gamma(-i as f64).is_nan());
        }
        assert!(gamma(0.0).is_nan());
        assert!(gamma(f64::NAN).is_nan());

        // Test the half-integers
        assert_almost_eq!(gamma(-2.5), -8.0 / 15.0 * SQRT_PI, PRECISION); // OEIS: A019707
        assert_almost_eq!(gamma(-1.5), (4.0 / 3.0) * SQRT_PI, PRECISION); // OEIS: A245886
        assert_almost_eq!(gamma(-0.5), -2.0 * SQRT_PI, PRECISION); // OEIS: A245887
        assert_almost_eq!(gamma(0.5), SQRT_PI, PRECISION); // OEIS: A002161
        assert_almost_eq!(gamma(1.5), SQRT_PI / 2.0, PRECISION); // OEIS: A019704
        assert_almost_eq!(gamma(2.5), 0.75 * SQRT_PI, PRECISION); // OEIS: A245884
        assert_almost_eq!(gamma(3.5), 15.0 / 8.0 * SQRT_PI, PRECISION); // OEIS: A245885
        assert_almost_eq!(gamma(4.5), 105.0 / 16.0 * SQRT_PI, PRECISION);
        assert_almost_eq!(gamma(5.5), 945.0 / 32.0 * SQRT_PI, PRECISION);

        // Rational Values
        assert_almost_eq!(gamma(1.0 / 3.0), 2.6789385347077476337, PRECISION); // OEIS: A073005
        assert_almost_eq!(gamma(0.25), 3.6256099082219083119, PRECISION); // OEIS: A068466
        assert_almost_eq!(gamma(0.2), 4.5908437119988030532, PRECISION); // OEIS: A175380
        assert_almost_eq!(gamma(1.0 / 6.0), 5.5663160017802352043, PRECISION); // OEIS: A175379
        assert_almost_eq!(gamma(1.0 / 7.0), 6.5480629402478244377, PRECISION); // OEIS: A220086
        assert_almost_eq!(gamma(1.0 / 8.0), 7.5339415987976119047, PRECISION); // OEIS: A203142

        // Other Important Values
        assert_almost_eq!(gamma(f64::PI()), 2.2880377953400324179, PRECISION); // OEIS: A269545

        assert_almost_eq!(
            gamma(1.000001e-35),
            9.9999900000099999900000099999899999522784235098567139293e+34,
            PRECISION
        );
        assert_almost_eq!(
            gamma(1.000001e-10),
            9.99998999943278432519738283781280989934496494539074049002e+9,
            PRECISION
        );
        assert_almost_eq!(gamma(1.000001e-5), 99999.32279432557746387, PRECISION);
        assert_almost_eq!(gamma(1.000001e-2), 99.43248512896257405886, PRECISION);
        assert_almost_eq!(gamma(1.62123), 0.896081923385351, PRECISION);

        assert_almost_eq!(gamma(-4.8), -0.062423361354759553, PRECISION);

        assert_almost_eq!(
            gamma(1.0e-5 + 1.0e-16),
            99999.42279322556767360213300482199406241771308740302819426480,
            1e-9
        );
        assert_almost_eq!(
            gamma(0.1),
            9.513507698668731836292487177265402192550578626088377343050000,
            1e-14
        );
        assert_almost_eq!(
            gamma(1.0 - 1.0e-14),
            1.000000000000005772156649015427511664653698987042926067639529,
            PRECISION
        );
        assert_almost_eq!(
            gamma(1.0 + 1.0e-14),
            0.99999999999999422784335098477029953441189552403615306268023,
            PRECISION
        );
        assert_almost_eq!(
            gamma(f64::PI() / 2.0),
            0.890560890381539328010659635359121005933541962884758999762766,
            PRECISION
        );

        assert_almost_eq!(gamma(5.0 - 1.0e-14), 23.999999999999652, PRECISION);

        assert_almost_eq!(gamma(10.1), 454760.7514415855, PRECISION);
        assert_almost_eq!(
            gamma(150.0 + 1.0e-12),
            3.8089226376496421386707466577615064443807882167327097140e+260,
            1e248
        );
    }
}
