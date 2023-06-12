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

use crate::constants::f64::{E, LOG_PI, PI, SQRT_TAU};
use num_traits::{One, Zero};
use std::ops::{Add, Mul};

/// Evaluate an $n$-degree polynomial at a specific value $x$.
///
/// Evaluates an $n$-degree polynomial where the coefficients are in
/// reversed order. That is if $\text{coeffs}\[i\] = c_i$, then evaluate
/// $$
/// c_0x^n + c_1x^{n-1} + \ldots + c_n
/// $$
fn eval_poly<T>(x: T, coeffs: &[T]) -> T
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
fn stirling_series(x: f64) -> f64 {
    // Taken from OEIS: A001164
    // Values pre-computed in rust
    const STIR_COEFFICIENTS: [f64; 5] = [
        7.84039221720066615423E-4,  // 163879/209018880
        -2.29472093621399167830E-4, // -571/2488320
        -2.68132716049382727186E-3, // -139/51840
        3.47222222222222202948E-3,  // 1/288
        8.33333333333333287074E-2,  // 1/12
    ];

    let series = 1.0 / x;
    let series = f64::one() + series * eval_poly(series, &STIR_COEFFICIENTS);
    let prefactor = (x / E).powf(x);
    SQRT_TAU / x.sqrt() * prefactor * series
}

#[inline]
fn euler_reflection_prefactor(x_abs: f64, x_floor: f64) -> f64 {
    let z = if (x_abs - x_floor) > 0.5 {
        (x_floor + 1.0) - x_abs
    } else {
        x_abs - x_floor
    };
    x_abs * (PI * z).sin()
}

/// The Gamma function for real-values arguments.
///
/// The Gamma function is defined as
/// $$
/// \Gamma(z) = \int^{\infty}_0 t^{z-1}e^{-t}dt
/// $$
/// where $\Re (z) > 0$. It is defined for the entire complex plane
/// through analytic continuation. It is a generalization of the
/// [crate::special::Factorial::factorial] function to integer values.
/// $$
/// \Gamma(x+1) = x!
/// $$
/// For a more thorough explanation of the Gamma function and all of its
/// properties, it is recommended to read through the [DLMF] or [wiki]
/// page.
///
/// # Examples
/// ```
/// use sci_rs::special::gamma;
/// assert_eq!(gamma(4.0), 6.0); // Gamma(4) = 3!
/// assert!(gamma(0.0).is_nan()); // Gamma(0) is undefined
/// assert!((gamma(4.5) - 11.6317283).abs() <  1e-5);
/// ```
/// ## Notes
/// The implementation uses a few different methods. Firstly, if
/// $|x| > 33$, then we utilize the Stirling series which is given by
/// $$
/// \sqrt{\frac{2\pi}{x}} \left(\frac{x}{e}\right)^x \left(1 + \frac{1}{12 x} + \frac{1}{288 x^2} - \frac{139}{51840 x^3} - \frac{571}{2488320 x^4} + \ldots \right)
/// $$
/// We additionally utilize the Euler's reflection formula for the Gamma
/// function to relate negative values to positive values.
/// $$
/// \Gamma(-z)\Gamma(z) = -\frac{\pi}{z\sin\pi z}
/// $$
/// Otherwise we recursively put the value into the range of $(2,3)$ using
/// $$
/// \Gamma(z+1) =z\Gamma(z)
/// $$
/// Then we use 2 rational functions of degree 6 and 7 to approximate the
/// Gamma function in this interval. This implementation is based off of
/// the implementation of Scipy which comes from the
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
pub fn gamma(x: f64) -> f64 {
    fn gamma_singularity() -> f64 {
        return f64::NAN;
    }

    if x.is_zero() {
        return gamma_singularity();
    }
    if !x.is_finite() {
        return x;
    }

    const MIN_TO_USE_STIRLING: f64 = 33.0;

    let x_abs = x.abs();
    if x_abs > MIN_TO_USE_STIRLING {
        if x.is_sign_positive() {
            return stirling_series(x);
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
        return PI / (z.abs() * stirling_series(x_abs)) * if is_positive_sign { 1.0 } else { -1.0 };
    }

    fn small(x: f64, z: f64) -> f64 {
        if x.is_zero() {
            return gamma_singularity();
        } else {
            z / ((1.0 + 0.5772156649015329 * x) * x)
        }
    }

    let mut x = x;
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

    const P: [f64; 7] = [
        1.60119522476751861407E-4,
        1.19135147006586384913E-3,
        1.04213797561761569935E-2,
        4.76367800457137231464E-2,
        2.07448227648435975150E-1,
        4.94214826801497100753E-1,
        9.99999999999999996796E-1,
    ];
    const Q: [f64; 8] = [
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

/// Natural logarithm of the Gamma function for
/// real-valued arguments.
///
/// Natural logarithm of the absolute value of the [gamma] function for
/// real-valued arguments.
///
/// $$
/// \ln \left|\Gamma(x)\right|
/// $$
///
/// # Examples
/// // TODO
///
/// ## Notes
/// Implementation is taken from the [cephes implementation] in the
/// Scipy package.
///
/// # References
/// - [cephes implementation]
///
/// [gamma]: crate::special::gamma()
/// [cephes implementation]: https://github.com/scipy/scipy/blob/main/scipy/special/cephes/gamma.c
pub fn gammaln(x: f64) -> f64 {
    fn gammaln_singularity() -> f64 {
        return f64::INFINITY;
    }
    if x < -34.0 {
        // Utilize Euler reflection and compute gammaln at positive value.
        let x_positive = -x;
        let p = x_positive.floor();
        if p == x_positive {
            return gammaln_singularity();
        }

        let z = euler_reflection_prefactor(x_positive, p);
        if z == 0.0 {
            return gammaln_singularity();
        }
        return LOG_PI - z.ln() - gammaln(x_positive);
    }

    if x < 13.0 {
        let mut z = 1.0;
        let mut x = x;

        while x >= 3.0 {
            x -= 1.0;
            z *= x;
        }

        while x < 2.0 {
            if x == 0.0 {
                return gammaln_singularity();
            }
            z /= x;
            x += 1.0;
        }
        z = z.abs();
        if x == 2.0 {
            return z.ln();
        }

        x -= 2.0;
        const B: [f64; 6] = [
            -1.37825152569120859100E3,
            -3.88016315134637840924E4,
            -3.31612992738871184744E5,
            -1.16237097492762307383E6,
            -1.72173700820839662146E6,
            -8.53555664245765465627E5,
        ];
        const C: [f64; 7] = [
            1.00000000000000000000E0,
            -3.51815701436523470549E2,
            -1.70642106651881159223E4,
            -2.20528590553854454839E5,
            -1.13933444367982507207E6,
            -2.53252307177582951285E6,
            -2.01889141433532773231E6,
        ];

        let p = x * eval_poly(x, &B) / eval_poly(x, &C);
        return z.ln() + p;
    }

    x
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::constants::f64::SQRT_PI;
    use crate::special::Factorial;

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
        assert_almost_eq!(gamma(PI), 2.2880377953400324179, PRECISION); // OEIS: A269545

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
            gamma(PI / 2.0),
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

    #[test]
    fn test_gammaln() {
        for i in -10..0 {
            assert_eq!(gammaln(i as f64), f64::INFINITY);
        }

        assert_almost_eq!(gammaln(12.5), 18.73434751193644570163, PRECISION)
    }
}
