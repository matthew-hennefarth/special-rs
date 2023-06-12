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

use crate::constants::{E, PI, SQRT_2_PI};
use num_traits::Zero;
use std::ops::{Add, Mul};

/// Evaluate an $n$-degree polynomial at a specific value $x$.
///
/// Evaluates an $n$-degree polynomial where the coefficients are in
/// reversed order. That is if $\text{coeffs}[i] = c_i$, then evaluate
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

/// Stirlings formula
///
/// Compute the Stirling series for a given real-valued $x$.
/// $$
/// \sqrt{\frac{2\pi}{x}} \left(\frac{x}{e}\right)^n \left(1 + \frac{1}{12 x} + \frac{1}{288 x^2} - \frac{139}{51840 x^3} - \frac{571}{2488320 x^4} + \ldots \right)
/// $$
/// See [here](https://dlmf.nist.gov/5.11) for a detailed explanation of
/// the Stirling series and its relationship to the Gamma function.
fn stirling_series(x: f64) -> f64 {
    const STIR_COEFFICIENTS: [f64; 5] = [
        7.84039221720066615423E-4,
        -2.29472093621399167830E-4,
        -2.68132716049382727186E-3,
        3.47222222222222202948E-3,
        8.33333333333333287074E-2,
    ];

    let series = 1.0 / x;
    let series = 1.0 + series * eval_poly(series, &STIR_COEFFICIENTS);
    let prefactor = (x / E).powf(x);
    SQRT_2_PI / x.sqrt() * prefactor * series
}

/// The Gamma function
///
/// The Gamma function is defined as
/// $$
/// \Gamma(z) = \int^{\infty}_0 t^{z-1}e^{-t}dt
/// $$
/// where $\Re (z) > 0$. It is defined for the entire complex plane
/// through analytic continuation.
///
/// # Examples
/// ```
/// use sci_rs::special::gamma;
/// assert_eq!(gamma(4.0), 6.0); // Gamma(4) = 3!
/// assert_eq!(gamma(0.0), f64::INFINITY); // Gamma(0) is undefined
/// assert!((gamma(4.5) - 11.6317283).abs() <  1e-5);
/// ```
/// ## Notes
/// The implementation uses a few different methods. Firstly, if
/// $|x| > 33$, then we utilize the Stirling series which is given by
/// $$
/// \sqrt{\frac{2\pi}{x}} \left(\frac{x}{e}\right)^x \left(1 + \frac{1}{12 x} + \frac{1}{288 x^2} - \frac{139}{51840 x^3} - \frac{571}{2488320 x^4} + \ldots \right)
/// $$
///
/// Otherwise we recursively put the value into the range of \(2,3\) using
/// $$
/// \Gamma(x+1) =x\Gamma(x)
/// $$
/// Then we use 2 ration functions of degree 6 and 7 to approximate the
/// Gamma function in this interval. This implementation is based off of
/// the implementation of Scipy which comes from cephes (see
/// [here](https://github.com/scipy/scipy/blob/main/scipy/special/cephes/gamma.c)).
pub fn gamma(x: f64) -> f64 {
    if x.is_zero() {
        return f64::INFINITY;
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
            return f64::INFINITY;
        }

        let is_positive_sign = (x_abs_floor as usize) % 2 == 1;

        // Utilize the Euler's reflection formula for the gamma function
        // Gamma(-z)Gamma(z) = -\frac{\pi}{z\sin\pi z}
        let z = if (x_abs - x_abs_floor) > 0.5 {
            (x_abs_floor + 1.0) - x_abs
        } else {
            x_abs - x_abs_floor
        };
        let z = x_abs * (PI * z).sin();
        if z.is_zero() {
            return if is_positive_sign {
                f64::INFINITY
            } else {
                f64::NEG_INFINITY
            };
        }
        return PI / (z.abs() * stirling_series(x_abs)) * if is_positive_sign { 1.0 } else { -1.0 };
    }

    let mut x = x;
    let mut z = 1.0;

    fn small(x: f64, z: f64) -> f64 {
        if x.is_zero() {
            f64::NAN
        } else {
            z / ((1.0 + 0.5772156649015329 * x) * x)
        }
    }

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
    x -= 2.0;

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

    let p = eval_poly(x, &P);
    let q = eval_poly(x, &Q);

    z * p / q
}

#[cfg(test)]
mod tests {
    use super::*;
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
            1.0002224601164145, 1.0000024896493827, 2.000000286800711, 6.0000001000594825,
            24.0000000672158, 120.00000006770058, 720.0000000819629, 5040.000000066589,
            40319.999999667416, 362879.99999610556,
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
        assert!(gamma(0.0).is_infinite());
        assert!(gamma(f64::NAN).is_nan());

        assert_almost_eq!(gamma(-4.8), -0.062423361354759553, PRECISION);
        assert_almost_eq!(gamma(-1.5), 2.3632718012073547030, PRECISION);
        assert_almost_eq!(gamma(-0.5), -3.544907701811032054, PRECISION);

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
        assert_almost_eq!(gamma(0.5), 1.7724538509055160, PRECISION);
        assert_almost_eq!(gamma(1.62123), 0.896081923385351, PRECISION);

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
        assert_almost_eq!(gamma(1.0), 1.0, 1e-15);
        assert_almost_eq!(
            gamma(1.0 + 1.0e-14),
            0.99999999999999422784335098477029953441189552403615306268023,
            PRECISION
        );
        assert_almost_eq!(
            gamma(1.5),
            0.886226925452758013649083741670572591398774728061193564106903,
            PRECISION
        );
        assert_almost_eq!(
            gamma(PI / 2.0),
            0.890560890381539328010659635359121005933541962884758999762766,
            PRECISION
        );
        assert_eq!(gamma(2.0), 1.0);
        assert_almost_eq!(
            gamma(2.5),
            1.329340388179137020473625612505858887098162092091790346160355,
            PRECISION
        );
        assert_almost_eq!(gamma(3.0), 2.0, 1e-14);
        assert_almost_eq!(
            gamma(PI),
            2.288037795340032417959588909060233922889688153356222441199380,
            PRECISION
        );
        assert_almost_eq!(
            gamma(3.5),
            3.323350970447842551184064031264647217745405230229475865400889,
            PRECISION
        );
        assert_almost_eq!(gamma(4.0), 6.0, 1e-13);
        assert_almost_eq!(
            gamma(4.5),
            11.63172839656744892914422410942626526210891830580316552890311,
            PRECISION
        );
        assert_almost_eq!(gamma(5.0 - 1.0e-14), 23.999999999999652, PRECISION);
        assert_almost_eq!(gamma(5.0), 24.0, 1e-12);
        assert_almost_eq!(
            gamma(5.0 + 1.0e-14),
            24.00000000000036146824042363510111050137786752408660789873592,
            PRECISION
        );
        assert_almost_eq!(
            gamma(5.5),
            52.34277778455352018114900849241819367949013237611424488006401,
            PRECISION
        );
        assert_almost_eq!(gamma(10.1), 454760.7514415855, PRECISION);
        assert_almost_eq!(
            gamma(150.0 + 1.0e-12),
            3.8089226376496421386707466577615064443807882167327097140e+260,
            1e248
        );
    }
}
