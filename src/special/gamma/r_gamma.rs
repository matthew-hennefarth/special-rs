//**********************************************************************
// This file is part of sci-rs                                         *
// Copyright 2023 Matthew R. Hennefarth                                *
//**********************************************************************

use crate::special::gamma::gamma_util::is_gamma_pole;
use crate::special::gamma::gamma_util::{
    euler_reflection_prefactor, gamma_stirling_series, StirlingSeriesCoefficients,
};
use crate::special::tools::eval_poly;
use crate::traits::FloatSciConst;
use num_traits::{cast, Float};
use std::ops::{AddAssign, DivAssign, MulAssign, SubAssign};

pub(crate) trait RealGammaConsts: Sized {
    const MIN_TO_USE_STIRLING: Self;
    const P: [Self; 7];
    const Q: [Self; 8];
    const INTERVAL: [Self; 2];
    const MIN_USE_SMALL: Self;
}

macro_rules! impl_realgammaconsts {
    ($($T: ty)*) => ($(
        impl RealGammaConsts for $T {
            const MIN_TO_USE_STIRLING: Self = 33.0;
            const P: [Self; 7] = [
                1.60119522476751861407E-4,
                1.19135147006586384913E-3,
                1.04213797561761569935E-2,
                4.76367800457137231464E-2,
                2.07448227648435975150E-1,
                4.94214826801497100753E-1,
                9.99999999999999996796E-1,
            ];
            const Q: [Self; 8] = [
                -2.31581873324120129819E-5,
                5.39605580493303397842E-4,
                -4.45641913851797240494E-3,
                1.18139785222060435552E-2,
                3.58236398605498653373E-2,
                -2.34591795718243348568E-1,
                7.14304917030273074085E-2,
                1.00000000000000000320E0,
            ];
            const INTERVAL: [Self; 2] = [2.0, 3.0];
            const MIN_USE_SMALL: Self = 1.0E-7;
        }
)*)
}

impl_realgammaconsts! {f32 f64}

/// What to return if we are at a singularity of the real-valued Gamma function
#[inline(always)]
fn gamma_singularity<T>() -> T
where
    T: Float,
{
    return T::nan();
}

/// Value of the Gamma function when |x| < `RealGammaConsts::MIN_TO_USE_SMALL`
///
/// We use the Laurant series of $\Gamma(z)$ around $z = 0$.
/// $$
/// \frac{1}{z} - \gamma + \left(\frac{\pi^2}{6} + \gamma^2\right)\frac{z}{2}
/// $$
#[inline(always)]
fn r_gamma_small<T>(x: T) -> T
where
    T: Float + FloatSciConst,
{
    x.recip() - T::GAMMA()
        + x / (T::one() + T::one())
            * (T::PI() * T::PI() / cast::<f32, T>(6.0).unwrap() + T::GAMMA() * T::GAMMA())
}

/// Implementation of the `gamma()` for real-valued inputs.
/// $$
/// \Gamma(x) = \int^{\infty}_0 t^{x-1}e^{-t}dt
/// $$
/// where $x$ is real-valued.
pub(crate) fn r_gamma<T>(x: T) -> T
where
    T: Float
        + Sized
        + RealGammaConsts
        + StirlingSeriesCoefficients
        + FloatSciConst
        + SubAssign
        + MulAssign
        + DivAssign
        + AddAssign,
{
    if !x.is_finite() {
        return x;
    }
    if is_gamma_pole(x) {
        return gamma_singularity();
    }
    if x.is_sign_negative() {
        return euler_reflection_prefactor(x) / r_gamma(-x);
    }

    if x > T::MIN_TO_USE_STIRLING {
        return gamma_stirling_series(x);
    }

    if x < T::MIN_USE_SMALL {
        return r_gamma_small(x);
    }

    let mut x = x;
    let mut z = T::one();

    while x >= T::INTERVAL[1] {
        x -= T::one();
        z *= x;
    }
    while x < T::INTERVAL[0] {
        z /= x;
        x += T::one();
    }

    if x == T::INTERVAL[0] {
        return z;
    }

    x -= T::INTERVAL[0];
    let p = eval_poly(x, &T::P);
    let q = eval_poly(x, &T::Q);

    z * p / q
}

#[cfg(test)]
mod tests {
    use super::*;

    use crate::constants::f64::{PI, SQRT_PI};
    use crate::special::Factorial;

    const PRECISION: f64 = 1E-14;

    #[test]
    fn test_r_gamma() {
        for i in 1..10 {
            assert_eq!(r_gamma(i as f64), (i - 1).factorial() as f64);
            assert!(r_gamma(-i as f64).is_nan());
        }
        assert!(r_gamma(0.0_f32).is_nan());
        assert!(r_gamma(f64::NAN).is_nan());

        const KNOWN_VALUES: [[f64; 2]; 28] = [
            // Test the half-integers
            [-2.5, -8.0 / 15.0 * SQRT_PI], // OEIS: A019707
            [-1.5, (4.0 / 3.0) * SQRT_PI], // OEIS: A245886
            [-0.5, -2.0 * SQRT_PI],        // OEIS: A245887
            [0.5, SQRT_PI],                // OEIS: A002161
            [1.5, SQRT_PI / 2.0],          // OEIS: A019704
            [2.5, 0.75 * SQRT_PI],         // OEIS: A245884
            [3.5, 15.0 / 8.0 * SQRT_PI],   // OEIS: A245885
            [4.5, 105.0 / 16.0 * SQRT_PI],
            [5.5, 945.0 / 32.0 * SQRT_PI],
            // Rational Values
            [1.0 / 3.0, 2.6789385347077476337], // OEIS: A073005
            [1.0 / 4.0, 3.6256099082219083119], // OEIS: A068466
            [1.0 / 5.0, 4.5908437119988030532], // OEIS: A175380
            [1.0 / 6.0, 5.5663160017802352043], // OEIS: A175379
            [1.0 / 7.0, 6.5480629402478244377], // OEIS: A220086
            [1.0 / 8.0, 7.5339415987976119047], // OEIS: A203142
            // Other Important Values
            [PI, 2.2880377953400324179], // OEIS: A269545
            [
                1.000001e-35,
                9.9999900000099999900000099999899999522784235098567139293e+34,
            ],
            [1.000001e-10, 9.99998999943278432519738283781e+9],
            [1.000001e-5, 99999.32279432557746387],
            [1.000001e-2, 99.43248512896257405886],
            [1.62123, 0.896081923385351],
            [-4.8, -0.062423361354759553],
            [0.1, 9.51350769866873183629],
            [1.0 - 1.0e-14, 1.00000000000000577215],
            [1.0 + 1.0e-14, 0.99999999999999422784],
            [PI / 2.0, 0.89056089038153932801],
            [5.0 - 1.0e-14, 23.999999999999652],
            [10.1, 454760.7514415855],
        ];

        for value in KNOWN_VALUES {
            assert_almost_eq!(r_gamma(value[0]), value[1], PRECISION);
        }
    }
}
