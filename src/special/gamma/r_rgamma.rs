//**********************************************************************
// This file is part of sci-rs                                         *
// Copyright 2023 Matthew R. Hennefarth                                *
//**********************************************************************

use std::ops::{AddAssign, DivAssign, MulAssign, SubAssign};

use crate::special::gamma::gamma_util::{euler_reflection_prefactor, eval_cheby, is_gamma_pole};
use crate::special::gamma::real_gamma_impl::*;
use crate::traits::FloatSciConst;
use num_traits::Float;

pub(crate) trait RGammaConsts: Sized {
    const MIN_VALUE_FOR_EXP: Self;
    const R: [Self; 16];
}

macro_rules! impl_rgammaconsts {
    ($($T: ty)*) => ($(
        impl RGammaConsts for $T {
            const MIN_VALUE_FOR_EXP: Self = 34.84425627277176174;
            const R: [Self; 16] = [
                3.13173458231230000000E-17,
                -6.70718606477908000000E-16,
                2.20039078172259550000E-15,
                2.47691630348254132600E-13,
                -6.60074100411295197440E-12,
                5.13850186324226978840E-11,
                1.08965386454418662084E-9,
                -3.33964630686836942556E-8,
                2.68975996440595483619E-7,
                2.96001177518801696639E-6,
                -8.04814124978471142852E-5,
                4.16609138709688864714E-4,
                5.06579864028608725080E-3,
                -6.41925436109158228810E-2,
                -4.98558728684003594785E-3,
                6.37730078052619747675E-2,
            ];
        }
)*)
}

impl_rgammaconsts! {f32 f64}

/// Implementation of the `rgamma()` for real-valued inputs
/// $$
/// \frac{1}{\Gamma(x)}
/// $$
/// where $x$ is real-valued. For arguments between -MIN_VALUE_FOR_EXP, MIN_VALUE_FOR_EXP, we use recursion to move the argument to within (0,1). For negative integers and `0.0` we return `0.0`. We use a Chebyshev expansion in the interval of `(0,1)`. For arguments around the MIN_VALUE_FOR_EXP, we use
/// $$
/// \frac{1}{\Gamma(x)} = \exp{-\ln{\Gamma(x)}}
/// $$
///
/// Overflow and underflow can accor despite the function not having any singularities.
pub(crate) fn r_rgamma<T>(x: T) -> T
where
    T: Float
        + RGammaConsts
        + FloatSciConst
        + SubAssign
        + MulAssign
        + DivAssign
        + AddAssign
        + RealGammaLnConsts,
{
    if is_gamma_pole(x) {
        return T::zero();
    }

    if x > T::MIN_VALUE_FOR_EXP {
        return (-r_lgamma(x)).exp();
    }

    if x < -T::MIN_VALUE_FOR_EXP {
        let y = r_lgamma(-x) - euler_reflection_prefactor(x).abs().ln();
        return r_gammasgn(x) * y.exp();
    }

    let mut z = T::one();
    let mut w = x;

    while w > T::one() {
        w -= T::one();
        z *= w;
    }
    while w < T::zero() {
        z /= w;
        w += T::one();
    }
    if w == T::one() {
        return z.recip();
    }
    w * (T::one() + eval_cheby((T::one() + T::one()) * w - T::one(), &T::R)) / z
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::constants::f64::{PI, SQRT_PI};
    use crate::special::Factorial;

    const PRECISION: f64 = 1E-14;

    #[test]
    fn test_r_rgamma() {
        for i in 1..10 {
            assert_eq!(r_rgamma(i as f64), ((i - 1).factorial() as f64).recip());
            assert_eq!(r_rgamma(-i as f64), 0.0);
        }

        assert_eq!(r_rgamma(0.0), 0.0);
        assert!(r_rgamma(f64::NAN).is_nan());

        // Test the half-integers
        assert_almost_eq!(r_rgamma(-2.5), -15.0 / 8.0 / SQRT_PI, PRECISION); // OEIS: A019707
        assert_almost_eq!(r_rgamma(-1.5), 0.75 / SQRT_PI, PRECISION); // OEIS: A245886
        assert_almost_eq!(r_rgamma(-0.5), -0.5 / SQRT_PI, PRECISION); // OEIS: A245887
        assert_almost_eq!(r_rgamma(0.5), SQRT_PI.recip(), PRECISION); // OEIS: A002161
        assert_almost_eq!(r_rgamma(1.5), 2.0 / SQRT_PI, PRECISION); // OEIS: A019704
        assert_almost_eq!(r_rgamma(2.5), 4.0 / 3.0 / SQRT_PI, PRECISION); // OEIS: A245884
        assert_almost_eq!(r_rgamma(3.5), 8.0 / 15.0 / SQRT_PI, PRECISION); // OEIS: A245885
        assert_almost_eq!(r_rgamma(4.5), 16.0 / 105.0 / SQRT_PI, PRECISION);
        assert_almost_eq!(r_rgamma(5.5), 32.0 / 945.0 / SQRT_PI, PRECISION);

        // Rational Values
        assert_almost_eq!(
            r_rgamma(1.0 / 3.0),
            2.6789385347077476337.recip(),
            PRECISION
        ); // OEIS: A073005
        assert_almost_eq!(r_rgamma(0.25), 3.6256099082219083119.recip(), PRECISION); // OEIS: A068466
        assert_almost_eq!(r_rgamma(0.2), 4.5908437119988030532.recip(), PRECISION); // OEIS: A175380
        assert_almost_eq!(
            r_rgamma(1.0 / 6.0),
            5.5663160017802352043.recip(),
            PRECISION
        ); // OEIS: A175379
        assert_almost_eq!(
            r_rgamma(1.0 / 7.0),
            6.5480629402478244377.recip(),
            PRECISION
        ); // OEIS: A220086
        assert_almost_eq!(
            r_rgamma(1.0 / 8.0),
            7.5339415987976119047.recip(),
            PRECISION
        ); // OEIS: A203142

        // Other Important Values
        assert_almost_eq!(r_rgamma(PI), 2.2880377953400324179.recip(), PRECISION); // OEIS: A269545

        assert_almost_eq!(
            r_rgamma(1.000001e-35),
            9.9999900000099999900000099999899999522784235098567139293e+34.recip(),
            PRECISION
        );
        assert_almost_eq!(
            r_rgamma(1.000001e-10),
            9.99998999943278432519738283781280989934496494539074049002e+9.recip(),
            PRECISION
        );
        assert_almost_eq!(
            r_rgamma(1.000001e-5),
            99999.32279432557746387.recip(),
            PRECISION
        );
        assert_almost_eq!(
            r_rgamma(1.000001e-2),
            99.43248512896257405886.recip(),
            PRECISION
        );
        assert_almost_eq!(r_rgamma(1.62123), 0.896081923385351.recip(), PRECISION);

        assert_almost_eq!(r_rgamma(-4.8), -0.062423361354759553.recip(), PRECISION);

        assert_almost_eq!(
            r_rgamma(1.0e-5 + 1.0e-16),
            99999.42279322556767360213300482199406241771308740302819426480.recip(),
            PRECISION
        );
        assert_almost_eq!(
            r_rgamma(0.1),
            9.513507698668731836292487177265402192550578626088377343050000.recip(),
            PRECISION
        );
        assert_almost_eq!(
            r_rgamma(1.0 - 1.0e-14),
            1.000000000000005772156649015427511664653698987042926067639529.recip(),
            PRECISION
        );
        assert_almost_eq!(
            r_rgamma(1.0 + 1.0e-14),
            0.99999999999999422784335098477029953441189552403615306268023.recip(),
            PRECISION
        );
        assert_almost_eq!(
            r_rgamma(PI / 2.0),
            0.890560890381539328010659635359121005933541962884758999762766.recip(),
            PRECISION
        );

        assert_almost_eq!(
            r_rgamma(5.0 - 1.0e-14),
            23.999999999999652.recip(),
            PRECISION
        );

        assert_almost_eq!(r_rgamma(10.1), 454760.7514415855.recip(), PRECISION);
        assert_almost_eq!(
            r_rgamma(150.0 + 1.0e-12),
            3.8089226376496421386707466577615064443807882167327097140e+260.recip(),
            PRECISION
        );
    }
}
