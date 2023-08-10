//**********************************************************************
// This file is part of sci-rs                                         *
// Copyright 2023 Matthew R. Hennefarth                                *
//**********************************************************************

use crate::special::gamma::gamma_util::{
    euler_reflection_prefactor, is_gamma_pole, lngamma_stirling, LnGammaStirlingConsts,
};
use crate::special::gamma::r_gamma::RealGammaConsts;
use crate::special::tools::eval_poly;
use crate::traits::FloatSciConst;
use num_traits::Float;
use std::ops::{AddAssign, DivAssign, MulAssign, SubAssign};

#[inline(always)]
fn gammaln_singularity<T>() -> T
where
    T: Float,
{
    T::infinity()
}

pub(crate) trait RealGammaLnConsts:
    Sized + RealGammaConsts + Float + LnGammaStirlingConsts
{
    const B: [Self; 6];
    const C: [Self; 7];
    const MAX_TO_RECURSE: Self;
}

macro_rules! impl_realgammalnconsts_coefficients {
    ($($T: ty)*) => ($(
        impl RealGammaLnConsts for $T {
            const B: [Self; 6] = [
                -1.37825152569120859100E3,
                -3.88016315134637840924E4,
                -3.31612992738871184744E5,
                -1.16237097492762307383E6,
                -1.72173700820839662146E6,
                -8.53555664245765465627E5,
            ];
            const C: [Self; 7] = [
                1.00000000000000000000E0,
                -3.51815701436523470549E2,
                -1.70642106651881159223E4,
                -2.20528590553854454839E5,
                -1.13933444367982507207E6,
                -2.53252307177582951285E6,
                -2.01889141433532773231E6,
            ];
            const MAX_TO_RECURSE: Self = 13.0;
        }
)*)
}

impl_realgammalnconsts_coefficients! {f32 f64}

/// Implementation of the `lgamma()` for real-valued inputs.
/// $$
/// \ln\left|\Gamma(x)right|
/// $$
/// where $x$ is real-valued.
pub(crate) fn r_lgamma<T>(x: T) -> T
where
    T: Float + FloatSciConst + SubAssign + MulAssign + DivAssign + AddAssign + RealGammaLnConsts,
{
    if is_gamma_pole(x) {
        return gammaln_singularity();
    }

    if x.is_sign_negative() {
        // Utilize Euler reflection and compute gammaln at positive value.
        return euler_reflection_prefactor(x).abs().ln() - r_lgamma(-x);
    }

    if x >= T::MAX_TO_RECURSE {
        return lngamma_stirling(x);
    }

    // Get into the range of (2, 3]
    let mut z = T::one();
    let mut x = x;

    while x >= T::INTERVAL[1] {
        x -= T::one();
        z *= x;
    }

    while x < T::INTERVAL[0] {
        z /= x;
        x += T::one();
    }

    // Gamma(2) = 1! = 1
    if x == T::INTERVAL[0] {
        return z.ln();
    }

    x -= T::INTERVAL[0];
    z.ln() + x * eval_poly(x, &T::B) / eval_poly(x, &T::C)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::constants::f64::PI;
    use crate::special::Factorial;

    const PRECISION: f64 = 1e-14;

    #[test]
    fn test_gammaln() {
        for i in 1..10 {
            assert_almost_eq!(
                r_lgamma(i as f64),
                ((i - 1).factorial() as f64).ln(),
                PRECISION
            );
            assert_eq!(r_lgamma(-i as f64), f64::INFINITY);
        }

        // Evaluated at Rationals
        assert_almost_eq!(
            r_lgamma(1.0 / 3.0),
            0.985420646927767069187174036977,
            PRECISION
        ); // OEIS: A256165
        assert_almost_eq!(r_lgamma(0.25), 1.28802252469807745737061044021, PRECISION); // OEIS: A256166
        assert_almost_eq!(r_lgamma(0.20), 1.52406382243078452488105649392, PRECISION); // OEIS: A256167
        assert_almost_eq!(
            r_lgamma(1.0 / 6.0),
            1.71673343507824046052784630958,
            PRECISION
        ); // OEIS: A255888
        assert_almost_eq!(
            r_lgamma(1.0 / 7.0),
            1.87916927159583583645595640934,
            PRECISION
        ); // OEIS: A256609
        assert_almost_eq!(
            r_lgamma(1.0 / 8.0),
            2.01941835755379634532029052116,
            PRECISION
        ); // OEIS: 255306
        assert_almost_eq!(
            r_lgamma(1.0 / 9.0),
            2.14273180037669310488040788489,
            PRECISION
        ); // OEIS: A256610
        assert_almost_eq!(r_lgamma(0.1), 2.25271265173420595986970164636, PRECISION); // OEIS: A256612
        assert_almost_eq!(
            r_lgamma(1.0 / 11.0),
            2.35193461079879276046368095764,
            PRECISION
        ); // OEIS: A256611
        assert_almost_eq!(
            r_lgamma(1.0 / 12.0),
            2.44229731118288975091554935219,
            PRECISION
        ); // OEIS: A256066

        // Other important
        assert_almost_eq!(
            r_lgamma(PI.recip()),
            1.03364612576558270648553745533,
            PRECISION
        ); // OEIS: A257957

        assert_almost_eq!(r_lgamma(12.5), 18.73434751193644570163, PRECISION);
        assert_almost_eq!(r_lgamma(13.34), 20.8506330413774776, PRECISION);
        assert_almost_eq!(r_lgamma(-34.5), -89.2102003799689880, 1e-13); // Taken from Wolframalpha
        assert_almost_eq!(r_lgamma(-40.2), -109.3852746800507908, 1e-13); // Taken from Wolframalpha
        assert_almost_eq!(r_lgamma(13.5), 21.2600761562447011, PRECISION); // Taken from Wolframalpha
        assert_almost_eq!(r_lgamma(14.5), 23.8627658416890849, PRECISION); // Taken from Wolframalpha
        assert_almost_eq!(r_lgamma(150.0 + 1.0e-12), 600.0094705553324354, PRECISION);
        // Taken from Wolframalpha

        assert_almost_eq!(r_lgamma(-1.72), 0.95458292505988099545337076566, PRECISION);
    }
}
