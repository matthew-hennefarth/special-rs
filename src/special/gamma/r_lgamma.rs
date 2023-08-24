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

        const KNOWN_VALUES: [[f64; 2]; 17] = [
            // Evaluated at Rationals
            [1.0 / 3.0, 0.985420646927767069187174036977], // OEIS: A256165
            [0.25, 1.28802252469807745737061044021],       // OEIS: A256166
            [0.20, 1.52406382243078452488105649392],       // OEIS: A256167
            [1.0 / 6.0, 1.71673343507824046052784630958],  // OEIS: A255888
            [1.0 / 7.0, 1.87916927159583583645595640934],  // OEIS: A256609
            [1.0 / 8.0, 2.01941835755379634532029052116],  // OEIS: 255306
            [1.0 / 9.0, 2.14273180037669310488040788489],  // OEIS: A256610
            [0.1, 2.25271265173420595986970164636],        // OEIS: A256612
            [1.0 / 11.0, 2.35193461079879276046368095764], // OEIS: A256611
            [1.0 / 12.0, 2.44229731118288975091554935219], // OEIS: A256066
            // Other Important
            [1.0 / PI, 1.03364612576558270648553745533], // OEIS: A257957
            [12.5, 18.73434751193644570163],
            [13.34, 20.8506330413774776],
            [13.5, 21.2600761562447011], // Taken from Wolframalpha
            [14.5, 23.8627658416890849], // Taken from Wolframalpha
            [150.0 + 1.0e-12, 600.0094705553324354], // Taken from Wolframalpha
            [-1.72, 0.95458292505988099545337076566],
        ];

        for value in KNOWN_VALUES {
            assert_almost_eq!(r_lgamma(value[0]), value[1], PRECISION);
        }
    }
}
