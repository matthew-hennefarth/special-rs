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

use crate::special::gamma::{r_gammaln, r_gammasgn, RealGammaLnConsts};
use crate::traits::FloatSciConst;
use num_traits::{cast, Float};
use std::ops::{AddAssign, DivAssign, MulAssign, SubAssign};

fn is_nonpositive_int<T>(x: T) -> bool
where
    T: Float,
{
    x.is_finite() && x <= T::zero() && x == x.ceil()
}

pub(crate) trait RealPochConsts {
    const MIN_FOR_EXP: Self;
}

macro_rules! impl_realpochconsts {
    ($($T: ty)*) => ($(
        impl RealPochConsts for $T {
            const MIN_FOR_EXP: Self = 1.0e4;
        }
)*)
}

impl_realpochconsts! {f32 f64}

/// Implementation of the Pochhammer symbols for real-valued arguments
/// $$
/// x^{(m)} = \frac{\Gamma(x+m)}{\Gamma(x)}
/// $$
/// where both $x$ and $m$ are real-valued.
pub(crate) fn r_poch<T>(x: T, mut m: T) -> T
where
    T: Float
        + SubAssign
        + AddAssign
        + MulAssign
        + DivAssign
        + RealPochConsts
        + FloatSciConst
        + RealGammaLnConsts,
{
    let mut r = T::one();

    while m >= T::one() {
        if x + m == T::one() {
            break;
        }
        m -= T::one();
        r *= x + m;
        if !r.is_finite() || r.is_zero() {
            break;
        }
    }

    while m <= -T::one() {
        if x + m == T::one() {
            break;
        }
        r /= x + m;
        m += T::one();
        if !r.is_finite() || r.is_zero() {
            break;
        }
    }

    if m.is_zero() {
        return r;
    }
    if x > T::MIN_FOR_EXP && m.abs() <= T::one() {
        let two = T::one() + T::one();
        let three = two + T::one();
        return r
            * x.powf(m)
            * (T::one()
                + m * (m - T::one()) / (two * x)
                + m * (m - T::one()) * (m - two) * (three * m - T::one())
                    / (cast::<u8, T>(24).unwrap() * x * x)
                + m * m * (m - T::one()) * (m - T::one()) * (m - two) * (m - three)
                    / (cast::<u8, T>(48).unwrap() * x * x * x));
    }

    // Check for infinite
    if is_nonpositive_int(x + m) && !is_nonpositive_int(x) && x + m != m {
        return T::nan();
    }

    // Check for zero
    if !is_nonpositive_int(x + m) && is_nonpositive_int(x) {
        return T::zero();
    }

    r * (r_gammaln(x + m) - r_gammaln(x)).exp() * r_gammasgn(x + m) * r_gammasgn(x)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::constants::f64::PI;
    use crate::special::Factorial;

    const PRECISION: f64 = 1e-14;

    #[test]
    fn test_rpoch() {
        for i in 1..10 {
            for m in 0..5 {
                assert_eq!(
                    r_poch(i as f64, m as f64),
                    (i + m - 1).factorial() as f64 / (i - 1).factorial() as f64
                );
            }
        }

        // Weird edge values
        assert_eq!(r_poch(0.0, 0.0), 1.0);
        assert_eq!(r_poch(0.0, 0.25), 0.0);
        assert_eq!(r_poch(-1.0, 0.0), 1.0);
        assert_eq!(r_poch(-1.0, 0.25), 0.0);
        // This feels wrong since this is gamma(-1)/gamma(-2) which should be undefined. But technically you can use the rules of the gamma function to write as -2 gamma(-2)/gamma(-2). So I am not sure how I feel about it. But it agrees with SciPy currently.
        assert_eq!(r_poch(-2.0, 1.0), -2.0);
        assert!(r_poch(-2.2, 0.2).is_nan());

        // Following values taken from SciPy version 1.10.1
        assert_almost_eq!(r_poch(0.5, 0.5), 0.56418958354775639030, PRECISION);
        assert_almost_eq!(r_poch(0.5, 1.0), 0.5, PRECISION);
        assert_almost_eq!(r_poch(1.0, 0.5), 0.88622692545275794096, PRECISION);
        assert_almost_eq!(r_poch(1.5, 1.0), 1.5, PRECISION);
        assert_almost_eq!(r_poch(1.0, 1.5), 1.32934038817913702246, PRECISION);
        assert_almost_eq!(r_poch(1.5, 1.5), 2.25675833419102511712, PRECISION);
        assert_almost_eq!(r_poch(2.5, 2.5), 18.05406667352819738426, PRECISION);
        assert_almost_eq!(
            r_poch(150.00001, PI),
            7015772.59900571219623088837,
            PRECISION
        );

        // Add some larger values to compare to!
        // From Wolframalpha
        assert_almost_eq!(r_poch(PI / 2.0, PI), 17.63940522158362397144, PRECISION);
        assert_almost_eq!(
            r_poch(-34.54, -2.4),
            0.000959271152790991576995792318463,
            PRECISION
        );

        assert_almost_eq!(r_poch(-1.5, -0.22), 1.099148632503722270901806, PRECISION);
    }
}
