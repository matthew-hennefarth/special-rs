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
use crate::special::Gamma;
use crate::traits::FloatSciConst;

/// Gamma related functions which only accept real-valued arguments.
pub trait RealGamma: Gamma {
    /// Natural logarithm of the absolute vale of the [gamma] function.
    /// $$
    /// \ln \left|\Gamma(x)\right|
    /// $$
    ///
    /// # Examples
    /// ```
    /// use sci_rs::special::RealGamma;
    /// assert_eq!(1.0.gammaln(), 0.0); // ln(1) = 0
    /// assert_eq!(2.0.gammaln(), 0.0); // ln(1) = 0
    /// assert!((14.5_f64.gammaln() - 23.86276584168908_f64).abs() < 1e-10 );
    /// ```
    ///
    /// ## Notes
    /// Implementation is taken from the [cephes implementation] in the Scipy package.
    ///
    /// # References
    /// - [cephes implementation]
    ///
    /// [gamma]: crate::special::gamma()
    /// [cephes implementation]: https://github.com/scipy/scipy/blob/main/scipy/special/cephes/gamma.c
    fn gammaln(self) -> Self;
}

/// Natural log of the Gamma function evaluated at $z$.
///
/// Has the same semantics as [gammaln] in the [Gamma trait].
///
/// # Examples
/// ```
/// use sci_rs::special::{RealGamma, gammaln};
/// assert_eq!(14.5_f64.gammaln(), gammaln(14.5_f64));
/// ```
/// [gammaln]: crate::special::RealGamma::gammaln
/// [Gamma trait]: crate::special::Gamma
pub fn gammaln<T>(x: T) -> T
where
    T: RealGamma,
{
    x.gammaln()
}

macro_rules! float_realgamma_impl {
    ($($T: ty)*) => ($(
        impl RealGamma for $T {
            fn gammaln(self) -> Self {
                fn gammaln_singularity() -> $T {
                    return <$T>::INFINITY;
                }
                if self < -34.0 {
                    // Utilize Euler reflection and compute gammaln at positive value.
                    let x_positive = -self;
                    let p = x_positive.floor();
                    if p == x_positive {
                        return gammaln_singularity();
                    }

                    let z = euler_reflection_prefactor(x_positive, p);
                    if z == 0.0 {
                        return gammaln_singularity();
                    }
                    return Self::LOG_PI() - z.ln() - x_positive.gammaln();
                }

                if self < 13.0 {
                    let mut z = 1.0;
                    let mut x = self;

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
                    const B: [$T; 6] = [
                        -1.37825152569120859100E3,
                        -3.88016315134637840924E4,
                        -3.31612992738871184744E5,
                        -1.16237097492762307383E6,
                        -1.72173700820839662146E6,
                        -8.53555664245765465627E5,
                    ];
                    const C: [$T; 7] = [
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

                let mut q = (self - 0.5) * self.ln() - self + Self::LOG_SQRT_2_PI();
                if (self > 1.0e8) {
                    return q;
                }
                let p = 1.0/(self * self);
                if self >= 1.0e3 {
                    q += ((7.9365079365079365079365e-4 * p
                        - 2.7777777777777777777778e-3) * p
	                    + 0.0833333333333333333333) / self;
                } else {
                    const A: [$T; 5] = [
                        8.11614167470508450300E-4,
                        -5.95061904284301438324E-4,
                        7.93650340457716943945E-4,
                        -2.77777777730099687205E-3,
                        8.33333333333331927722E-2
                    ];
                    q += eval_poly(p, &A) / self;
                }

                q
            }
        }
    )*)
}

float_realgamma_impl! {f32 f64}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::special::Factorial;

    const PRECISION: f64 = 1E-14;

    #[test]
    fn test_gammaln() {
        for i in 1..10 {
            assert_almost_eq!(
                gammaln(i as f64),
                ((i - 1).factorial() as f64).ln(),
                PRECISION
            );
            assert_eq!(gammaln(-i as f64), f64::INFINITY);
        }

        // Evaluated at Rationals
        assert_almost_eq!(
            gammaln(1.0 / 3.0),
            0.985420646927767069187174036977,
            PRECISION
        ); // OEIS: A256165
        assert_almost_eq!(gammaln(0.25), 1.28802252469807745737061044021, PRECISION); // OEIS: A256166
        assert_almost_eq!(gammaln(0.20), 1.52406382243078452488105649392, PRECISION); // OEIS: A256167
        assert_almost_eq!(
            gammaln(1.0 / 6.0),
            1.71673343507824046052784630958,
            PRECISION
        ); // OEIS: A255888
        assert_almost_eq!(
            gammaln(1.0 / 7.0),
            1.87916927159583583645595640934,
            PRECISION
        ); // OEIS: A256609
        assert_almost_eq!(
            gammaln(1.0 / 8.0),
            2.01941835755379634532029052116,
            PRECISION
        ); // OEIS: 255306
        assert_almost_eq!(
            gammaln(1.0 / 9.0),
            2.14273180037669310488040788489,
            PRECISION
        ); // OEIS: A256610
        assert_almost_eq!(gammaln(0.1), 2.25271265173420595986970164636, PRECISION); // OEIS: A256612
        assert_almost_eq!(
            gammaln(1.0 / 11.0),
            2.35193461079879276046368095764,
            PRECISION
        ); // OEIS: A256611
        assert_almost_eq!(
            gammaln(1.0 / 12.0),
            2.44229731118288975091554935219,
            PRECISION
        ); // OEIS: A256066

        // Other important
        assert_almost_eq!(
            gammaln(f64::PI().recip()),
            1.03364612576558270648553745533,
            PRECISION
        ); // OEIS: A257957

        assert_almost_eq!(gammaln(12.5), 18.73434751193644570163, PRECISION);
        assert_almost_eq!(gammaln(13.34), 20.8506330413774776, PRECISION);
        assert_almost_eq!(gammaln(-34.5), -89.2102003799689880, 1e-13); // Taken from Wolframalpha
        assert_almost_eq!(gammaln(-40.2), -109.3852746800507908, 1e-13); // Taken from Wolframalpha
        assert_almost_eq!(gammaln(13.5), 21.2600761562447011, PRECISION); // Taken from Wolframalpha
        assert_almost_eq!(gammaln(14.5), 23.8627658416890849, PRECISION); // Taken from Wolframalpha
        assert_almost_eq!(gammaln(150.0 + 1.0e-12), 600.0094705553324354, PRECISION);
        // Taken from Wolframalpha
    }
}
