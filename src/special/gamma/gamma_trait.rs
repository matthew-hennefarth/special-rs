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

use crate::special::gamma::{r_gamma, r_gammaln, r_gammasgn};

/// Implementation of the Gamma and related functions for both real and complex-valued inputs.
pub trait Gamma {
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

    /// Natural logarithm of the absolute vale of the [gamma] function.
    ///
    /// $$
    /// \ln \left|\Gamma(x)\right|
    /// $$
    /// Useful since it avoids the need of choosing a particular branch of the complex log function. For real-valued arguments $x$, we have the additional relationship:
    /// $$
    /// \Gamma(x) = \text{gammasgn}(x)e^{\ln{\left|\Gamma(x)\right|}}
    /// $$
    /// where [gammasgn] is the sign of the [gamma] function
    ///
    /// # Examples
    /// ```
    /// use sci_rs::special::Gamma;
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
    /// [gamma]: crate::special::Gamma::gamma()
    /// [gammasgn]: crate::special::Gamma::gammasgn()
    /// [cephes implementation]: https://github.com/scipy/scipy/blob/main/scipy/special/cephes/gamma.c
    fn gammaln(self) -> Self;

    /// Sign of the [gamma] function.
    ///
    /// $$
    /// \text{gammasgn}(x) = \begin{cases}
    /// +1.0 & \Gamma(x) > 0 \\\\
    /// -1.0 & \Gamma(x) < 0
    /// \end{cases}
    /// $$
    /// The [gamma] function is never zero and so this is a well-defined function on the gamma function domain.
    ///
    /// # Examples
    /// ```
    /// use sci_rs::special::Gamma;
    /// assert_eq!(1.23.gammasgn(), 1.0);
    /// assert_eq!((-0.23).gammasgn(), -1.0);
    /// assert_eq!((-1.5).gammasgn(), 1.0);
    /// ```
    ///
    /// ## Notes
    /// We return $0.0$ if $\Gamma(x)$  is undefined (where [gamma] returns `NaN` or `Inf`). This is $x=0.0, -1, -2, \ldots$.
    ///
    /// [gamma]: crate::special::Gamma::gamma()
    fn gammasgn(self) -> Self;
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

/// Natural log of the Gamma function evaluated at $z$.
///
/// Has the same semantics as [gammaln] in the [Gamma trait].
///
/// # Examples
/// ```
/// use sci_rs::special::{Gamma, gammaln};
/// assert_eq!(14.5_f64.gammaln(), gammaln(14.5_f64));
/// ```
/// [gammaln]: crate::special::Gamma::gammaln
/// [Gamma trait]: crate::special::Gamma
pub fn gammaln<T>(x: T) -> T
where
    T: Gamma,
{
    x.gammaln()
}

/// Sign of the Gamma function.
///
/// Has the same semantics as [gammasgn] in the [Gamma trait]
///
/// # Examples
/// ```
/// use sci_rs::special::{Gamma, gammasgn};
/// assert_eq!(5.2.gammasgn(), gammasgn(5.2));
/// ```
///
/// [gammasgn]: crate::special::Gamma::gammasgn
/// [Gamma trait]: crate::special::Gamma
pub fn gammasgn<T>(x: T) -> T
where
    T: Gamma,
{
    x.gammasgn()
}

macro_rules! float_gamma_impl {
    ($($T: ty)*) => ($(
        impl Gamma for $T {
            #[inline(always)]
            fn gamma(self) -> Self {
                r_gamma(self)
            }

            #[inline(always)]
            fn gammaln(self) -> Self {
                r_gammaln(self)
            }

            fn gammasgn(self) -> Self {
                r_gammasgn(self)
            }
        }
    )*)
}

float_gamma_impl! {f32 f64}
