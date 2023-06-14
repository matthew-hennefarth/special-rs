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

use crate::special::gamma::{r_gamma, r_gammaln, r_gammasgn, r_poch, r_rgamma};

use num_complex::ComplexFloat;

/// Implementation of the Gamma and related functions for both real and complex-valued inputs.
///
/// For most implementations, a few properties are exploited to simplify the approximations. Firstly, the reflection property is used to always translate the value to the positive real axis.
/// $$
/// \Gamma(-z)\Gamma(z) = -\frac{\pi}{z\sin(\pi z)}
/// $$
/// Additionally, the recursive nature of the Gamma function is often used to move the value into some desired region to then approximate.
/// $$
/// \Gamma(z+1) = z\Gamma(z)
/// $$
pub trait Gamma: ComplexFloat {
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
    /// ## Notes on Real-Valued Implementation
    /// The implementation is loosely based on the [cephes implementation]. For large $x>33$, we utilize the Stirling series approximations which is given by
    /// $$
    /// \sqrt{\frac{2\pi}{x}} \left(\frac{x}{e}\right)^x \left(1 + \frac{1}{12 x} + \frac{1}{288 x^2} - \frac{139}{51840 x^3} - \frac{571}{2488320 x^4} + \ldots \right)
    /// $$
    /// For $x < 1.0\times 10^{-7}$, we utilize the Laurant series around $x=0$ as
    /// $$
    /// \Gamma(x) \approx \frac{1}{x} - \gamma + \left(\frac{\pi^2}{6} + \gamma^2\right)\frac{x}{2}
    /// $$
    /// where $\gamma$ is the [Euler-Mascheroni constant]. Otherwise, we move the value into the interval $(2,3)$ and use 2 rational functions of degree $6$ and $7$ to approximate the Gamma function in this interval.
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
    /// [Euler-Mascheroni constant]: crate::constants::f64::GAMMA
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
    /// Implementation is based on the [cephes implementation] in the Scipy package.
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
    /// The [gamma] function, for real-valued arguments $x$, is never zero and so this is a well-defined function. This is not well defined for complex-values arguments though.
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

    /// Reciprocal of the [gamma] function
    /// $$
    /// \frac{1}{\Gamma(z)}
    /// $$
    ///
    /// # Examples
    /// ```
    /// use sci_rs::special::Gamma;
    /// assert_eq!(1.5.rgamma(), 1.0/(1.5.gamma()));
    /// assert_eq!((-2.0).rgamma(), 0.0);
    /// assert_eq!((4.0).rgamma(), 1.0/(4.0.gamma())); // Should be 1/24
    /// ```
    ///
    /// ## Notes
    /// Since the [gamma] function is never zero, this is a well-defined function. Where $\Gamma(z)$ is undefined (negative integers and $0$), we return `0.0`. The implementation here is based off of the [cephes implementation] in the Scipy package.
    ///
    /// ## Implementation Details for Real-Valued Arguments
    /// Like the [cephes implementation], we use a Chebyshev series to order 16 to approximate values between $(0,1)$. For values outside this region, but $|x| < 34$, we use recursion to move the value into this interval. For $x > 34$, we use
    /// $$
    /// \frac{1}{\Gamma(x)} = e^{-\ln\Gamma(x)}
    /// $$
    /// Of course, overflow and underflows may occur for large enough values despite the function not having any singularities.
    ///
    /// [gamma]: crate::special::Gamma::gamma()
    /// [cephes implementation]: https://github.com/scipy/scipy/blob/46081a85c3a6ca4c45610f4207abf791985e17e0/scipy/special/cephes/rgamma.c
    fn rgamma(self) -> Self;

    /// Pochhammer symbol
    ///
    /// The Pochhammer symbol is defined as
    /// $$
    /// z^{(m)} = \frac{\Gamma(z + m)}{\Gamma(z)}
    /// $$
    /// This is a generalization of the rising factorial which, for non-negative integers $n$ and $m$ is
    /// $$
    /// n^{(m)} = n(n+1)(n+2)\ldots(n+m-1) = \prod_{k = 1}^{m} (n+k-1)
    /// $$
    /// See the [dlmf] or [wiki] page for more details.
    ///
    /// # Examples
    /// ```
    /// use sci_rs::special::Gamma;
    /// assert_eq!(2.0.poch(2.0), 4.0.gamma()/2.0.gamma());
    /// assert!(((-1.5_f32).poch(0.25) - (-1.25).gamma()/(-1.5).gamma()).abs() < 1e-6);
    /// ```
    ///
    /// ## Notes
    /// The real-value implementation here is based on that from the [cephes] library.
    ///
    /// [dlmf]: https://dlmf.nist.gov/5.2#iii
    /// [wiki]: https://en.wikipedia.org/wiki/Falling_and_rising_factorials
    /// [cephes]: https://github.com/scipy/scipy/blob/46081a85c3a6ca4c45610f4207abf791985e17e0/scipy/special/cephes/poch.c
    fn poch(self, m: Self) -> Self;
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
pub fn gammaln<T>(z: T) -> T
where
    T: Gamma,
{
    z.gammaln()
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
pub fn gammasgn<T>(z: T) -> T
where
    T: Gamma,
{
    z.gammasgn()
}

/// Reciprocal of the Gamma function.
///
/// Has the same semantics as [rgamma] in the [Gamma trait].
///
/// # Examples
/// ```
/// use sci_rs::special::{Gamma, rgamma};
/// assert_eq!(rgamma(1.5), 1.0/(1.5.gamma()));
/// assert_eq!(rgamma(-2.0), 0.0);
/// assert_eq!(rgamma(4.0), 4.0.rgamma()); // Should be 1/24
/// ```
///
/// [rgamma]: crate::special::Gamma::rgamma
/// [Gamma trait]: crate::special::Gamma
pub fn rgamma<T>(z: T) -> T
where
    T: Gamma,
{
    z.rgamma()
}

/// Pochhammer symbol
///
/// Has the same semantics as [poch] in the [Gamma trait].
///
/// # Examples
/// ```
/// use sci_rs::special::{Gamma, poch};
/// assert_eq!(poch(2.0, 3.0), 2.0.poch(3.0));
/// assert_eq!(poch(-1.5, -0.22), (-1.5).poch(-0.22));
/// ```
///
/// [poch]: crate::special::Gamma::poch
/// [Gamma trait]: crate::special::Gamma
pub fn poch<T>(z: T, m: T) -> T
where
    T: Gamma,
{
    z.poch(m)
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

            #[inline(always)]
            fn gammasgn(self) -> Self {
                r_gammasgn(self)
            }

            fn rgamma(self) -> Self {
                r_rgamma(self)
            }

            fn poch(self, m: Self) -> Self {
                r_poch(self, m)
            }
        }
    )*)
}

float_gamma_impl! {f32 f64}
