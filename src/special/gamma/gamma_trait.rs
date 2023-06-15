//**********************************************************************
// This file is part of sci-rs                                         *
// Copyright 2023 Matthew R. Hennefarth                                *
//**********************************************************************
use crate::special::gamma::complex_gamma_impl::*;
use crate::special::gamma::real_gamma_impl::*;

use num_complex::Complex;

/// Gamma and related functions for both real and complex-valued valued arguments.
///
/// # Implementation Notes
/// For most implementations, a few properties are exploited to simplify the approximations. Firstly, the reflection property is used to always translate the value to the positive real axis.
/// $$
/// \Gamma(-z)\Gamma(z) = -\frac{\pi}{z\sin(\pi z)}
/// $$
/// Additionally, the recursive nature of the Gamma function is often used to move the value into some desired region to then approximate.
/// $$
/// \Gamma(z+1) = z\Gamma(z)
/// $$
pub trait Gamma {
    /// The Gamma function.
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
    /// For real-valued inputs:
    /// ```
    /// use sci_rs::special::Gamma;
    /// assert_eq!(4.0_f32.gamma(), 6.0); // Gamma(4) = 3!
    /// assert!((0.0_f64).gamma().is_nan()); // Gamma(0) is undefined
    /// assert!(((4.5_f64).gamma() - 11.6317283).abs() <  1e-5);
    /// ```
    /// Also for complex-valued inputs:
    /// ```
    /// use sci_rs::special::Gamma;
    /// use num_complex::Complex32;
    /// let z = Complex32{re: 1.2, im: -0.4};
    /// println!("{}", z.gamma()); // 0.828 + 0.0836 j
    /// ```
    /// # Notes
    /// The real-valued implementation is loosely based on the [cephes implementation]. For large $x>33$, we utilize the Stirling series approximations which is given by
    /// $$
    /// \sqrt{\frac{2\pi}{x}} \left(\frac{x}{e}\right)^x \left(1 + \frac{1}{12 x} + \frac{1}{288 x^2} - \frac{139}{51840 x^3} - \frac{571}{2488320 x^4} + \ldots \right)
    /// $$
    /// For $x < 1.0\times 10^{-7}$, we utilize the Laurant series around $x=0$ as
    /// $$
    /// \Gamma(x) \approx \frac{1}{x} - \gamma + \left(\frac{\pi^2}{6} + \gamma^2\right)\frac{x}{2}
    /// $$
    /// where $\gamma$ is the [Euler-Mascheroni constant]. Otherwise, we move the value into the interval $(2,3)$ and use 2 rational functions of degree $6$ and $7$ to approximate the Gamma function in this interval.
    ///
    /// For complex inputs, we simply use
    /// $$
    /// \Gamma(z) = e^{\ln(\Gamma(z)}
    /// $$
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
    /// assert_eq!(1.0.lgamma(), 0.0); // ln(1) = 0
    /// assert_eq!(2.0.lgamma(), 0.0); // ln(1) = 0
    /// assert!((14.5_f64.lgamma() - 23.86276584168908_f64).abs() < 1e-10 );
    /// ```
    ///
    /// # Notes
    /// Implementation is based on the [cephes implementation] in the Scipy package. Note however though that the SciPy package does not implement this function for complex-valued arguments. Here we have implemented it simply as
    /// $$
    /// \ln\left|\Gamma(z)\right| = \Re \left(\ln\Gamma(z)\right)8
    /// $$
    ///
    /// # References
    /// - [cephes implementation]
    ///
    /// [gamma]: crate::special::Gamma::gamma()
    /// [gammasgn]: crate::special::RealGamma::gammasgn()
    /// [cephes implementation]: https://github.com/scipy/scipy/blob/main/scipy/special/cephes/gamma.c
    fn lgamma(self) -> Self;

    /// Principle branch of the natural log of the Gamma function.
    /// $$
    /// \ln\Gamma(z)
    /// $$
    /// The function has a single branch cut on the negative real axis.
    ///
    /// # Examples
    /// Real-valued inputs must be greater than $0.0$.
    /// ```
    /// use sci_rs::special::Gamma;
    /// assert_eq!(1.0.lngamma(), 0.0);
    /// assert_eq!(1.5.lngamma(), 1.5.lgamma());
    /// assert!((-1.5_f32).lngamma().is_nan()); // Result would be imaginary
    /// ```
    /// To get a number, they must be explicitly cast to a complex number.
    /// ```
    /// use sci_rs::special::Gamma;
    /// use num_complex::Complex32;
    /// // We have to explicitly cast negative real numbers to complex-values.
    /// println!("{}", Complex32{re: -1.5, im: 0.0}.lngamma()); // 0.8600 - 6.283j
    /// println!("{}", Complex32{re: 0.0, im: 1.0}.lngamma()); // -0.6509 - 1.872j
    /// ```
    ///
    /// # Notes
    /// Implementation is based on that from the SciPy (v1.10.1) package. For real-valued arguments, returns `NaN` when $x \leq 0.0$ since the log of a negative number is complex-valued.
    fn lngamma(self) -> Self;

    /// Reciprocal of the [Gamma] function.
    /// $$
    /// \frac{1}{\Gamma(z)}
    /// $$
    ///
    /// # Examples
    /// For real-valued arguments:
    /// ```
    /// use sci_rs::special::Gamma;
    /// assert_eq!(1.5.rgamma(), 1.0/(1.5.gamma()));
    /// assert_eq!((-2.0).rgamma(), 0.0);
    /// assert_eq!((4.0).rgamma(), 1.0/(4.0.gamma())); // Should be 1/24
    /// ```
    /// For complex-valued arguments:
    /// ```
    /// use sci_rs::special::Gamma;
    /// use num_complex::Complex32;
    /// let z = Complex32{re: 1.0, im: 1.0};
    /// assert_eq!(z.rgamma(), 1.0/z.gamma());
    /// ```
    ///
    /// # Notes
    /// Since the [Gamma] function is never zero, this is a well-defined function. Where $\Gamma(z)$ is undefined (negative integers and $0$), we return `0.0`. The implementation here is based off of the [cephes implementation] in the Scipy package.
    ///
    /// ## Implementation Details for Real-Valued Arguments
    /// Like the [cephes implementation], we use a Chebyshev series to order 16 to approximate values between $(0,1)$. For values outside this region, but $|x| < 34$, we use recursion to move the value into this interval. For $x > 34$, we use
    /// $$
    /// \frac{1}{\Gamma(x)} = e^{-\ln\Gamma(x)}
    /// $$
    /// Of course, overflow and underflows may occur for large enough values despite the function not having any singularities.
    ///
    /// [Gamma]: crate::special::Gamma::gamma()
    /// [cephes implementation]: https://github.com/scipy/scipy/blob/46081a85c3a6ca4c45610f4207abf791985e17e0/scipy/special/cephes/rgamma.c
    fn rgamma(self) -> Self;
}

/// Gamma related functions which only make sense, or are only currently supported for real-valued arguments
pub trait RealGamma: Gamma {
    /// Sign of the [gamma] function.
    ///
    /// $$
    /// \text{gammasgn}(x) = \begin{cases}
    /// +1.0 & \Gamma(x) > 0 \\\\
    /// -1.0 & \Gamma(x) < 0
    /// \end{cases}
    /// $$
    /// The [gamma] function, for real-valued arguments $x$, is never zero and so this is a well-defined function.
    ///
    /// # Examples
    /// ```
    /// use sci_rs::special::RealGamma;
    /// assert_eq!(1.23.gammasgn(), 1.0);
    /// assert_eq!((-0.23).gammasgn(), -1.0);
    /// assert_eq!((-1.5).gammasgn(), 1.0);
    /// ```
    /// # Notes
    /// We return $0.0$ if $\Gamma(x)$ is undefined (where [gamma] returns `NaN` or `Inf`). This is $x=0.0, -1, -2, \ldots$.
    ///
    /// [gamma]: crate::special::Gamma::gamma()
    fn gammasgn(self) -> Self;

    /// Pochhammer symbol.
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
    /// use sci_rs::special::{RealGamma, Gamma};
    /// assert_eq!(2.0.poch(2.0), 4.0.gamma()/2.0.gamma());
    /// assert!(((-1.5_f32).poch(0.25) - (-1.25).gamma()/(-1.5).gamma()).abs() < 1e-6);
    /// ```
    ///
    /// # Notes
    /// The real-value implementation here is based on that from the [cephes] library.
    ///
    /// [dlmf]: https://dlmf.nist.gov/5.2#iii
    /// [wiki]: https://en.wikipedia.org/wiki/Falling_and_rising_factorials
    /// [cephes]: https://github.com/scipy/scipy/blob/46081a85c3a6ca4c45610f4207abf791985e17e0/scipy/special/cephes/poch.c
    fn poch(self, m: Self) -> Self;
}

macro_rules! float_gamma_impl {
    ($($T: ty)*) => ($(
        impl Gamma for $T {
            #[inline(always)]
            fn gamma(self) -> Self {
                r_gamma(self)
            }

            #[inline(always)]
            fn lgamma(self) -> Self {
                r_lgamma(self)
            }

            #[inline(always)]
            fn lngamma(self) -> Self {
                if self < 0.0 {
                    return Self::NAN;
                }
                self.lgamma()
            }

            #[inline(always)]
            fn rgamma(self) -> Self {
                r_rgamma(self)
            }
        }

        impl RealGamma for $T {
            #[inline(always)]
            fn gammasgn(self) -> Self {
                r_gammasgn(self)
            }

            #[inline(always)]
            fn poch(self, m: Self) -> Self {
                r_poch(self, m)
            }
        }
    )*)
}

float_gamma_impl! {f32 f64}

macro_rules! float_complexgamma_impl {
    ($($T: ty)*) => ($(
        impl Gamma for Complex<$T> {
            #[inline(always)]
            fn gamma(self) -> Self {
                c_gamma(self)
            }

            #[inline(always)]
            fn lgamma(self) -> Self {
                c_lgamma(self).into()
            }

            #[inline(always)]
            fn lngamma(self) -> Self {
                c_lngamma(self)
            }

            #[inline(always)]
            fn rgamma(self) -> Self {
                c_rgamma(self)
            }
        }
    )*)
}

float_complexgamma_impl! {f32 f64}
