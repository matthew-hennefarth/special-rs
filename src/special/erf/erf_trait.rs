//**********************************************************************
// This file is part of sci-rs                                         *
// Copyright 2023 Matthew R. Hennefarth                                *
//**********************************************************************

use crate::special::erf::real_erf_impl::*;

/// Error and related functions.
pub trait Erf {
    /// The Error Function.
    /// $$
    /// \mathrm{erf} z = \frac{2}{\pi}\int_0^z e^{-t^2}dt
    /// $$
    /// It is an entire complex function with no singularities (except at infinity). It is also called the Gauss error function. For a more thorough explanation of the error function and all of its properties, it is recommended to read through the [DLMF] and [wiki] pages.
    ///
    /// # Examples
    /// For real-valued inputs
    /// ```
    /// use sci_rs::special::Erf;
    /// assert!((2.0_f32.erf() - 0.9953222650189527341).abs() < 1e-14);
    /// assert_eq!((-3.2).erf(), -(3.2.erf())); // Error Function is odd!
    /// ```
    /// # Notes
    /// The real-valued implementation is based on the [Boost] library version (v 1.82.0) ([documentation](https://www.boost.org/doc/libs/1_82_0/libs/math/doc/html/math_toolkit/sf_erf/error_function.html)). Largely, it uses various rational polynomial approximations. Additionally, it uses the fact that the Error function is odd to always evaluate at positive arguments. At times, it will use the [Compliment Error function] to evaluate.
    ///
    /// # References
    /// - [Boost Error function implementation][Boost] version 1.82.0
    /// - [DLMF]
    /// - [Wikipedia][wiki]
    ///
    /// [comment]: <> (Reference hyperlinks)
    /// [Boost]: https://www.boost.org/doc/libs/1_82_0/boost/math/special_functions/erf.hpp
    /// [DLMF]: https://dlmf.nist.gov/7.2
    /// [wiki]: https://en.wikipedia.org/wiki/Error_function
    /// [Compliment Error function]: crate::special::erf::Erf::erfc
    fn erf(self) -> Self;

    /// The Complementary Error Function.
    /// $$
    /// \mathrm{erfc}z = 1 - \mathrm{erf}z
    /// $$
    /// See the [wiki] and [DLMF] pages for a more detailed explanation of the [Error function] and its properties.
    ///
    /// # Examples
    /// For real-valued inputs
    /// ```
    /// use sci_rs::special::Erf;
    /// assert!((2.0_f32.erfc() - 0.004677734981047265837).abs() < 1e-7);
    /// assert_eq!((-3.2).erfc(), 1.0 + (3.2.erf())); // Error and Complementary are related
    /// assert_eq!(3.0.erfc() + 3.0.erf(), 1.0);
    /// ```
    /// # Notes
    /// The real-valued implementation is based on the [Boost] library (v 1.82.0) ([documentation](https://www.boost.org/doc/libs/1_82_0/libs/math/doc/html/math_toolkit/sf_erf/error_function.html)).
    ///
    /// [comment]: <> (Reference hyperlinks)
    /// [Boost]: https://www.boost.org/doc/libs/1_82_0/boost/math/special_functions/erf.hpp
    /// [wiki]: https://en.wikipedia.org/wiki/Error_function
    /// [DLMF]: https://dlmf.nist.gov/7.2
    /// [Error function]: crate::special::erf::Erf::erf
    fn erfc(self) -> Self;
}

macro_rules! float_erf_impl {
    ($($T: ty)*) => ($(
        impl Erf for $T {
            #[inline(always)]
            fn erf(self) -> Self {
                r_erf(self, false)
            }

            #[inline(always)]
            fn erfc(self) -> Self {
                r_erf(self, true)
            }
        }
    )*)
}

float_erf_impl! {f32 f64}
