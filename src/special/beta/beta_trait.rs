//**********************************************************************
// This file is part of sci-rs                                         *
// Copyright 2023 Matthew R. Hennefarth                                *
//**********************************************************************

use crate::special::beta::real_beta_impl::*;

/// Beta and related functions for real-valued arguments.
pub trait Beta {
    /// The Beta function.
    /// $$
    /// B(a,b) = \int^{1}_{0} t^{a-1}(1-t)^{b-1}dt = \frac{\Gamma(a)\Gamma(b)}{\Gamma(a+b)}
    /// $$
    /// where $\Gamma$ is the [Gamma] function. For a more detailed discussion on the Beta function, see the [dlmf] or [wiki] page.
    ///
    /// # Examples
    /// ```
    /// use sci_rs::special::{Beta, Gamma};
    /// assert_eq!(2.0.beta(2.0), 1.0/6.0); // 1! * 1! / 3! = 1/6
    /// assert_eq!(3.2.beta(-1.2), 3.2.gamma() * (-1.2).gamma() / (3.2-1.2).gamma());
    /// ```
    /// Note, when negative integers are included, one cannot simply use the Gamma functions directly! The Beta function exploits the Gamma function to remove computing the Gamma function at undefined points.
    /// ```
    /// use sci_rs::special::{Beta, Gamma};
    /// // Gamma function is undefined for negative integers!
    /// assert!(((-2.0_f32).gamma() * 1.0.gamma() / (-1.0).gamma()).is_nan());
    /// // Still possible to compute the Beta function though
    /// assert_eq!((-2.0).beta(1.0), -0.5);
    /// ```
    ///
    /// # Notes
    /// For real-valued inputs, the implementation is based on that from the [cephes] library in SciPy (v 1.10.1). For small enough inputs, the normal [Gamma] function is used, though some care is made to prevent overflows. For large values, an asymptotic expansion for $\ln\left|B(a,b)\right|$ is used.
    ///
    /// [Gamma]: crate::special::Gamma::gamma
    /// [dlmf]: https://dlmf.nist.gov/5.12
    /// [wiki]: https://en.wikipedia.org/wiki/Beta_function
    /// [cephes]: https://github.com/scipy/scipy/blob/main/scipy/special/cephes/beta.c
    fn beta(self, b: Self) -> Self;
}

macro_rules! float_beta_impl {
    ($($T: ty)*) => ($(
        impl Beta for $T {
            fn beta(self, b: Self) -> Self {
                r_beta(self, b)
            }
        }
    )*)
}

float_beta_impl! {f32 f64}
