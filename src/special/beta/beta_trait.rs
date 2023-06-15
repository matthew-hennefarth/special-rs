//**********************************************************************
// This file is part of sci-rs                                         *
// Copyright 2023 Matthew R. Hennefarth                                *
//**********************************************************************

use crate::special::beta::real_beta_impl::*;

/// Beta and related functions for real-valued arguments.
pub trait Beta {
    /// Beta function.
    /// $$
    /// B(a,b) = \int^{1}_{0} t^{a-1}(1-t)^{b-1}dt = \frac{\Gamma(a)\Gamma(b)}{\Gamma(a+b)}
    /// $$
    /// where $\Gamma$ is the [Gamma] function. For a more detailed discussion on the Beta function, see the [dlmf] or [wiki] page.
    ///
    /// # Notes
    /// For real-valued inputs, the implementation is based on that from the [cephes] library in SciPy (v 1.10.1).
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
