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
    /// [comment]: <> (Reference hyperlinks)
    /// [DLMF]: https://dlmf.nist.gov/7.2
    /// [wiki]: https://en.wikipedia.org/wiki/Error_function
    fn erf(self) -> Self;

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
