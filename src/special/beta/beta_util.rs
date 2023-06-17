//**********************************************************************
// This file is part of sci-rs                                         *
// Copyright 2023 Matthew R. Hennefarth                                *
//**********************************************************************

use crate::special::RealGamma;
use num_traits::{cast, Float, FloatConst};
use std::ops::{AddAssign, SubAssign};

pub(crate) trait RealBetaConsts {
    const ASYMP_FACTOR: Self;
    const MAX_GAMMA: Self;
}

macro_rules! float_rbetaconsts_impl {
    ($($T: ty)*) => ($(
        impl RealBetaConsts for $T {
            const ASYMP_FACTOR: Self = 1.0E6;
            const MAX_GAMMA: Self = 171.624376956302725;
        }
    )*)
}

float_rbetaconsts_impl! {f32 f64}

/// Asymptotic expansion for  ln(|B(a, b)|) for a > ASYMP_FACTOR*max(|b|, 1).
/// Taken from the Cephes library, and an unknown source.
pub(crate) fn r_lbeta_asymp<T>(a: T, b: T) -> T
where
    T: Float + RealGamma + FloatConst + SubAssign + AddAssign,
{
    let mut r = b.lgamma();

    let two = T::one() + T::one();
    let twelve = cast::<usize, T>(12).unwrap();

    r -= b * a.ln();
    r += b * (T::one() - b) / (two * a);
    r += b * (T::one() - b) * (T::one() - two * b) / (twelve * a * a);
    r -= b * b * (T::one() - b) * (T::one() - b) / (twelve * a * a * a);

    r
}
