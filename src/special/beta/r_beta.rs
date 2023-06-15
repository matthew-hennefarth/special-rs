//**********************************************************************
// This file is part of sci-rs                                         *
// Copyright 2023 Matthew R. Hennefarth                                *
//**********************************************************************

use crate::special::gamma_util::is_gamma_pole;
use std::mem;
use std::ops::{AddAssign, SubAssign};

use crate::special::RealGamma;
use num_traits::{cast, Float, FloatConst};

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

/// For large arguments, the log of the function is evaluated using lgam() then exponentiated.
///
/// Implementation comes from the cephes library in SciPy: https://github.com/scipy/scipy/blob/main/scipy/special/cephes/beta.c#L55
pub(crate) fn r_beta<T>(mut a: T, mut b: T) -> T
where
    T: Float + RealBetaConsts + RealGamma + FloatConst + SubAssign + AddAssign,
{
    if is_gamma_pole(a) {
        return r_beta_negint(a, b);
    }
    if is_gamma_pole(b) {
        return r_beta_negint(b, a);
    }

    if a.abs() < b.abs() {
        mem::swap(&mut a, &mut b);
    }

    // a > ASYMP_FACTOR * max(|b|, 1)
    if a.abs() > T::ASYMP_FACTOR * b.abs() && a.abs() > T::ASYMP_FACTOR {
        // Avoid loss of precision in lgamma(a+b) - lgamma(a)
        return r_lbeta_asymp(a, b);
    }

    let y = a + b;
    if y.abs() > T::MAX_GAMMA || a.abs() > T::MAX_GAMMA || b.abs() > T::MAX_GAMMA {
        let sign = y.gammasgn() * a.gammasgn() * b.gammasgn();
        return sign * (a.lgamma() + b.lgamma() - y.lgamma()).exp();
    }

    let y = y.gamma();
    let a = a.gamma();
    let b = b.gamma();

    // Some care to avoid overflow..
    if (a.abs() - y.abs()).abs() > (b.abs() - y.abs()).abs() {
        return (b / y) * a;
    }
    (a / y) * b
}

/// Swap some things around to see if we can still do the computation.
fn r_beta_negint<T>(a: T, b: T) -> T
where
    T: Float + RealBetaConsts + RealGamma + FloatConst + SubAssign + AddAssign,
{
    // We are assuming that a is a negative integer.
    debug_assert!(a == a.floor() && a <= T::zero());

    if b != b.floor() || T::one() - a - b <= T::zero() {
        return T::infinity();
    }

    let sign = if cast::<T, usize>(b.floor()).unwrap() % 2 == 0 {
        T::one()
    } else {
        -T::one()
    };
    sign * r_beta(T::one() - a - b, b)
}

/// Asymptotic expansion for  ln(|B(a, b)|) for a > ASYMP_FACTOR*max(|b|, 1).
/// Taken from the Cephes library, and an unknown source.
fn r_lbeta_asymp<T>(a: T, b: T) -> T
where
    T: Float + RealGamma + RealGamma + FloatConst + SubAssign + AddAssign,
{
    let mut r = b.lgamma();
    let sign = b.gammasgn();

    let two = T::one() + T::one();
    let twelve = cast::<usize, T>(12).unwrap();

    r -= b * a.ln();
    r += b * (T::one() - b) / (two * a);
    r += b * (T::one() - b) * (T::one() - two * b) / (twelve * a * a);
    r -= b * b * (T::one() - b) * (T::one() - b) / (twelve * a * a * a);

    sign * r.exp()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::special::Factorial;

    const PRECISION: f64 = 1.0e-14;

    #[test]
    fn test_beta() {
        for a in 1_usize..10 {
            for b in 1..10 {
                assert_almost_eq!(
                    r_beta(a as f64, b as f64),
                    ((a - 1).factorial() * (b - 1).factorial()) as f64
                        / (a + b - 1).factorial() as f64,
                    PRECISION
                );
            }
        }

        // check asymptote (From SciPy v 1.10.1)
        assert_almost_eq!(r_beta(1.1E6, 0.5), 0.0016899686300445731797, PRECISION);
        assert_almost_eq!(r_beta(1.1e6, 0.7), 0.0000766158047285244695, PRECISION);
        assert_almost_eq!(r_beta(122000.64, 1.1), 0.0000024173646023778890, PRECISION);

        // Check Max Gamma?
        assert_almost_eq!(
            r_beta(165.4, 7.1),
            0.0000000000001352937563602587676443787892,
            PRECISION
        );

        // Check for negative integers
        assert_eq!(r_beta(-1.0, 1.0), -1.0);
        assert_eq!(r_beta(-2.0, 1.0), -0.5);
        assert_eq!(r_beta(-1.0, 3.0), f64::INFINITY);
    }
}
