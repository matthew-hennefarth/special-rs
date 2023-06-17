//**********************************************************************
// This file is part of sci-rs                                         *
// Copyright 2023 Matthew R. Hennefarth                                *
//**********************************************************************

use crate::special::beta_util::{r_lbeta_asymp, RealBetaConsts};
use crate::special::gamma_util::is_gamma_pole;
use crate::special::RealGamma;
use num_traits::{Float, FloatConst};
use std::mem;
use std::ops::{AddAssign, SubAssign};

pub(crate) fn r_lbeta<T>(mut a: T, mut b: T) -> T
where
    T: Float + FloatConst + RealBetaConsts + RealGamma + SubAssign + AddAssign,
{
    if is_gamma_pole(a) {
        return r_lbeta_negint(a, b);
    }
    if is_gamma_pole(b) {
        return r_lbeta_negint(b, a);
    }

    if a.abs() < b.abs() {
        mem::swap(&mut a, &mut b);
    }

    if a.abs() > T::ASYMP_FACTOR * b.abs() && a > T::ASYMP_FACTOR {
        return r_lbeta_asymp(a, b);
    }

    let y = a + b;
    if y.abs() > T::MAX_GAMMA || a.abs() > T::MAX_GAMMA || b.abs() > T::MAX_GAMMA {
        let sign = y.gammasgn() * a.gammasgn() * b.gammasgn();
        return sign * (a.lgamma() + b.lgamma() - y.lgamma());
    }

    let y = y.gamma();
    a = a.gamma();
    b = b.gamma();
    if y.is_zero() {
        return T::infinity();
    }
    if (a.abs() - y.abs()).abs() > (b.abs() - y.abs()).abs() {
        ((b / y) * a).abs().ln()
    } else {
        ((a / y) * b).abs().ln()
    }
}

fn r_lbeta_negint<T>(a: T, b: T) -> T
where
    T: Float + FloatConst + RealBetaConsts + RealGamma + SubAssign + AddAssign,
{
    debug_assert!(a.floor() == a && a <= T::zero());
    if b != b.floor() || T::one() - a - b <= T::zero() {
        return T::infinity();
    }
    r_lbeta(T::one() - a - b, b)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::special::Factorial;

    const PRECISION: f64 = 1.0e-14;

    #[test]
    fn test_r_lbeta() {
        for a in 1_usize..10 {
            for b in 1..10 {
                assert_almost_eq!(
                    r_lbeta(a as f64, b as f64),
                    (((a - 1).factorial() * (b - 1).factorial()) as f64).ln()
                        - ((a + b - 1).factorial() as f64).ln(),
                    PRECISION
                );
            }

            // check asymptote (From SciPy v 1.10.1)
            assert_almost_eq!(r_lbeta(1.1E6, 0.5), -6.3830453123232357981465, PRECISION);
            assert_almost_eq!(r_lbeta(0.5, 1.1E6), -6.3830453123232357981465, PRECISION);
            assert_almost_eq!(r_lbeta(1.1e6, 0.7), -9.4767071744518105447241, PRECISION);
            assert_almost_eq!(
                r_lbeta(122000.64, 1.1),
                -12.9328326184768229722977,
                PRECISION
            );
            assert_almost_eq!(
                r_lbeta(142346.1124, 11551.51221),
                -41022.7405061810277402400970,
                PRECISION
            );

            // Check for negative integers
            assert_eq!(r_lbeta(-1.0, 1.0), 0.0);
            assert_almost_eq!(r_lbeta(-2.0, 1.0), -0.6931471805599452862268, PRECISION);
            assert_almost_eq!(r_lbeta(1.0, -2.0), -0.6931471805599452862268, PRECISION);
            assert_eq!(r_lbeta(-1.0, 3.0), f64::INFINITY);
        }
    }
}
