//**********************************************************************
// This file is part of sci-rs                                         *
// Copyright 2023 Matthew R. Hennefarth                                *
//**********************************************************************

use crate::special::gamma::gamma_util::is_gamma_pole;
use crate::traits::FloatSciConst;
use num_traits::{cast, Float};

pub(crate) fn r_gammasgn<T>(x: T) -> T
where
    T: Float + FloatSciConst,
{
    if !x.is_finite() {
        return x;
    }
    if is_gamma_pole(x) {
        return T::zero();
    }

    if x.is_sign_positive() {
        return T::one();
    }
    if cast::<T, usize>(x.abs().floor()).unwrap() % 2_usize == 1 {
        T::one()
    } else {
        -T::one()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_gammasgn() {
        for i in 1..20 {
            assert_eq!(r_gammasgn(i as f32), 1.0);
            assert_eq!(r_gammasgn(i as f32 + 0.5), 1.0);
            assert_eq!(r_gammasgn(i as f32 + 0.25), 1.0);
        }
        // Should be negative
        assert_eq!(r_gammasgn(-0.5), -1.0);
        assert_eq!(r_gammasgn(-2.5), -1.0);
        assert_eq!(r_gammasgn(-4.5), -1.0);
        assert_eq!(r_gammasgn(-6.5), -1.0);

        // Should be positive
        assert_eq!(r_gammasgn(-1.5), 1.0);
        assert_eq!(r_gammasgn(-3.5), 1.0);
        assert_eq!(r_gammasgn(-5.5), 1.0);
        assert_eq!(r_gammasgn(-7.5), 1.0);
    }
}
