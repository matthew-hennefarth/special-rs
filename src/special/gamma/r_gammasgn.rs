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

            let neg_even = -(2 * i) as f32;
            assert_eq!(r_gammasgn(neg_even - 0.5), -1.0);
            assert_eq!(r_gammasgn(neg_even - 0.25), -1.0);

            let neg_odd = -(2 * i + 1) as f32;
            assert_eq!(r_gammasgn(neg_odd - 0.5), 1.0);
            assert_eq!(r_gammasgn(neg_odd - 0.25), 1.0);
        }
    }
}
