//**********************************************************************
// This file is part of Sci-rs                                         *
// Copyright 2023 Matthew R. Hennefarth                                *
//**********************************************************************

use num_traits::{cast, Float};

pub(crate) fn r_gammasgn<T>(x: T) -> T
where
    T: Float,
{
    if x > T::zero() {
        return T::one();
    }
    if x.is_zero() {
        return T::zero();
    }

    let q = x.abs();
    let p = q.floor();
    if p == q {
        return T::zero();
    }

    if cast::<T, usize>(p).unwrap() % 2_usize == 1 {
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
