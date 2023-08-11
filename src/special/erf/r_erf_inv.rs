//**********************************************************************
// This file is part of sci-rs                                         *
// Copyright 2023 Matthew R. Hennefarth                                *
//**********************************************************************

use num_traits::Float;

fn r_erf_impl<T>(p: T, q: T) -> T
where
    T: Float,
{
    todo!();
}

pub(crate) fn r_erf_inv<T>(x: T) -> T
where
    T: Float,
{
    if x < -T::one() || x > T::one() {
        panic!("Domain error");
    }
    if x.is_one() {
        return T::infinity();
    }
    if x == -T::one() {
        return T::neg_infinity();
    }
    if x.is_zero() {
        return T::zero();
    }

    if x.is_sign_negative() {
        -r_erf_impl(-x, T::one() + x)
    } else {
        r_erf_impl(x, T::one() - x)
    }
}

pub(crate) fn r_erfc_inv<T>(x: T) -> T
where
    T: Float,
{
    let two = T::from(2.0).unwrap();
    if x < T::zero() || x > two {
        panic!("Domain error");
    }
    if x.is_zero() {
        return T::infinity();
    }
    if x == two {
        return T::neg_infinity();
    }
    if x.is_one() {
        return T::zero();
    }

    if x > T::one() {
        -r_erf_impl(x - T::one(), two - x)
    } else {
        r_erf_impl(T::one() - x, x)
    }
}
