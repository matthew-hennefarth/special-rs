//**********************************************************************
// This file is part of Sci-rs                                         *
// Copyright 2023 Matthew R. Hennefarth                                *
//**********************************************************************

#![warn(missing_docs)]
//! Sci-rs
//!
//! A scientific library written in pure Rust inspired by [SciPy].
//!
//! [SciPy]: https://scipy.org/

use num_complex::ComplexFloat;

// TODO put into some precision module file, then remove the warnings
#[allow(dead_code)]
fn is_close<T>(x: T, y: T, epsilon: <T as ComplexFloat>::Real) -> bool
where
    T: ComplexFloat,
{
    if x.is_finite() {
        return (x - y).abs() < epsilon;
    }
    if x.is_infinite() {
        return x == y;
    }
    // NaN != Nan apparently
    x.is_nan() && y.is_nan()
}

#[allow(unused_macros)]
macro_rules! assert_almost_eq {
    ($a:expr, $b:expr, $prec:expr) => {
        if !$crate::is_close($a, $b, $prec) {
            panic!(
                "assertion failed: `abs(left - right) < {:e}`, (left:
`{}`, right: `{}`)",
                $prec, $a, $b
            );
        }
    };
}

pub mod constants;
pub mod special;
pub mod traits;
