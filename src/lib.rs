//**********************************************************************
// This file is part of sci-rs                                         *
//                                                                     *
// See the LICENSE files at the top-level directory of this            *
// distribution.                                                       *
//                                                                     *
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or   *
// http://www.apache.org/licenses/LICENSE-2.0> or the MIT license      *
// <LICENSE-MIT or http://opensource.org/licenses/MIT>, at your        *
// option. This file may not be copied, modified, or distributed       *
// except according to those terms.                                    *
//                                                                     *
// Copyright 2023 Matthew R. Hennefarth                                *
//**********************************************************************

#![warn(missing_docs)]
#![doc(test(attr(deny(warnings))))]

//! # sci_rs
//!
//! A scientific library written in pure Rust inspired by [SciPy].
//!
//! # Features
//! Working on!
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
