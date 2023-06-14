//**********************************************************************
// This file is part of Sci-rs                                         *
//                                                                     *
// Sci-rs is licensed under the Apache License, Version 2.0 (the       *
// "License"); you may not use this file except in compliance with the *
// License. You may obtain a copy of the License at                    *
//                                                                     *
//     http://www.apache.org/licenses/LICENSE-2.0                      *
//                                                                     *
// Unless required by applicable law or agreed to in writing, software *
// distributed under the License is distributed on an "AS IS" BASIS,   *
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or     *
// implied. See the License for the specific language governing        *
// permissions and limitations under the License.                      *
//                                                                     *
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
