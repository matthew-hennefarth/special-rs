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
//! Rust library for various scientific applications. Inspiration taken
//! from Scipy. The focus is to implement everything in pure rust with
//! minimal outside dependencies.

// fn is_close<T>(x: T, y: T, epsilon: T) -> bool
// where
//     T: num_traits::Float,
// {
//     if x.is_finite() {
//         return (x - y).abs() < epsilon;
//     }
//     if x.is_infinite() {
//         return x == y;
//     }
//     // NaN != Nan apparently
//     x.is_nan() && y.is_nan()
// }

// macro_rules! assert_almost_eq {
//     ($a:expr, $b:expr, $prec:expr) => {
//         if !$crate::is_close($a, $b, $prec) {
//             panic!(
//                 "assertion failed: `abs(left - right) < {:e}`, (left:
// `{}`, right: `{}`)",
//                 $prec, $a, $b
//             );
//         }
//     };
// }

//pub(crate) use assert_almost_eq;

pub mod special;

// mod preamble {
//     pub use crate::*;
// }

//use crate::preamble::*;
