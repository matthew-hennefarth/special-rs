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
