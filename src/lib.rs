#![warn(missing_docs)]
#![warn(missing_doc_code_examples)]

fn is_close<T>(x: T, y: T, epsilon: T) -> bool
where
    T: num_traits::Float,
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

#[macro_use]
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

pub mod special;

// mod preamble {
//     pub use crate::*;
// }

//use crate::preamble::*;
