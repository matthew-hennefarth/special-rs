//! Special mathematical functions
//!
//! # Available Functions
//! - Factorial, double factorial, and k-factorial
//! - Combinatorics
//!
mod combinatorics;
mod factorial;

pub use combinatorics::Combinatorics;
pub use factorial::Factorial;

/// Allows me to know if a number is negative irregardless of signed or
/// unsigned
pub(crate) trait IsNegative {
    fn is_negative(&self) -> bool;
}

macro_rules! impl_isnegative_signed {
    ($($T: ty)*) => ($(
        impl IsNegative for $T {
            #[inline(always)]
            fn is_negative(&self) -> bool {
                (*self) < 0
            }
        }
)*)
}

macro_rules! impl_isnegative_unsigned {
    ($($T: ty)*) => ($(
        impl IsNegative for $T {
            #[inline(always)]
            fn is_negative(&self) -> bool { false }
        }
)*)
}

impl_isnegative_signed! {i8 i16 i32 i64 isize}
impl_isnegative_unsigned! {u8 u16 u32 u64 usize}
#[cfg(has_i128)]
impl_isnegative_signed! {i128}
#[cfg(has_i128)]
impl_isnegative_unsigned! {u128}
