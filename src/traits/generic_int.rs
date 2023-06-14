//**********************************************************************
// This file is part of Sci-rs                                         *
// Copyright 2023 Matthew R. Hennefarth                                *
//**********************************************************************
use num_traits::{PrimInt, Signed};

/// Allows me to know if a number is negative irregardless of signed or
/// unsigned
pub(crate) trait GenericInt: PrimInt {
    fn is_negative(&self) -> bool;
}

macro_rules! impl_genericint_signed {
    ($($T: ty)*) => ($(
        impl GenericInt for $T {
            #[inline(always)]
            fn is_negative(&self) -> bool {
                Signed::is_negative(self)
            }
        }
)*)
}

macro_rules! impl_genericint_unsigned {
    ($($T: ty)*) => ($(
        impl GenericInt for $T {
            #[inline(always)]
            fn is_negative(&self) -> bool { false }
        }
)*)
}

impl_genericint_signed! {i8 i16 i32 i64 isize}
impl_genericint_unsigned! {u8 u16 u32 u64 usize}
#[cfg(has_i128)]
impl_genericint_signed! {i128}
#[cfg(has_u128)]
impl_genericint_unsigned! {u128}
