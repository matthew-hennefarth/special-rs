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

//! Special mathematical functions
//!
//! # Available Functions
//! - Factorial, double factorial, and $k$-factorial
//! - Combinatorics (choice and permutations)
//! - Gamma and related functions
mod combinatorics;
mod factorial;
mod gamma;

pub use combinatorics::*;
pub use factorial::*;
pub use gamma::*;

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
