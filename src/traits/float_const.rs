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

use num_traits::FloatConst;

macro_rules! constant {
    ($( $method:ident () -> $ret:expr ; )*)
        => {$(
            #[inline]
            fn $method() -> Self {
                $ret
            }
        )*};
}

macro_rules! float_const_impl {
    ($(#[$doc:meta] $constant:ident,)+) => (
        #[allow(non_snake_case)]
        /// Numeric traits for floating points which implement constant values from [crate::constants].
        ///
        /// Different implementation of the [num_traits::FloatConst] trait to include the additional constants in the [sci_rs] library.
        ///
        /// [sci_rs]: crate
        pub trait FloatSciConst : FloatConst {
            $(#[$doc] fn $constant() -> Self;)+
        }
        float_const_impl! { @float f32, $($constant,)+ }
        float_const_impl! { @float f64, $($constant,)+ }
    );
    (@float $T:ident, $($constant:ident,)+) => (
        impl FloatSciConst for $T {
            constant! {
                $( $constant() -> crate::constants::$T::$constant; )+
            }
        }
    );
}

float_const_impl! {
    #[doc = "Returns $\\sqrt{\\tau}$"]
    SQRT_TAU,
    #[doc = "Returns $\\sqrt{\\pi}$"]
    SQRT_PI,
    #[doc = "Returns $\\ln{\\pi}$"]
    LOG_PI,
    #[doc= "Returns $\\ln\\sqrt{2\\pi}$"]
    LOG_SQRT_2_PI,
}
