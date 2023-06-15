//**********************************************************************
// This file is part of sci-rs                                         *
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
    #[doc= "Returns $\\ln\\sqrt{2\\pi}$"]
    LOG_SQRT_2_PI,
    #[doc= "Returns $\\gamma$"]
    GAMMA,
}
