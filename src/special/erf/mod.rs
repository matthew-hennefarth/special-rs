//**********************************************************************
// This file is part of sci-rs                                         *
// Copyright 2023 Matthew R. Hennefarth                                *
//**********************************************************************

//! Various functions related to the Error function.

mod erf_trait;

pub use erf_trait::*;

mod r_erf;

mod real_erf_impl {
    pub(crate) use super::r_erf::*;
}
