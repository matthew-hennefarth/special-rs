//**********************************************************************
// This file is part of sci-rs                                         *
// Copyright 2023 Matthew R. Hennefarth                                *
//**********************************************************************

//! Various functions related to the Beta function.
//!
mod beta_trait;
pub(crate) mod beta_util;

pub use beta_trait::*;

mod r_beta;
mod r_lbeta;

mod real_beta_impl {
    pub(crate) use super::r_beta::*;
    pub(crate) use super::r_lbeta::*;
}
