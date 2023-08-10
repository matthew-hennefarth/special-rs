//**********************************************************************
// This file is part of sci-rs                                         *
// Copyright 2023 Matthew R. Hennefarth                                *
//**********************************************************************

//! Various functions related to the Gamma function.

mod gamma_trait;
pub(crate) mod gamma_util;

pub use gamma_trait::*;

mod c_gamma;
mod r_gamma;
mod r_gammainc;
mod r_gammasgn;
mod r_lgamma;
mod r_poch;
mod r_rgamma;

mod real_gamma_impl {
    pub(crate) use super::r_gamma::*;
    pub(crate) use super::r_gammainc::*;
    pub(crate) use super::r_gammasgn::*;
    pub(crate) use super::r_lgamma::*;
    pub(crate) use super::r_poch::*;
    pub(crate) use super::r_rgamma::*;
}
mod complex_gamma_impl {
    pub(crate) use super::c_gamma::*;
}
