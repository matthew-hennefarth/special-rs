//**********************************************************************
// This file is part of sci-rs                                         *
// Copyright 2023 Matthew R. Hennefarth                                *
//**********************************************************************

//! Various functions related to the Beta function.
//!
mod beta_trait;

pub use beta_trait::*;

mod r_beta;

mod real_beta_impl {
    pub(crate) use super::r_beta::*;
}
