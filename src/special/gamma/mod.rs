//**********************************************************************
// This file is part of Sci-rs                                         *
// Copyright 2023 Matthew R. Hennefarth                                *
//**********************************************************************

//! Various functions related to the gamma function.

mod gamma_trait;
mod gamma_util;

pub use gamma_trait::*;
use gamma_util::*;

// Functions which assume a real-valued input
mod r_gamma;
mod r_gammasgn;
mod r_lgamma;
mod r_poch;
mod r_rgamma;

use r_gamma::*;
use r_gammasgn::*;
use r_lgamma::*;
use r_poch::*;
use r_rgamma::*;

mod c_gamma;
use c_gamma::*;
