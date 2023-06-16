//**********************************************************************
// This file is part of sci-rs                                         *
// Copyright 2023 Matthew R. Hennefarth                                *
//**********************************************************************

//! Special mathematical functions
//!
//! # Available Functions
//! - Factorial, double factorial, and $k$-factorial
//! - Combinatorics (choice and permutations)
//! - Gamma and related functions
mod bernoulli;
mod beta;
mod combinatorics;
mod factorial;
mod gamma;
mod tangent_num;

pub use beta::*;
pub use combinatorics::*;
pub use factorial::*;
pub use gamma::*;
pub use tangent_num::*;
