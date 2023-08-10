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
mod erf;
mod factorial;
mod gamma;
mod tools;
mod zigzag;

pub use bernoulli::*;
pub use beta::*;
pub use combinatorics::*;
pub use erf::*;
pub use factorial::*;
pub use gamma::*;
pub(crate) use tools::*;
pub use zigzag::*;
