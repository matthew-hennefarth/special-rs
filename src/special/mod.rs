//**********************************************************************
// This file is part of Sci-rs                                         *
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
