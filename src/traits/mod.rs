//**********************************************************************
// This file is part of Sci-rs                                         *
// Copyright 2023 Matthew R. Hennefarth                                *
//**********************************************************************

//! Various numerical traits which often extends those present in
//! [num_traits].

mod float_const;
mod generic_int;

pub use float_const::*;
pub(crate) use generic_int::*;
