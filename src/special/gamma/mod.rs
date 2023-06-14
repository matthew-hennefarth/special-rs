//**********************************************************************
// This file is part of Sci-rs                                         *
//                                                                     *
// Sci-rs is licensed under the Apache License, Version 2.0 (the       *
// "License"); you may not use this file except in compliance with the *
// License. You may obtain a copy of the License at                    *
//                                                                     *
//     http://www.apache.org/licenses/LICENSE-2.0                      *
//                                                                     *
// Unless required by applicable law or agreed to in writing, software *
// distributed under the License is distributed on an "AS IS" BASIS,   *
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or     *
// implied. See the License for the specific language governing        *
// permissions and limitations under the License.                      *
//                                                                     *
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
