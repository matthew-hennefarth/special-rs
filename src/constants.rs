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

//! Various important physical and mathematical constants.
//!
//! Values are taken primarily from the On-Line Encyclopedia of Integer Sequences ([OEIS]) and implemented for both `f64` and `f32` types. Some constants may overlap with those in the standard library but are present here to reduce dependence on the standard library.
//!
//! [OEIS]: https://oeis.org/

pub mod f64 {
    //! Extended mathematical constants for the `f64` double-precision
    //! floating point type.

    pub use std::f64::consts::*;

    /// $\sqrt{2\pi} = \sqrt{\tau}$
    ///
    /// Value taken from OEIS: [A019727]
    ///
    /// [A019727]: https://oeis.org/A019727
    pub const SQRT_TAU: f64 = 2.50662827463100050241576528481_f64;

    /// $\sqrt{\pi}$
    ///
    /// Value taken from OEIS: [A002161]
    ///
    /// [A002161]: https://oeis.org/A002161
    pub const SQRT_PI: f64 = 1.77245385090551602729816748334_f64;

    /// $\ln\sqrt{2\pi}$
    ///
    /// Value taken from [Wolframalpha].
    ///
    /// [Wolframalpha]: https://www.wolframalpha.com/input?i=ln+sqrt%282pi%29%29
    pub const LOG_SQRT_2_PI: f64 = 0.91893853320467274178032973640_f64;

    /// Euler-Mascheroni constant ($\gamma$)
    ///
    /// Value taken from OEIS: [A001620]
    ///
    /// [A001620]: https://oeis.org/A001620
    pub const GAMMA: f64 = 0.577215664901532860606512090082402431042_f64;
}

pub mod f32 {
    //! Extended mathematical constants for the `f32` double-precision
    //! floating point type.

    pub use std::f32::consts::*;

    /// $\sqrt{2\pi} = \sqrt{\tau}$
    ///
    /// Value taken from OEIS: [A019727]
    ///
    /// [A019727]: https://oeis.org/A019727
    pub const SQRT_TAU: f32 = 2.50662827463100050241576528481_f32;

    /// $\sqrt{\pi}$
    ///
    /// Value taken from OEIS: [A002161]
    ///
    /// [A002161]: https://oeis.org/A002161
    pub const SQRT_PI: f32 = 1.77245385090551602729816748334_f32;

    /// $\ln\sqrt{2\pi}$
    ///
    /// Value taken from [Wolframalpha].
    ///
    /// [Wolframalpha]: https://www.wolframalpha.com/input?i=ln+sqrt%282pi%29%29
    pub const LOG_SQRT_2_PI: f32 = 0.91893853320467274178032973640_f32;

    /// Euler-Mascheroni constant ($\gamma$)
    ///
    /// Value taken from OEIS: [A001620]
    ///
    /// [A001620]: https://oeis.org/A001620
    pub const GAMMA: f32 = 0.577215664901532860606512090082402431042_f32;
}
