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
//! Includes all standard library constants, but includes additional
//! important quantities.
//!
//! - $\sqrt{\tau}$
//! - $\sqrt{\pi}$
//! - $\ln{\pi}$

pub mod f64 {
    //! Extended mathematical constants for the `f64` double-precision
    //! floating point type.

    /// Euler's number $e$
    ///
    /// Sometimes called Napier's constant. Value taken from OEIS: [A001113]
    ///
    /// [A001113]: https://oeis.org/A001113
    pub const E: f64 = 2.71828182845904523536028747135_f64;

    /// Archimedes's constant $\pi$
    ///
    /// Value taken from OEIS: [A000796]
    ///
    /// [A000796]: https://oeis.org/A000796
    pub const PI: f64 = 3.14159265358979323846264338327_f64;

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

    /// $\ln\pi$
    ///
    /// Value taken from OEIS: [A053510]
    ///
    /// [A053510]: https://oeis.org/A053510
    pub const LOG_PI: f64 = 1.14472988584940017414342735135_f64;

    /// $\ln\sqrt{2\pi}$
    ///
    /// Value taken from [Wolframalpha].
    ///
    /// [Wolframalpha]: https://www.wolframalpha.com/input?i=ln+sqrt%282pi%29%29
    pub const LOG_SQRT_2_PI: f64 = 0.91893853320467274178032973640_f64;
}

pub mod f32 {
    //! Extended mathematical constants for the `f32` double-precision
    //! floating point type.

    /// Euler's number $e$
    ///
    /// Sometimes called Napier's constant. Value taken from OEIS: [A001113]
    ///
    /// [A001113]: https://oeis.org/A001113
    pub const E: f32 = 2.71828182845904523536028747135_f32;

    /// Archimedes's constant $\pi$
    ///
    /// Value taken from OEIS: [A000796]
    ///
    /// [A000796]: https://oeis.org/A000796
    pub const PI: f32 = 3.14159265358979323846264338327_f32;

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

    /// $\log\pi$
    ///
    /// Value taken from OEIS: [A053510]
    ///
    /// [A053510]: https://oeis.org/A053510
    pub const LOG_PI: f32 = 1.14472988584940017414342735135_f32;

    /// $\ln\sqrt{2\pi}$
    ///
    /// Value taken from [Wolframalpha].
    ///
    /// [Wolframalpha]: https://www.wolframalpha.com/input?i=ln+sqrt%282pi%29%29
    pub const LOG_SQRT_2_PI: f32 = 0.91893853320467274178032973640_f32;
}
