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

use num_traits::{cast, Float};

pub(crate) fn r_gammasgn<T>(x: T) -> T
where
    T: Float,
{
    if x > T::zero() {
        return T::one();
    }
    if x.is_zero() {
        return T::zero();
    }

    let q = x.abs();
    let p = q.floor();
    if p == q {
        return T::zero();
    }

    if cast::<T, usize>(p).unwrap() % 2_usize == 1 {
        T::one()
    } else {
        -T::one()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_gammasgn() {
        for i in 1..20 {
            assert_eq!(r_gammasgn(i as f32), 1.0);
            assert_eq!(r_gammasgn(i as f32 + 0.5), 1.0);
            assert_eq!(r_gammasgn(i as f32 + 0.25), 1.0);
        }
        // Should be negative
        assert_eq!(r_gammasgn(-0.5), -1.0);
        assert_eq!(r_gammasgn(-2.5), -1.0);
        assert_eq!(r_gammasgn(-4.5), -1.0);
        assert_eq!(r_gammasgn(-6.5), -1.0);

        // Should be positive
        assert_eq!(r_gammasgn(-1.5), 1.0);
        assert_eq!(r_gammasgn(-3.5), 1.0);
        assert_eq!(r_gammasgn(-5.5), 1.0);
        assert_eq!(r_gammasgn(-7.5), 1.0);
    }
}
