use num_traits::{FromPrimitive, ToPrimitive, Zero};
use std::ops;

pub trait Factorial<Output = Self> {
    fn factorial2(&self) -> Self;
    fn factorial(&self) -> Self;
}

impl<T> Factorial<T> for T
where
    T: Copy
        + FromPrimitive
        + ToPrimitive
        + PartialOrd
        + Zero
        + ops::Sub<Output = T>
        + ops::Mul<Output = T>,
{
    fn factorial2(&self) -> Self {
        assert!(*self >= T::zero());
        const CACHE_SIZE: usize = 33;
        const FACTORIAL2_CACHE: [usize; CACHE_SIZE] = [
            1, 1, 2, 3, 8, 15, 48, 105, 384, 945, 3840, 10395, 46080, 135135, 645120, 2027025,
            10321920, 34459425, 185794560, 654729075, 3715891200, 13749310575, 81749606400,
            316234143225, 1961990553600, 7905853580625, 51011754393600, 213458046676875,
            1428329123020800, 6190283353629375, 42849873690624000, 191898783962510625,
            1371195958099968000,
        ];
        if *self < T::from_usize(CACHE_SIZE).unwrap() {
            return T::from_usize(FACTORIAL2_CACHE[T::to_usize(self).unwrap()]).unwrap();
        }
        partial_product(*self - T::from_usize(CACHE_SIZE - 1).unwrap(), *self, 2)
            * Self::factorial2(&(*self - T::from_usize(CACHE_SIZE).unwrap()))
    }

    fn factorial(&self) -> Self {
        assert!(*self >= T::zero());
        const CACHE_SIZE: usize = 17;
        const FACTORIAL_CACHE: [usize; CACHE_SIZE] = [
            1, 1, 2, 6, 24, 120, 720, 5040, 40320, 362880, 3628000, 39916800, 479001600,
            6227020800, 87178291200, 1307674368000, 20922789888000,
        ];
        if *self < T::from_usize(CACHE_SIZE).unwrap() {
            return T::from_usize(FACTORIAL_CACHE[T::to_usize(self).unwrap()]).unwrap();
        }

        partial_product(*self - T::from_usize(CACHE_SIZE - 1).unwrap(), *self, 1)
            * Self::factorial(&(*self - T::from_usize(CACHE_SIZE).unwrap()))
    }
}

#[inline]
fn partial_product<T>(start: T, stop: T, step: usize) -> T
where
    T: Copy + PartialOrd + FromPrimitive + ToPrimitive,
{
    assert!(start < stop);
    T::from_usize(
        (T::to_usize(&start).unwrap()..T::to_usize(&stop).unwrap() + 1)
            .step_by(step)
            .product::<usize>(),
    )
    .unwrap()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn factorial_usize() {
        assert_eq!(0_usize.factorial(), 1);
        for i in 1_usize..3 {
            assert_eq!(i.factorial(), i);
        }
        assert_eq!(4_usize.factorial(), 24);
        assert_eq!(13_usize.factorial(), 6_227_020_800);
        assert_eq!(15_usize.factorial(), 1_307_674_368_000);
        assert_eq!(20_usize.factorial(), 2_432_902_008_176_640_000);
    }

    #[test]
    fn factorial_i64() {
        assert_eq!(0_i64.factorial(), 1);
        for i in 1_isize..3 {
            assert_eq!(i.factorial(), i);
        }
        assert_eq!(4_i64.factorial(), 24);
        assert_eq!(13_i64.factorial(), 6_227_020_800);
        assert_eq!(15_i64.factorial(), 1_307_674_368_000);
        assert_eq!(20_i64.factorial(), 2_432_902_008_176_640_000);
    }

    #[test]
    fn factorial2_usize() {
        assert_eq!(0_usize.factorial2(), 1);
        for i in 1_usize..4 {
            assert_eq!(i.factorial2(), i);
        }
        assert_eq!(4_usize.factorial2(), 8);
        assert_eq!(13_usize.factorial2(), 135_135); // Taken from WolframAlpha
        assert_eq!(22_usize.factorial2(), 81_749_606_400); // Taken from WolframAlpha
        assert_eq!(23_usize.factorial2(), 316_234_143_225); // Taken from WolframAlpha
    }

    #[test]
    fn factorial2_i64() {
        assert_eq!(0_i64.factorial2(), 1);
        for i in 1_i64..4 {
            assert_eq!(i.factorial2(), i);
        }
        assert_eq!(4_i64.factorial2(), 8);
        assert_eq!(13_i64.factorial2(), 135_135); // Taken from WolframAlpha
        assert_eq!(22_i64.factorial2(), 81_749_606_400); // Taken from WolframAlpha
        assert_eq!(23_i64.factorial2(), 316_234_143_225); // Taken from WolframAlpha
    }
}