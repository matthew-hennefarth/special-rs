use num_traits::{FromPrimitive, One, PrimInt, ToPrimitive, Zero};
use std::ops;

const MAX_MULTIPLICATIONS: usize = 16;

/// Defines the `factorial`, `factorial2`, and `factorialk` functions for `Self`. Implemented for
/// primitive integer types (usize and isize).
pub trait Factorial<Output = Self> {
    /// # Factorial function
    /// The factorial function is defined as the product of all positive integers less than or equal
    /// to $n$.
    /// $$
    /// n! = n (n-1) (n-2) \times \ldots \times 1
    /// $$
    /// Additionally we have that $0!=1$.
    ///
    /// ## Examples
    /// ```
    /// # use sci_rs::special::Factorial;
    /// assert!(0.factorial(), 1);
    /// assert!(1.factorial(), 1);
    /// assert!(2.factorial(), 2);
    /// assert!(3.factorial(), 6);
    /// ```
    fn factorial(&self) -> Self;

    /// # Factorial2 function
    /// The double factorial $n!!$ is defined as the product of all positive integers up to $n$ that have the same parity as $n$.
    /// $$
    /// n!! = n (n-2) (n-4) \times \ldots
    /// $$
    /// In the case that $n$ is even ($n = 2k$), then
    /// $$
    /// (2k)!! = (2k)(2k-2)(2k-4)\ldots(4)(2)
    /// $$
    /// In the case that $n$ is odd ($n=2k+1$), then
    /// $$
    /// (2k+1)!! = (2k+1)(2k-1)(2k-3)\ldots(3)(1)
    /// $$
    /// Additionally we have that $0!!=1$.
    fn factorial2(&self) -> Self;
    fn factorialk(&self, k: usize) -> Self;
}

impl<T> Factorial<T> for T
where
    T: Copy
        + FromPrimitive
        + ToPrimitive
        + PrimInt
        + PartialOrd
        + Zero
        + One
        + ops::Sub<Output = T>
        + ops::Mul<Output = T>,
{
    fn factorial(&self) -> Self {
        assert!(*self >= T::zero());
        const CACHE_SIZE: usize = MAX_MULTIPLICATIONS + 1;
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

    fn factorial2(&self) -> Self {
        assert!(*self >= T::zero());
        const CACHE_SIZE: usize = 2 * MAX_MULTIPLICATIONS + 1;
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

    fn factorialk(&self, k: usize) -> Self {
        assert!(k > 0);
        assert!(*self >= T::zero());
        if *self == Self::zero() {
            return Self::one();
        }

        let max_window = T::from_usize(k * MAX_MULTIPLICATIONS).unwrap();

        if *self > max_window {
            partial_product(*self - max_window, *self, k)
                * Self::factorialk(&(*self - max_window - Self::one()), k)
        } else {
            let k_as_t = T::from(k).unwrap();
            let window = k_as_t * (*self / k_as_t);
            let window = if window == *self {
                window - k_as_t
            } else {
                window
            };
            partial_product(*self - window, *self, k)
        }
    }
}

#[inline]
fn partial_product<T>(start: T, stop: T, step: usize) -> T
where
    T: Copy + PartialOrd + FromPrimitive + ToPrimitive,
{
    assert!(start <= stop);
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

    #[test]
    fn test_factorial1() {
        let k = 1;
        for i in 0..10 {
            assert_eq!(i.factorialk(k), i.factorial());
        }
    }

    #[test]
    fn test_factorial2() {
        let k = 2;
        for i in 0..15 {
            assert_eq!(i.factorialk(k), i.factorial2());
        }
    }

    #[test]
    fn test_factorial3() {
        let k = 3;
        assert_eq!(0.factorialk(k), 1);
        for i in 1..k {
            assert_eq!(i.factorialk(k), i);
        }
        assert_eq!(3.factorialk(k), 3);
        assert_eq!(4.factorialk(k), 4);
        assert_eq!(5.factorialk(k), 10);
        assert_eq!(6.factorialk(k), 18);
        assert_eq!(7.factorialk(k), 28);
    }

    #[test]
    fn test_factorial4() {
        let k = 4;
        assert_eq!(0.factorialk(k), 1);
        for i in 1..k {
            assert_eq!(i.factorialk(k), i);
        }
        assert_eq!(4.factorialk(k), 4);
        assert_eq!(5.factorialk(k), 5);
        assert_eq!(6.factorialk(k), 12);
        assert_eq!(7.factorialk(k), 21);
        assert_eq!(8.factorialk(k), 32);
    }
}
