use crate::special::IsNegative;
use num_traits::{CheckedAdd, CheckedMul, FromPrimitive, One, ToPrimitive};

/// Number of multiplications to perform at once.
/// This number should be optimized to do the most amount of multiplications in a single CPU
const MAX_MULTIPLICATIONS: usize = 16;

/// Factorial functions for integer-types.
///
/// Defines the `factorial`, `factorial2`, and `factorialk` functions
/// for `Self`. Implemented for primitive integer types (usize, isize,
/// etc).
///
/// ## Implementation for Primitive Integers
/// For primitive integer types, a different implementation is
/// provided for the regular `factorial`, `factorial2`, and `factorialk`
/// which does not rely on checked operations for increased speed at the
/// risk of panicking. For other types which satisfy the Factorial type,
/// the un-checked methods simply unwrap the checked results.
pub trait Factorial: Sized + CheckedMul + CheckedAdd {
    /// The factorial function is defined as the product of all positive integers less than or equal
    /// to $n$.
    /// $$
    /// n! = n (n-1) (n-2) \times \ldots \times 1
    /// $$
    /// Additionally we have that $0!=1$.
    ///
    /// # Examples
    /// ## For `usize`
    /// ```
    /// use sci_rs::special::Factorial;
    ///
    /// assert_eq!(0.factorial(), 1);
    /// assert_eq!(3.factorial(), 6);
    /// assert_eq!(5.factorial(), 120);
    /// ```
    /// ## For `isize`
    /// ```
    /// use sci_rs::special::Factorial;
    ///
    /// assert_eq!(0_isize.factorial(), 1);
    /// assert_eq!(3_isize.factorial(), 6);
    /// assert_eq!(5_isize.factorial(), 120);
    /// ```
    /// # Negative Values
    /// To avoid panicking, and to be consistent with Scipy, if provided
    /// a negative number, `0` is returned.
    fn factorial(self) -> Self {
        self.checked_factorial().unwrap()
    }

    /// Checked version of the factorial function which catches overflow.
    ///
    /// # Examples
    /// ```
    /// use sci_rs::special::Factorial;
    ///
    /// assert_eq!(2_u8.checked_factorial(), Some(2));
    /// assert_eq!(10_u8.checked_factorial(), None); // 10! overflows a u8
    /// ```
    fn checked_factorial(self) -> Option<Self>;

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
    ///
    /// # Examples
    /// ## For `usize`
    /// ```
    /// use sci_rs::special::Factorial;
    ///
    /// assert_eq!(0.factorial2(), 1);
    /// assert_eq!(3.factorial2(), 3);
    /// assert_eq!(5.factorial2(), 15);
    /// ```
    /// ## For `isize`
    /// ```
    /// use std::collections::hash_map::Values;
    /// use sci_rs::special::Factorial;
    ///
    /// assert_eq!(0_isize.factorial2(), 1);
    /// assert_eq!(3_isize.factorial2(), 3);
    /// assert_eq!(6_isize.factorial2(), 48);
    /// ```
    /// # Negative Values
    /// To avoid panicking, and to be consistent with Scipy, if provided
    /// a negative number, `0` is returned.
    fn factorial2(self) -> Self {
        self.checked_factorial2().unwrap()
    }

    /// Checked version of the `factorial2` function which catches overflow.
    ///
    /// # Examples
    /// ```
    /// use sci_rs::special::Factorial;
    ///
    /// assert_eq!(2_u8.checked_factorial2(), Some(2));
    /// assert_eq!(5_u8.checked_factorial2(), Some(15));
    /// assert_eq!(33_u8.checked_factorial2(), None); // 33!! overflows a u8
    /// ```
    fn checked_factorial2(self) -> Option<Self>;

    /// Generalized $k$-factorial.
    ///
    /// Generalization of the factorial to steps of size $k$.
    /// $$
    /// n(!!\ldots!) = (n)(n-k)(n-2k)\ldots
    /// $$
    /// We always have that $0(!!\ldots!) = 1$ irregardless of $k$.
    ///
    /// In particular for any integer $n$, we have that
    /// ```
    /// # use sci_rs::special::Factorial;
    /// # let n = 1;
    ///
    /// assert_eq!(n.factorialk(1), n.factorial());
    /// assert_eq!(n.factorialk(1), n.factorial2());
    /// ```
    ///
    /// # Examples
    /// ```
    /// use sci_rs::special::Factorial;
    ///
    /// assert_eq!(5.factorialk(3), 10); // 5 * 2
    /// assert_eq!(10.factorialk(5), 50); // 10 * 5
    /// ```
    /// # Negative Values
    /// To avoid panicking, and to be consistent with Scipy, if provided
    /// a negative number, `0` is returned.
    fn factorialk(self, k: Self) -> Self {
        self.checked_factorialk(k).unwrap()
    }

    /// Checked version of the `factorialk` function which catches overflow.
    ///
    /// # Examples
    /// ```
    /// use sci_rs::special::Factorial;
    ///
    /// assert_eq!(2_u8.checked_factorialk(2), Some(2));
    /// assert_eq!(10_u8.checked_factorialk(4), Some(120));
    /// assert_eq!(31_u8.checked_factorialk(3), None); // 33!! overflows a u8
    /// ```
    fn checked_factorialk(self, k: Self) -> Option<Self>;
}

// Cached factorial and factorial2 values
const FACTORIAL_CACHE_LEN: usize = MAX_MULTIPLICATIONS + 1;
const FACTORIAL_CACHE: [u64; FACTORIAL_CACHE_LEN] = [
    1, 1, 2, 6, 24, 120, 720, 5040, 40320, 362880, 3628000, 39916800, 479001600, 6227020800,
    87178291200, 1307674368000, 20922789888000,
];

const FACTORIAL2_CACHE_LEN: usize = 2 * MAX_MULTIPLICATIONS + 1;
const FACTORIAL2_CACHE: [u64; FACTORIAL2_CACHE_LEN] = [
    1, 1, 2, 3, 8, 15, 48, 105, 384, 945, 3840, 10395, 46080, 135135, 645120, 2027025, 10321920,
    34459425, 185794560, 654729075, 3715891200, 13749310575, 81749606400, 316234143225,
    1961990553600, 7905853580625, 51011754393600, 213458046676875, 1428329123020800,
    6190283353629375, 42849873690624000, 191898783962510625, 1371195958099968000,
];

macro_rules! factorial_primint_impl {
    ($($T: ty)*) => ($(
        impl Factorial for $T {
            fn factorial(self) -> Self {
                if self.is_negative() {
                    return 0;
                }
                let cache_as_type = FACTORIAL_CACHE_LEN as $T;
                if self < cache_as_type {
                    return FACTORIAL_CACHE[self as usize].try_into().unwrap();
                }
                partial_product(self - (cache_as_type - 1), self, 1) * (self - cache_as_type).factorial()
            }

            fn checked_factorial(self) -> Option<Self> {
                if self.is_negative() {
                    return Some(0);
                }
                let cache_as_type = FACTORIAL_CACHE_LEN as $T;
                if self < cache_as_type {
                    return FACTORIAL_CACHE[self as usize].try_into().ok();
                }
                (checked_partial_product(self - (cache_as_type - 1), self, 1)?).checked_mul(
                    (self - cache_as_type).checked_factorial()?)
            }

            fn factorial2(self) -> Self {
                if self.is_negative() {
                    return 0;
                }
                let cache_as_type = FACTORIAL2_CACHE_LEN as $T;
                if self < cache_as_type {
                    return FACTORIAL2_CACHE[self as usize].try_into().unwrap();
                }
                partial_product(self - (cache_as_type - 1), self, 2) * (self - cache_as_type).factorial2()
            }

            fn checked_factorial2(self) -> Option<Self> {
                if self.is_negative() {
                    return Some(0);
                }
                let cache_as_type = FACTORIAL2_CACHE_LEN as $T;
                if self < cache_as_type {
                    return FACTORIAL2_CACHE[self as usize].try_into().ok();
                }
                (checked_partial_product(self - (cache_as_type - 1), self, 2)?).checked_mul(
                    (self - cache_as_type).factorial2())
            }

            fn factorialk(self, k: Self) -> Self {
                if self.is_negative() {
                    return 0;
                }
                assert!(k > 0);
                if self == 0 {
                    return 1;
                }
                let max_window = k * MAX_MULTIPLICATIONS as $T;

                if self > max_window {
                    return partial_product(self - max_window, self, k) * (self - max_window - 1).factorialk(k);
                }

                let k_as_t = k as $T;
                let window = k_as_t * (self / k_as_t);
                let window = if window == self {
                    window - k_as_t
                } else {
                    window
                };
                partial_product(self - window, self, k)

            }

            fn checked_factorialk(self, k: $T) -> Option<Self> {
                if self.is_negative() {
                    return Some(0);
                }
                assert!(k > 0);
                if self == 0 {
                    return Some(1);
                }

                let max_window = k.checked_mul(MAX_MULTIPLICATIONS as $T)?;

                if self > max_window {
                    return (checked_partial_product(self - max_window, self, k)?).checked_mul(
                        (self - max_window - 1).checked_factorialk(k)?);
                }
                let k_as_t = k as $T;
                let window = k_as_t * (self / k_as_t);
                let window = if window == self {
                    window - k_as_t
                } else {
                    window
                };
                checked_partial_product(self - window, self, k)
            }
        }
    )*)
}

factorial_primint_impl! {u8 u16 u32 u64 usize i8 i16 i32 i64 isize}
#[cfg(has_i128)]
factorial_primint_impl! {u128 i128}

/// Computes the checked product between `start` and `stop` stepping
/// with `step`. Catches overflow from multiplication or addition
///
/// $$
/// n_0\times n_1\times\ldots \times n_k
/// $$
/// where $n_0 = $ `start`, $n_{i+1} = \text{step}\times n_i$, and $n_k$
/// is the largest value such that $n_k < $ `stop + 1`./
#[inline(always)]
fn checked_partial_product<T>(mut start: T, stop: T, step: T) -> Option<T>
where
    T: Sized + PartialOrd + FromPrimitive + ToPrimitive + One + CheckedMul + CheckedAdd + Copy,
{
    let mut result = Some(T::one());
    while start <= stop {
        result = (result?).checked_mul(&start);
        start = start.checked_add(&step)?;
    }
    result
}

/// Computes the product between `start` and `stop` stepping with `step`.
///
/// $$
/// n_0\times n_1\times\ldots \times n_k
/// $$
/// where $n_0 = $ `start`, $n_{i+1} = \text{step}\times n_i$, and $n_k$
/// is the largest value such that $n_k < $ `stop + 1`./
#[inline(always)]
fn partial_product<T>(start: T, stop: T, step: T) -> T
where
    T: PartialOrd + FromPrimitive + ToPrimitive,
{
    assert!(start <= stop);
    T::from_isize(
        (start.to_isize().unwrap()..stop.to_isize().unwrap() + 1)
            .step_by(step.to_usize().unwrap())
            .product::<isize>(),
    )
    .unwrap()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn factorial_u8() {
        assert_eq!(0_u8.factorial(), 1);
        for i in 1_u8..3 {
            assert_eq!(i.factorial(), i);
            assert_eq!(i.checked_factorial(), Some(i.factorial()));
        }
        assert_eq!(4_u8.factorial(), 24);
        assert_eq!(5_u8.factorial(), 120);
        assert_eq!(6_u8.checked_factorial(), None);
        assert_eq!(7_u8.checked_factorial(), None);
    }

    #[test]
    fn factorial_usize() {
        assert_eq!(0_usize.factorial(), 1);
        for i in 1_usize..3 {
            assert_eq!(i.factorial(), i);
            assert_eq!(i.checked_factorial(), Some(i.factorial()));
        }
        assert_eq!(4_usize.factorial(), 24);
        assert_eq!(13_usize.factorial(), 6_227_020_800);
        assert_eq!(15_usize.factorial(), 1_307_674_368_000);
        assert_eq!(20_usize.factorial(), 2_432_902_008_176_640_000);
        assert_eq!(5_usize.checked_factorial(), Some(120));
        assert_eq!(25_usize.checked_factorial(), None);
    }

    #[test]
    fn factorial_i64() {
        assert_eq!(0_i64.factorial(), 1);
        for i in 1_isize..3 {
            assert_eq!(i.factorial(), i);
            assert_eq!(i.checked_factorial(), Some(i.factorial()));
        }
        assert_eq!(4_i64.factorial(), 24);
        assert_eq!(13_i64.factorial(), 6_227_020_800);
        assert_eq!(15_i64.factorial(), 1_307_674_368_000);
        assert_eq!(20_i64.factorial(), 2_432_902_008_176_640_000);
    }

    #[test]
    fn factorial2_u8() {
        assert_eq!(0_u8.factorial2(), 1);
        for i in 1_u8..4 {
            assert_eq!(i.factorial2(), i);
            assert_eq!(i.checked_factorial2(), Some(i.factorial2()));
        }
        assert_eq!(5_u8.checked_factorial2(), Some(15));
        assert_eq!(7_u8.checked_factorial2(), Some(105));
        assert_eq!(8_u8.checked_factorial2(), None);
        assert_eq!(33_u8.checked_factorial2(), None);
    }

    #[test]
    fn factorial2_usize() {
        assert_eq!(0_usize.factorial2(), 1);
        for i in 1_usize..4 {
            assert_eq!(i.factorial2(), i);
            assert_eq!(i.checked_factorial2(), Some(i.factorial2()));
        }
        assert_eq!(4_usize.factorial2(), 8);
        assert_eq!(13_usize.factorial2(), 135_135); // Taken from WolframAlpha
        assert_eq!(22_usize.factorial2(), 81_749_606_400); // Taken from WolframAlpha
        assert_eq!(23_usize.factorial2(), 316_234_143_225); // Taken from WolframAlpha
        assert_eq!(100_usize.checked_factorial2(), None);
        assert_eq!(200_usize.checked_factorial2(), None);
    }

    #[test]
    fn factorial2_i64() {
        assert_eq!(0_i64.factorial2(), 1);
        for i in 1_i64..4 {
            assert_eq!(i.factorial2(), i);
            assert_eq!(i.checked_factorial2(), Some(i.factorial2()));
        }
        assert_eq!(4_i64.factorial2(), 8);
        assert_eq!(13_i64.factorial2(), 135_135); // Taken from WolframAlpha
        assert_eq!(22_i64.factorial2(), 81_749_606_400); // Taken from WolframAlpha
        assert_eq!(23_i64.factorial2(), 316_234_143_225); // Taken from WolframAlpha
    }

    #[test]
    fn factorial1() {
        let k = 1;
        for i in 0..10 {
            assert_eq!(i.factorialk(k), i.factorial());
            assert_eq!(i.checked_factorialk(k), Some(i.factorialk(k)));
        }
    }

    #[test]
    fn factorial2() {
        let k = 2;
        for i in 0..15 {
            assert_eq!(i.factorialk(k), i.factorial2());
            assert_eq!(i.checked_factorialk(k), Some(i.factorialk(k)));
        }
    }

    #[test]
    fn factorial3() {
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
    fn factorial4() {
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
