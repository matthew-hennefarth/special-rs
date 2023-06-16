//**********************************************************************
// This file is part of sci-rs                                         *
// Copyright 2023 Matthew R. Hennefarth                                *
//**********************************************************************

use num_traits::{FromPrimitive, PrimInt};

/// Number of multiplications to perform at once.
/// This number should be optimized to do the most amount of multiplications in a single CPU
const MAX_MULTIPLICATIONS: usize = 16;

/// Factorial and related functions.
///
/// Note that for primitive integer types, it is know the maximum value one can calculate the factorial for before overflowing. See the following table.
///
/// | Type | n! | n!! |
/// |------|----|-----|
/// | u8   | 5  | 7   |
/// | u16  | 8  | 12  |
/// | u32  | 12 | 20  |
/// | u64  | 20 | 33  |
/// | i8   | 5  | 7   |
/// | i16  | 7  | 11  |
/// | i32  | 12 | 19  |
/// | i64  | 20 | 33  |
///
/// The [CheckedFactorial] functions should be used if there is need to care about overflow possibility.
///
/// [CheckedFactorial]: crate::special::CheckedFactorial
pub trait Factorial: Sized {
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
    /// # Notes
    /// If $n < 0$ returns $0$.
    fn factorial(self) -> Self;

    /// The double factorial $n!!$ is defined as the product of all
    /// positive integers up to $n$ that have the same parity as $n$.
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
    /// use sci_rs::special::Factorial;
    ///
    /// assert_eq!(0_isize.factorial2(), 1);
    /// assert_eq!(3_isize.factorial2(), 3);
    /// assert_eq!(6_isize.factorial2(), 48);
    /// ```
    /// # Notes
    /// If $n < 0$ returns $0$.
    fn factorial2(self) -> Self;

    /// Generalized $k$-factorial.
    ///
    /// Generalization of the [factorial] to steps of size $k$.
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
    /// # Notes
    /// If $n < 0$ returns $0$.
    ///
    /// [factorial]: crate::special::Factorial::factorial
    fn factorialk(self, k: Self) -> Self;
}

// Cached factorial and factorial2 values
const FACTORIAL_CACHE_LEN: usize = MAX_MULTIPLICATIONS + 1;
const FACTORIAL_CACHE: [u64; FACTORIAL_CACHE_LEN] = [
    1, 1, 2, 6, 24, 120, 720, 5040, 40320, 362880, 3628800, 39916800, 479001600, 6227020800,
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
            #[inline(always)]
            fn factorial(self) -> Self {
                primint_factorial(self)
            }

            #[inline(always)]
            fn factorial2(self) -> Self {
                primint_factorial2(self)
            }

            #[inline(always)]
            fn factorialk(self, k: Self) -> Self {
                primint_factorialk(self, k)
            }
        }
    )*)
}

factorial_primint_impl! {u8 u16 u32 u64 u128 usize i8 i16 i32 i64 i128 isize}

fn primint_factorial<T>(n: T) -> T
where
    T: PrimInt + FromPrimitive,
{
    if n < T::zero() {
        return T::zero();
    }
    if n < T::from_usize(FACTORIAL_CACHE_LEN).unwrap() {
        return T::from_u64(FACTORIAL_CACHE[n.to_usize().unwrap()]).unwrap();
    }
    partial_product(
        n - T::from_usize(FACTORIAL_CACHE_LEN - 1).unwrap(),
        n,
        T::one(),
    ) * primint_factorial(n - T::from_usize(FACTORIAL_CACHE_LEN).unwrap())
}

fn primint_factorial2<T>(n: T) -> T
where
    T: PrimInt + FromPrimitive,
{
    if n < T::zero() {
        return T::zero();
    }
    if n < T::from_usize(FACTORIAL2_CACHE_LEN).unwrap() {
        return T::from_u64(FACTORIAL2_CACHE[n.to_usize().unwrap()]).unwrap();
    }
    partial_product(
        n - T::from_usize(FACTORIAL2_CACHE_LEN - 1).unwrap(),
        n,
        T::from_usize(2).unwrap(),
    ) * primint_factorial2(n - T::from_usize(FACTORIAL2_CACHE_LEN).unwrap())
}

fn primint_factorialk<T>(n: T, k: T) -> T
where
    T: PrimInt + FromPrimitive,
{
    if n < T::zero() {
        return T::zero();
    }
    assert!(k > T::zero());
    if n.is_zero() {
        return T::one();
    }
    let max_window = k * T::from_usize(MAX_MULTIPLICATIONS).unwrap();

    if n > max_window {
        return partial_product(n - max_window, n, k)
            * primint_factorialk(n - max_window - T::one(), k);
    }
    let window = k * (n / k);
    let window = if window == n { window - k } else { window };
    partial_product(n - window, n, k)
}

/// Checked Factorial functions to prevent overflow.
pub trait CheckedFactorial: Factorial {
    /// Checked version of the [factorial] function which catches overflow.
    ///
    /// # Examples
    /// ```
    /// use sci_rs::special::CheckedFactorial;
    ///
    /// assert_eq!(2_u8.checked_factorial(), Some(2));
    /// assert_eq!(10_u8.checked_factorial(), None); // 10! overflows a u8
    /// ```
    /// [factorial]: crate::special::Factorial::factorial
    fn checked_factorial(self) -> Option<Self>;

    /// Checked version of the [factorial2] function which catches overflow.
    ///
    /// # Examples
    /// ```
    /// use sci_rs::special::CheckedFactorial;
    ///
    /// assert_eq!(2_u8.checked_factorial2(), Some(2));
    /// assert_eq!(5_u8.checked_factorial2(), Some(15));
    /// assert_eq!(33_u8.checked_factorial2(), None); // 33!! overflows a u8
    /// ```
    /// [factorial2]: crate::special::Factorial::factorial2
    fn checked_factorial2(self) -> Option<Self>;

    /// Checked version of the [factorialk] function which catches overflow.
    ///
    /// # Examples
    /// ```
    /// use sci_rs::special::CheckedFactorial;
    ///
    /// assert_eq!(2_u8.checked_factorialk(2), Some(2));
    /// assert_eq!(10_u8.checked_factorialk(4), Some(120));
    /// assert_eq!(31_u8.checked_factorialk(3), None); // 33!! overflows a u8
    /// ```
    /// [factorialk]: crate::special::Factorial::factorialk
    fn checked_factorialk(self, k: Self) -> Option<Self>;
}

macro_rules! impl_slow_checkedfactorial {
    ($($T: ty)*) => ($(
        impl CheckedFactorial for $T {
            #[inline(always)]
            fn checked_factorial(self) -> Option<Self> {
                slow_checked_factorial(self)
            }

            #[inline(always)]
            fn checked_factorial2(self) -> Option<Self> {
                slow_checked_factorial2(self)
            }

            fn checked_factorialk(self, k: Self) -> Option<Self> {
                slow_checked_factorialk(self, k)
            }
        }
    )*)
}

impl_slow_checkedfactorial!(usize isize);

fn slow_checked_factorial<T>(n: T) -> Option<T>
where
    T: PrimInt + FromPrimitive,
{
    if n < T::zero() {
        return Some(T::zero());
    }

    if n < T::from_usize(FACTORIAL_CACHE_LEN)? {
        return T::from_u64(FACTORIAL_CACHE[n.to_usize()?]);
    }
    (checked_partial_product(n - T::from_usize(FACTORIAL_CACHE_LEN - 1)?, n, T::one())?)
        .checked_mul(&slow_checked_factorial(
            n - T::from_usize(FACTORIAL_CACHE_LEN)?,
        )?)
}

fn slow_checked_factorial2<T>(n: T) -> Option<T>
where
    T: PrimInt + FromPrimitive,
{
    if n < T::zero() {
        return Some(T::zero());
    }

    if n < T::from_usize(FACTORIAL2_CACHE_LEN)? {
        return T::from_u64(FACTORIAL2_CACHE[n.to_usize()?]);
    }
    checked_partial_product(
        n - T::from_usize(FACTORIAL2_CACHE_LEN - 1)?,
        n,
        T::from_usize(2).unwrap(),
    )?
    .checked_mul(&slow_checked_factorial2(
        n - T::from_usize(FACTORIAL2_CACHE_LEN)?,
    )?)
}

fn slow_checked_factorialk<T>(n: T, k: T) -> Option<T>
where
    T: PrimInt + FromPrimitive,
{
    if n < T::zero() {
        return Some(T::zero());
    }
    assert!(k > T::zero());
    if n.is_zero() {
        return Some(T::one());
    }

    let max_window = k.checked_mul(&T::from_usize(MAX_MULTIPLICATIONS).unwrap())?;

    if n > max_window {
        return (checked_partial_product(n - max_window, n, k)?)
            .checked_mul(&slow_checked_factorialk(n - max_window - T::one(), k)?);
    }
    let window = k * (n / k);
    let window = if window == n { window - k } else { window };
    checked_partial_product(n - window, n, k)
}

macro_rules! impl_maxfactorial {
    ($($T: ty, $mf: expr, $mf2: expr;)*) => ($(

        impl CheckedFactorial for $T {
            #[inline(always)]
            fn checked_factorial(self) -> Option<Self> {
                if self > $mf {
                    return None;
                }
                Some(self.factorial())
            }

            #[inline(always)]
            fn checked_factorial2(self) -> Option<Self> {
                if self > $mf2 {
                    return None;
                }
                Some(self.factorial2())
            }

            #[inline(always)]
            fn checked_factorialk(self, k: Self) -> Option<Self> {
                slow_checked_factorialk(self, k)
            }
        }
    )*)
}

impl_maxfactorial!(u8, 5, 7; u16, 8, 12; u32, 12, 20; u64, 20, 33; u128, 34, 56;);
impl_maxfactorial!(i8, 5, 7; i16, 7, 11; i32, 12, 19; i64, 20, 33; i128, 33, 55;);

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
    T: PrimInt,
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
fn partial_product<T>(mut start: T, stop: T, step: T) -> T
where
    T: PrimInt,
{
    let mut result = T::one();
    while start <= stop {
        result = result * start;
        start = start + step;
    }
    result
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
    fn factorial_u16() {
        assert_eq!(6_u16.factorial(), 720);
        assert_eq!(7_u16.factorial(), 5040);
        assert_eq!(8_u16.factorial(), 40320);
        assert_eq!(9_u16.checked_factorial(), None);
    }

    #[test]
    fn factorial_u32() {
        assert_eq!(9_u32.factorial(), 362880);
        assert_eq!(10_u32.factorial(), 3628800);
        assert_eq!(11_u32.factorial(), 39916800);
        assert_eq!(12_u32.factorial(), 479001600);
        assert_eq!(13_u32.checked_factorial(), None);
    }

    #[test]
    fn factorial_u64() {
        assert_eq!(13_u64.factorial(), 6227020800);
        assert_eq!(14_u64.factorial(), 87178291200);
        assert_eq!(15_u64.factorial(), 1307674368000);
        assert_eq!(16_u64.factorial(), 20922789888000);
        assert_eq!(17_u64.factorial(), 355687428096000);
        assert_eq!(18_u64.factorial(), 6402373705728000);
        assert_eq!(19_u64.factorial(), 121645100408832000);
        assert_eq!(20_u64.factorial(), 2432902008176640000);
        assert_eq!(21_u64.checked_factorial(), None);
    }

    #[test]
    fn factorial_u128() {
        assert_eq!(21_u128.factorial(), 51090942171709440000);
        assert_eq!(22_u128.factorial(), 1124000727777607680000);
        assert_eq!(23_u128.factorial(), 25852016738884976640000);
        assert_eq!(24_u128.factorial(), 620448401733239439360000);
        assert_eq!(25_u128.factorial(), 15511210043330985984000000);
        assert_eq!(26_u128.factorial(), 403291461126605635584000000);
        assert_eq!(27_u128.factorial(), 10888869450418352160768000000);
        assert_eq!(28_u128.factorial(), 304888344611713860501504000000);
        assert_eq!(29_u128.factorial(), 8841761993739701954543616000000);
        assert_eq!(30_u128.factorial(), 265252859812191058636308480000000);
        assert_eq!(31_u128.factorial(), 8222838654177922817725562880000000);
        assert_eq!(32_u128.factorial(), 263130836933693530167218012160000000);
        assert_eq!(33_u128.factorial(), 8683317618811886495518194401280000000);
        assert_eq!(34_u128.factorial(), 295232799039604140847618609643520000000);
        assert_eq!(35_u128.checked_factorial(), None);
    }

    #[test]
    fn factorial_i8() {
        assert_eq!(0_i8.factorial(), 1);
        assert_eq!(1_i8.factorial(), 1);
        assert_eq!(2_i8.factorial(), 2);
        assert_eq!(3_i8.factorial(), 6);
        assert_eq!(4_i8.factorial(), 24);
        assert_eq!(5_i8.factorial(), 120);
        assert_eq!(6_i8.checked_factorial(), None);
    }

    #[test]
    fn factorial_i16() {
        assert_eq!(6_i16.factorial(), 720);
        assert_eq!(7_i16.factorial(), 5040);
        assert_eq!(8_i16.checked_factorial(), None);
    }

    #[test]
    fn factorial_i32() {
        assert_eq!(8_i32.factorial(), 40320);
        assert_eq!(9_i32.factorial(), 362880);
        assert_eq!(10_i32.factorial(), 3628800);
        assert_eq!(11_i32.factorial(), 39916800);
        assert_eq!(12_i32.factorial(), 479001600);
        assert_eq!(13_i32.checked_factorial(), None);
    }

    #[test]
    fn factorial_i64() {
        assert_eq!(13_i64.factorial(), 6227020800);
        assert_eq!(14_i64.factorial(), 87178291200);
        assert_eq!(15_i64.factorial(), 1307674368000);
        assert_eq!(16_i64.factorial(), 20922789888000);
        assert_eq!(17_i64.factorial(), 355687428096000);
        assert_eq!(18_i64.factorial(), 6402373705728000);
        assert_eq!(19_i64.factorial(), 121645100408832000);
        assert_eq!(20_i64.factorial(), 2432902008176640000);
        assert_eq!(21_i64.checked_factorial(), None);
    }

    #[test]
    fn factorial_i128() {
        assert_eq!(21_i128.factorial(), 51090942171709440000);
        assert_eq!(22_i128.factorial(), 1124000727777607680000);
        assert_eq!(23_i128.factorial(), 25852016738884976640000);
        assert_eq!(24_i128.factorial(), 620448401733239439360000);
        assert_eq!(25_i128.factorial(), 15511210043330985984000000);
        assert_eq!(26_i128.factorial(), 403291461126605635584000000);
        assert_eq!(27_i128.factorial(), 10888869450418352160768000000);
        assert_eq!(28_i128.factorial(), 304888344611713860501504000000);
        assert_eq!(29_i128.factorial(), 8841761993739701954543616000000);
        assert_eq!(30_i128.factorial(), 265252859812191058636308480000000);
        assert_eq!(31_i128.factorial(), 8222838654177922817725562880000000);
        assert_eq!(32_i128.factorial(), 263130836933693530167218012160000000);
        assert_eq!(33_i128.factorial(), 8683317618811886495518194401280000000);
        assert_eq!(34_i128.checked_factorial(), None);
    }

    #[test]
    fn factorial2_u8() {
        assert_eq!(0_u8.factorial2(), 1);
        assert_eq!(1_u8.factorial2(), 1);
        assert_eq!(2_u8.factorial2(), 2);
        assert_eq!(3_u8.factorial2(), 3);
        assert_eq!(4_u8.factorial2(), 8);
        assert_eq!(5_u8.factorial2(), 15);
        assert_eq!(6_u8.factorial2(), 48);
        assert_eq!(7_u8.factorial2(), 105);
        assert_eq!(8_u8.checked_factorial2(), None);
    }

    #[test]
    fn factorial2_u16() {
        assert_eq!(8_u16.factorial2(), 384);
        assert_eq!(9_u16.factorial2(), 945);
        assert_eq!(10_u16.factorial2(), 3840);
        assert_eq!(11_u16.factorial2(), 10395);
        assert_eq!(12_u16.factorial2(), 46080);
        assert_eq!(13_u16.checked_factorial2(), None);
    }

    #[test]
    fn factorial2_u32() {
        assert_eq!(13_u32.factorial2(), 135135);
        assert_eq!(14_u32.factorial2(), 645120);
        assert_eq!(15_u32.factorial2(), 2027025);
        assert_eq!(16_u32.factorial2(), 10321920);
        assert_eq!(17_u32.factorial2(), 34459425);
        assert_eq!(18_u32.factorial2(), 185794560);
        assert_eq!(19_u32.factorial2(), 654729075);
        assert_eq!(20_u32.factorial2(), 3715891200);
        assert_eq!(21_u32.checked_factorial2(), None);
    }

    #[test]
    fn factorial2_u64() {
        assert_eq!(21_u64.factorial2(), 13749310575);
        assert_eq!(22_u64.factorial2(), 81749606400);
        assert_eq!(23_u64.factorial2(), 316234143225);
        assert_eq!(24_u64.factorial2(), 1961990553600);
        assert_eq!(25_u64.factorial2(), 7905853580625);
        assert_eq!(26_u64.factorial2(), 51011754393600);
        assert_eq!(27_u64.factorial2(), 213458046676875);
        assert_eq!(28_u64.factorial2(), 1428329123020800);
        assert_eq!(29_u64.factorial2(), 6190283353629375);
        assert_eq!(30_u64.factorial2(), 42849873690624000);
        assert_eq!(31_u64.factorial2(), 191898783962510625);
        assert_eq!(32_u64.factorial2(), 1371195958099968000);
        assert_eq!(33_u64.factorial2(), 6332659870762850625);
        assert_eq!(34_u64.checked_factorial2(), None);
    }

    #[test]
    fn factorial2_u128() {
        assert_eq!(21_u128.factorial2(), 13749310575);
        assert_eq!(22_u128.factorial2(), 81749606400);
        assert_eq!(23_u128.factorial2(), 316234143225);
        assert_eq!(24_u128.factorial2(), 1961990553600);
        assert_eq!(25_u128.factorial2(), 7905853580625);
        assert_eq!(26_u128.factorial2(), 51011754393600);
        assert_eq!(27_u128.factorial2(), 213458046676875);
        assert_eq!(28_u128.factorial2(), 1428329123020800);
        assert_eq!(29_u128.factorial2(), 6190283353629375);
        assert_eq!(30_u128.factorial2(), 42849873690624000);
        assert_eq!(31_u128.factorial2(), 191898783962510625);
        assert_eq!(32_u128.factorial2(), 1371195958099968000);
        assert_eq!(33_u128.factorial2(), 6332659870762850625);
        assert_eq!(34_u128.factorial2(), 46620662575398912000);
        assert_eq!(35_u128.factorial2(), 443286190953399543750);
        assert_eq!(36_u128.factorial2(), 2517515779071541248000);
        assert_eq!(37_u128.factorial2(), 21868785420367710825000);
        assert_eq!(38_u128.factorial2(), 119581999505898209280000);
        assert_eq!(39_u128.factorial2(), 1023459157673208866610000);
        assert_eq!(40_u128.factorial2(), 5580493310275249766400000);
        assert_eq!(41_u128.factorial2(), 47956371959544644035440000);
        assert_eq!(42_u128.factorial2(), 263678308910505551462400000);
        assert_eq!(43_u128.factorial2(), 2291248882511577437248800000);
        assert_eq!(44_u128.factorial2(), 12762030151268468690780160000);
        assert_eq!(45_u128.factorial2(), 112479490596022892374032000000);
        assert_eq!(46_u128.factorial2(), 635974502538212023090544640000);
        assert_eq!(47_u128.factorial2(), 5693192677860235629393312000000);
        assert_eq!(48_u128.factorial2(), 32707260130536618330370867200000);
        assert_eq!(49_u128.factorial2(), 297564203962828315562957107200000);
        assert_eq!(50_u128.factorial2(), 1737573194434757848800952320000000);
        assert_eq!(51_u128.factorial2(), 16068467013992729040399683788800000);
        assert_eq!(52_u128.factorial2(), 95373462005641153034185605120000000);
        assert_eq!(53_u128.factorial2(), 896451317622752251727561306112000000);
        assert_eq!(54_u128.factorial2(), 5407675295719853377038323810304000000);
        assert_eq!(55_u128.factorial2(), 51652671158263344028111865733120000000);
        assert_eq!(
            56_u128.factorial2(),
            316594808222144143164789139439616000000
        );
        assert_eq!(57_u128.checked_factorial2(), None);
    }

    #[test]
    fn factorial2_i8() {
        assert_eq!(0_i8.factorial2(), 1);
        assert_eq!(1_i8.factorial2(), 1);
        assert_eq!(2_i8.factorial2(), 2);
        assert_eq!(3_i8.factorial2(), 3);
        assert_eq!(4_i8.factorial2(), 8);
        assert_eq!(5_i8.factorial2(), 15);
        assert_eq!(6_i8.factorial2(), 48);
        assert_eq!(7_i8.factorial2(), 105);
        assert_eq!(8_i8.checked_factorial2(), None);
    }

    #[test]
    fn factorial2_i16() {
        assert_eq!(8_i16.factorial2(), 384);
        assert_eq!(9_i16.factorial2(), 945);
        assert_eq!(10_i16.factorial2(), 3840);
        assert_eq!(11_i16.factorial2(), 10395);
        assert_eq!(12_i16.checked_factorial2(), None);
    }

    #[test]
    fn factorial2_i32() {
        assert_eq!(12_i32.factorial2(), 46080);
        assert_eq!(13_i32.factorial2(), 135135);
        assert_eq!(14_i32.factorial2(), 645120);
        assert_eq!(15_i32.factorial2(), 2027025);
        assert_eq!(16_i32.factorial2(), 10321920);
        assert_eq!(17_i32.factorial2(), 34459425);
        assert_eq!(18_i32.factorial2(), 185794560);
        assert_eq!(19_i32.factorial2(), 654729075);
        assert_eq!(20_i32.checked_factorial2(), None);
    }

    #[test]
    fn factorial2_i64() {
        assert_eq!(20_i64.factorial2(), 3715891200);
        assert_eq!(21_i64.factorial2(), 13749310575);
        assert_eq!(22_i64.factorial2(), 81749606400);
        assert_eq!(23_i64.factorial2(), 316234143225);
        assert_eq!(24_i64.factorial2(), 1961990553600);
        assert_eq!(25_i64.factorial2(), 7905853580625);
        assert_eq!(26_i64.factorial2(), 51011754393600);
        assert_eq!(27_i64.factorial2(), 213458046676875);
        assert_eq!(28_i64.factorial2(), 1428329123020800);
        assert_eq!(29_i64.factorial2(), 6190283353629375);
        assert_eq!(30_i64.factorial2(), 42849873690624000);
        assert_eq!(31_i64.factorial2(), 191898783962510625);
        assert_eq!(32_i64.factorial2(), 1371195958099968000);
        assert_eq!(33_i64.factorial2(), 6332659870762850625);
        assert_eq!(34_i64.checked_factorial2(), None);
    }

    #[test]
    fn factorial2_i128() {
        assert_eq!(34_i128.factorial2(), 46620662575398912000);
        assert_eq!(35_i128.factorial2(), 443286190953399543750);
        assert_eq!(36_i128.factorial2(), 2517515779071541248000);
        assert_eq!(37_i128.factorial2(), 21868785420367710825000);
        assert_eq!(38_i128.factorial2(), 119581999505898209280000);
        assert_eq!(39_i128.factorial2(), 1023459157673208866610000);
        assert_eq!(40_i128.factorial2(), 5580493310275249766400000);
        assert_eq!(41_i128.factorial2(), 47956371959544644035440000);
        assert_eq!(42_i128.factorial2(), 263678308910505551462400000);
        assert_eq!(43_i128.factorial2(), 2291248882511577437248800000);
        assert_eq!(44_i128.factorial2(), 12762030151268468690780160000);
        assert_eq!(45_i128.factorial2(), 112479490596022892374032000000);
        assert_eq!(46_i128.factorial2(), 635974502538212023090544640000);
        assert_eq!(47_i128.factorial2(), 5693192677860235629393312000000);
        assert_eq!(48_i128.factorial2(), 32707260130536618330370867200000);
        assert_eq!(49_i128.factorial2(), 297564203962828315562957107200000);
        assert_eq!(50_i128.factorial2(), 1737573194434757848800952320000000);
        assert_eq!(51_i128.factorial2(), 16068467013992729040399683788800000);
        assert_eq!(52_i128.factorial2(), 95373462005641153034185605120000000);
        assert_eq!(53_i128.factorial2(), 896451317622752251727561306112000000);
        assert_eq!(54_i128.factorial2(), 5407675295719853377038323810304000000);
        assert_eq!(55_i128.factorial2(), 51652671158263344028111865733120000000);
        assert_eq!(56_i128.checked_factorial2(), None);
    }

    #[test]
    fn factorial1() {
        let k = 1_usize;
        for i in 0..10 {
            assert_eq!(i.factorialk(k), i.factorial());
        }
    }

    #[test]
    fn factorial2() {
        let k = 2_usize;
        for i in 0..15 {
            assert_eq!(i.factorialk(k), i.factorial2());
        }
    }

    #[test]
    fn factorial3() {
        let k = 3_usize;
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
