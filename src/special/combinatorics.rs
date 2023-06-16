//**********************************************************************
// This file is part of sci-rs                                         *
// Copyright 2023 Matthew R. Hennefarth                                *
//**********************************************************************

use num_traits::{FromPrimitive, PrimInt};
use std::cmp::min;

/// Various combinatorics functions for integer-types.
pub trait Comb {
    /// The number of combinations of $n$ taken $k$ at a time.
    ///
    /// This is also known as $n$ choose $k$ and is generally given by
    /// the formula
    /// $$
    /// \begin{pmatrix}
    /// n \\\\ k
    /// \end{pmatrix} = \frac{n!}{k!(n-k)!}
    /// $$
    /// # Examples
    /// ```
    /// use sci_rs::special::Comb;
    /// assert_eq!(5.choose(2), 10);
    /// assert_eq!(7.choose(3), 35);
    /// ```
    /// # Notes
    /// When $n$ < 0 or $k<0$ or $n<k$, the $0$ is returned.
    ///
    /// ## Implementation
    /// Does not actually compute using the factorial functions as this would likely lead to
    /// unnecessary overflows. Instead a different approach is taken which uses an iterative
    /// multiplicative formula.
    fn choose(self, k: Self) -> Self;

    /// Number of combinations with repetition.
    ///
    /// This is also known as a $k$-combination with repetitions or
    /// $k$-multicombinations. Formula is given as
    /// $$
    /// \left(\begin{pmatrix}
    ///    n \\\\ k
    /// \end{pmatrix}\right) = \begin{pmatrix} n + k -1 \\\\ k \end{pmatrix}
    /// $$
    /// For more detailed explanation see the [wiki] page.
    ///
    /// # Examples
    /// ```
    /// use sci_rs::special::Comb;
    /// assert_eq!(5.choose_rep(2), 15);
    /// assert_eq!(10.choose_rep(3), 220);
    /// ```
    ///
    /// # Notes
    /// When $n$ < 0 or $k<0$ or $n<k$, then $0$ is returned.
    ///
    /// # References
    /// - [wiki]
    ///
    /// [wiki]: https://en.wikipedia.org/wiki/Combination#Number_of_combinations_with_repetition
    fn choose_rep(self, k: Self) -> Self;

    /// Number of permutations of $n$ things taken $k$ at a time.
    ///
    /// Also known as the $k$-permutations of $n$.
    /// $$
    /// \text{Perm}(n, k) = \frac{n!}{(n-k)!}
    /// $$
    /// # Examples
    /// ```
    /// use sci_rs::special::Comb;
    /// assert_eq!(5.perm(5), 120); // should be 5!
    /// assert_eq!(5.perm(0), 1);
    /// assert_eq!(6.perm(3), 6*5*4);
    /// ```
    ///
    /// # Notes
    /// When $n$ < 0 or $k<0$, then $0$ is returned.
    fn perm(self, k: Self) -> Self;
}

macro_rules! comb_primint_impl {
    ($($T: ty)*) => ($(
        impl Comb for $T {
            #[inline(always)]
            fn choose(self, k: Self) -> Self {
                primint_choose(self, k)
            }

            #[inline(always)]
            fn choose_rep(self, k: Self) -> Self {
                primint_choose_rep(self, k)
            }

            #[inline(always)]
            fn perm(self, k: Self) -> Self {
                primint_perm(self, k)
            }
        }
    )*)
}

comb_primint_impl! {u8 u16 u32 u64 usize i8 i16 i32 i64 isize}
#[cfg(has_u128)]
comb_primint_impl! {u128}
#[cfg(has_i128)]
comb_primint_impl! {i128}

fn primint_choose<T>(n: T, k: T) -> T
where
    T: PrimInt + FromPrimitive,
{
    if k > n || n < T::zero() || k < T::zero() {
        return T::zero();
    }
    let m = n + T::one();
    let n_terms = (min(k, n - k) + T::one()).to_usize().unwrap();
    (1..n_terms).fold(T::one(), |result, i| {
        result * (m - T::from_usize(i).unwrap()) / T::from_usize(i).unwrap()
    })
}

fn primint_choose_rep<T>(n: T, k: T) -> T
where
    T: PrimInt + FromPrimitive,
{
    primint_choose(n + k - T::one(), k)
}

fn primint_perm<T>(n: T, k: T) -> T
where
    T: PrimInt + FromPrimitive,
{
    if k > n || n < T::zero() || k < T::zero() {
        return T::zero();
    }

    let start = (n - k + T::one()).to_usize().unwrap();
    let end = (n + T::one()).to_usize().unwrap();
    (start..end).fold(T::one(), |result, val| result * T::from_usize(val).unwrap())
}

/// Checked combinatorics functions.
pub trait CheckedComb: Sized {
    /// Checked version of the [choose] function.
    ///
    /// # Examples
    /// ```
    /// use sci_rs::special::CheckedComb;
    /// assert_eq!(5.checked_choose(2), Some(10));
    /// assert_eq!(7.checked_choose(3), Some(35));
    /// assert_eq!(50_u8.checked_choose(4), None); // Overflows a u8
    /// ```
    /// [choose]: crate::special::Comb::choose
    fn checked_choose(self, k: Self) -> Option<Self>;

    /// Checked version of [choose_rep] to prevent overflow panics.
    ///
    /// # Examples
    ///```
    /// use sci_rs::special::CheckedComb;
    /// assert_eq!(5.checked_choose_rep(2), Some(15));
    /// assert_eq!(10.checked_choose_rep(3), Some(220));
    /// assert_eq!(12_u8.checked_choose_rep(3), None); // Overflows a u8
    /// ```
    /// [choose_rep]: crate::special::Comb::choose_rep
    fn checked_choose_rep(self, k: Self) -> Option<Self>;

    /// Checked version of [perm] to prevent overflow panics.
    ///
    /// # Examples
    ///```
    /// use sci_rs::special::CheckedComb;
    /// assert_eq!(154.checked_perm(154), None); //Should overflow since 154!
    /// assert_eq!(4.checked_perm(3), Some(4*3*2*1));
    /// ```
    /// [perm]: crate::special::Comb::perm
    fn checked_perm(self, k: Self) -> Option<Self>;
}

macro_rules! checkedcomb_primint_impl {
    ($($T: ty)*) => ($(
        impl CheckedComb for $T {
            #[inline(always)]
            fn checked_choose(self, k: Self) -> Option<Self> {
                primint_checked_choose(self, k)
            }

            #[inline(always)]
            fn checked_choose_rep(self, k: Self) -> Option<Self> {
                primint_checked_choose_rep(self, k)
            }

            #[inline(always)]
            fn checked_perm(self, k: Self) -> Option<Self> {
                primint_checked_perm(self, k)
            }
        }
    )*)
}

checkedcomb_primint_impl! {u8 u16 u32 u64 usize i8 i16 i32 i64 isize}
#[cfg(has_u128)]
checkedcomb_primint_impl! {u128}
#[cfg(has_i128)]
checkedcomb_primint_impl! {i128}

fn primint_checked_choose<T>(n: T, k: T) -> Option<T>
where
    T: PrimInt + FromPrimitive,
{
    if k > n || n < T::zero() || k < T::zero() {
        return Some(T::zero());
    }

    let m = n.checked_add(&T::one())?;
    let n_terms = (min(k, n - k) + T::one()).to_usize()?;

    let mut result = T::one();
    for i in 1..n_terms {
        let i = T::from_usize(i)?;
        result = result.checked_mul(&(m - i))? / i;
    }
    Some(result)
}

fn primint_checked_choose_rep<T>(n: T, k: T) -> Option<T>
where
    T: PrimInt + FromPrimitive,
{
    primint_checked_choose(n.checked_add(&k)?.checked_sub(&T::one())?, k)
}

fn primint_checked_perm<T>(n: T, k: T) -> Option<T>
where
    T: PrimInt + FromPrimitive,
{
    if k > n || n < T::zero() || k < T::zero() {
        return Some(T::zero());
    }

    let start = (n - k).checked_add(&T::one())?.to_usize()?;
    let end = n.checked_add(&T::one())?.to_usize()?;

    let mut result = T::one();
    for i in start..end {
        result = result.checked_mul(&T::from_usize(i)?)?;
    }
    Some(result)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn check_values<T>(
        x: T,
        ref_values: &[T],
        func: fn(T, T) -> T,
        checked_func: fn(T, T) -> Option<T>,
    ) where
        T: PrimInt + FromPrimitive + std::fmt::Debug,
    {
        for (i, &val) in ref_values.iter().enumerate() {
            let i = T::from_usize(i).unwrap();
            assert_eq!(func(x, i), val);
            assert_eq!(checked_func(x, i), Some(val));
        }
    }

    #[test]
    fn choose() {
        assert_eq!(3_u8.choose(1), 3);
        assert_eq!(3_u8.choose(2), 3);
        assert_eq!(3_u8.choose(3), 1);
        assert_eq!(3_u8.checked_choose(1), Some(3_u8.choose(1)));
        assert_eq!(3_u8.checked_choose(2), Some(3_u8.choose(2)));
        assert_eq!(3_u8.checked_choose(3), Some(3_u8.choose(3)));
        assert_eq!(10_u8.checked_choose(5), None);

        let ref_values_5 = [1, 5, 10, 10, 5, 1, 0];
        let ref_values_10 = [1, 10, 45, 120, 210, 252, 210, 120, 45, 10, 1, 0];
        let ref_values_15 = [
            1, 15, 105, 455, 1365, 3003, 5005, 6435, 6435, 5005, 3003, 1365, 455, 105, 15, 1,
        ];
        check_values(5, &ref_values_5, i32::choose, i32::checked_choose);
        check_values(10, &ref_values_10, i32::choose, i32::checked_choose);
        check_values(15, &ref_values_15, i32::choose, i32::checked_choose);
    }

    #[test]
    fn choose_negatives() {
        for i in 0..4 {
            assert_eq!((-4).choose(i), 0);
            assert_eq!((-3).choose(i), 0);
            assert_eq!((-3241).choose(i), 0);
        }

        for i in -4..0 {
            assert_eq!(4.choose(i), 0);
            assert_eq!(2.choose(i), 0);
            assert_eq!(2341.choose(i), 0);
            assert_eq!((-2).choose(i), 0);
            assert_eq!((-4).choose(i), 0);
            assert_eq!((-5).choose(i), 0);
            assert_eq!((-3241).choose(i), 0);
        }
    }

    #[test]
    fn choose_edge() {
        assert_eq!(0.choose(0), 1);
    }

    #[test]
    fn choose_repl() {
        assert_eq!(0.choose_rep(0), 0);
        assert_eq!(1.choose_rep(0), 1);
        assert_eq!(0.choose_rep(1), 0);
        assert_eq!(1.choose_rep(1), 1);

        let ref_values_5 = [1, 5, 15, 35, 70, 126, 210, 330, 495, 715];
        let ref_values_7 = [1, 7, 28, 84, 210, 462, 924, 1716, 3003, 5005];
        let ref_values_10 = [
            1, 10, 55, 220, 715, 2002, 5005, 11440, 24310, 48620, 92378, 167960, 293930, 497420,
            817190,
        ];
        check_values(5, &ref_values_5, i32::choose_rep, i32::checked_choose_rep);
        check_values(7, &ref_values_7, i32::choose_rep, i32::checked_choose_rep);
        check_values(10, &ref_values_10, i32::choose_rep, i32::checked_choose_rep);
    }

    #[test]
    fn perm() {
        let ref_values_4 = [1, 4, 12, 24, 24, 0];
        let ref_values_7 = [1, 7, 42, 210, 840, 2520, 5040, 5040, 0];
        let ref_values_13 = [
            1, 13, 156, 1716, 17160, 154440, 1235520, 8648640, 51891840, 259459200, 1037836800,
        ];
        //3113510400, 6227020800, 6227020800,
        check_values(4, &ref_values_4, i32::perm, i32::checked_perm);
        check_values(7, &ref_values_7, i32::perm, i32::checked_perm);
        check_values(13, &ref_values_13, i32::perm, i32::checked_perm);
    }

    #[test]
    fn perm_edge() {
        assert_eq!(0.perm(0), 1);
        assert_eq!(1.perm(0), 1);
        assert_eq!(0.perm(1), 0);
    }

    #[test]
    fn perm_negative() {
        for i in 0..4 {
            assert_eq!((-4).perm(i), 0);
            assert_eq!((-3).perm(i), 0);
            assert_eq!((-3241).perm(i), 0);
        }

        for i in -4..0 {
            assert_eq!(4.perm(i), 0);
            assert_eq!(2.perm(i), 0);
            assert_eq!(2341.perm(i), 0);
            assert_eq!((-2).perm(i), 0);
            assert_eq!((-4).perm(i), 0);
            assert_eq!((-5).perm(i), 0);
            assert_eq!((-3241).perm(i), 0);
        }
    }
}
