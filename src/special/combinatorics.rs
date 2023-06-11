use crate::special::IsNegative;
use num_traits::{CheckedAdd, CheckedMul};
use std::cmp::min;

/// Defined the `comb` function for `Self`. Implemented for primitive
/// integer types.
pub trait Combinatorics: Sized + CheckedMul + CheckedAdd {
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
    /// use sci_rs::special::Combinatorics;
    /// assert_eq!(5.choose(2), 10);
    /// assert_eq!(7.choose(3), 35);
    /// ```
    /// ## Notes
    /// When $n$ < 0 or $k<0$ or $n<k$, the $0$ is returned.
    ///
    /// ## Implementation
    /// Does not actually compute using the factorial functions as this would likely lead to
    /// unnecessary overflows. Instead a different approach is taken which uses an iterative
    /// multiplicative formula.
    fn choose(self, k: Self) -> Self {
        self.checked_choose(k).unwrap()
    }

    /// Checked version of the `choose` function.
    ///
    /// # Examples
    /// ```
    /// use sci_rs::special::Combinatorics;
    /// assert_eq!(5.checked_choose(2), Some(10));
    /// assert_eq!(7.checked_choose(3), Some(35));
    /// assert_eq!(50_u8.checked_choose(4), None); // Overflows a u8
    /// ```
    fn checked_choose(self, k: Self) -> Option<Self>;

    /// Number of combinations with repetition.
    ///
    /// This is also known as a $k$-combination with repetitions or
    /// $k$-multicombinations. Formula is given as
    /// $$
    /// \left(\begin{pmatrix}
    ///    n \\\\ k
    /// \end{pmatrix}\right) = \begin{pmatrix} n + k -1 \\\\ k \end{pmatrix}
    /// $$
    /// For more detailed explanation see
    /// [Wiki](https://en.wikipedia.org/wiki/Combination#Number_of_combinations_with_repetition)
    ///
    /// # Examples
    /// ```
    /// use sci_rs::special::Combinatorics;
    /// assert_eq!(5.choose_rep(2), 15);
    /// assert_eq!(10.choose_rep(3), 220);
    /// ```
    ///
    /// ## Notes
    /// When $n$ < 0 or $k<0$ or $n<k$, the $0$ is returned.
    fn choose_rep(self, k: Self) -> Self {
        self.checked_choose_rep(k).unwrap()
    }

    /// Checked version of `choose_rep` to prevent overflow panics.
    ///
    /// # Examples
    ///```
    /// use sci_rs::special::Combinatorics;
    /// assert_eq!(5.checked_choose_rep(2), Some(15));
    /// assert_eq!(10.checked_choose_rep(3), Some(220));
    /// assert_eq!(12_u8.checked_choose_rep(3), None); // Overflows a u8
    /// ```
    fn checked_choose_rep(self, k: Self) -> Option<Self>;
}

macro_rules! combinatorics_primint_impl {
    ($($T: ty)*) => ($(
        impl Combinatorics for $T {
            fn choose(self, k: Self) -> Self {
                if k > self || self.is_negative() || k.is_negative() {
                    return 0;
                }
                let m = self + 1;
                let n_terms = min(k, self-k) + 1;

                (1..n_terms).fold(1, |result, i| (result * (m-i))/i)
            }

            fn checked_choose(self, k: Self) -> Option<Self> {
                if k > self || self.is_negative() || k.is_negative()  {
                    return Some(0);
                }

                let m = self.checked_add(1)?;
                let n_terms = min(k, self-k) + 1;

                let mut result = 1;
                for i in 1..n_terms {
                    result = result.checked_mul(&(m-i))? / i;
                }
                Some(result)
            }

            fn choose_rep(self, k: Self) -> Self {
                (self + k - 1).choose(k)
            }

            fn checked_choose_rep(self, k: Self) -> Option<Self> {
                self.checked_add(k)?.checked_sub(1)?.checked_choose(k)
            }
        }
    )*)
}

combinatorics_primint_impl! {u8 u16 u32 u64 usize i8 i16 i32 i64 isize}
#[cfg(has_i128)]
factorial_primint_impl! {u128 i128}

#[cfg(test)]
mod tests {
    use super::*;
    use num_traits::FromPrimitive;

    #[test]
    fn choose() {
        assert_eq!(3_u8.choose(1), 3);
        assert_eq!(3_u8.choose(2), 3);
        assert_eq!(3_u8.choose(3), 1);
        assert_eq!(3_u8.checked_choose(1), Some(3_u8.choose(1)));
        assert_eq!(3_u8.checked_choose(2), Some(3_u8.choose(2)));
        assert_eq!(3_u8.checked_choose(3), Some(3_u8.choose(3)));
        assert_eq!(10_u8.checked_choose(5), None);

        fn check_values<T>(x: T, ref_values: &[T])
        where
            T: Combinatorics + PartialEq + std::fmt::Debug + FromPrimitive + Copy,
        {
            for (i, &val) in ref_values.iter().enumerate() {
                let i = T::from_usize(i).unwrap();
                assert_eq!(x.choose(i), val);
                assert_eq!(x.checked_choose(i as T), Some(val));
            }
        }

        let ref_values_5 = [1, 5, 10, 10, 5, 1, 0];
        let ref_values_10 = [1, 10, 45, 120, 210, 252, 210, 120, 45, 10, 1, 0];
        let ref_values_15 = [
            1, 15, 105, 455, 1365, 3003, 5005, 6435, 6435, 5005, 3003, 1365, 455, 105, 15, 1,
        ];
        check_values(5, &ref_values_5);
        check_values(10, &ref_values_10);
        check_values(15, &ref_values_15);
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

        fn check_values<T>(x: T, ref_values: &[T])
        where
            T: Combinatorics + PartialEq + std::fmt::Debug + FromPrimitive + Copy,
        {
            for (i, &val) in ref_values.iter().enumerate() {
                let i = T::from_usize(i).unwrap();
                assert_eq!(x.choose_rep(i), val);
                assert_eq!(x.checked_choose_rep(i as T), Some(val));
            }
        }

        let ref_values_5 = [1, 5, 15, 35, 70, 126, 210, 330, 495, 715];
        let ref_values_7 = [1, 7, 28, 84, 210, 462, 924, 1716, 3003, 5005];
        let ref_values_10 = [
            1, 10, 55, 220, 715, 2002, 5005, 11440, 24310, 48620, 92378, 167960, 293930, 497420,
            817190,
        ];
        check_values(5, &ref_values_5);
        check_values(7, &ref_values_7);
        check_values(10, &ref_values_10);
    }
}
