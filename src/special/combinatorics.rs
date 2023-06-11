use crate::special::IsNegative;
use num_traits::{CheckedAdd, CheckedMul};
use std::cmp::min;

/// Defined the `comb` function for `Self`. Implemented for primitive
/// integer types.
pub trait Combinatorics: Sized + CheckedMul + CheckedAdd + Copy {
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
    fn choose(&self, k: Self) -> Self {
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
    fn checked_choose(&self, k: Self) -> Option<Self>;
}

macro_rules! combinatorics_primint_impl {
    ($($T: ty)*) => ($(
        impl Combinatorics for $T {
            fn choose(&self, k: Self) -> Self {
                if k > *self || self.is_negative() {
                    return 0;
                }
                let m = *self + 1;
                let n_terms = min(k, *self-k) + 1;

                (1..n_terms).fold(1, |result, i| (result * (m-i))/i)
            }

            fn checked_choose(&self, k: Self) -> Option<Self> {
                if k > *self || self.is_negative() {
                    return Some(0);
                }

                let m = (*self).checked_add(1)?;
                let n_terms = min(k, *self-k) + 1;

                let mut result = 1;
                for i in 1..n_terms {
                    result = result.checked_mul(&(m-i))? / i;
                }
                Some(result)
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
    fn comb() {
        assert_eq!(3_u8.choose(1), 3);
        assert_eq!(3_u8.choose(2), 3);
        assert_eq!(3_u8.choose(3), 1);
        assert_eq!(3_u8.checked_choose(1), Some(3_u8.choose(1)));
        assert_eq!(3_u8.checked_choose(2), Some(3_u8.choose(2)));
        assert_eq!(3_u8.checked_choose(3), Some(3_u8.choose(3)));
        assert_eq!(10_u8.checked_choose(5), None);

        fn check_values<T>(x: T, ref_values: &[T])
        where
            T: Combinatorics + std::cmp::PartialEq + std::fmt::Debug + FromPrimitive,
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
    fn comb_negatives() {
        for i in 0..4 {
            assert_eq!((-4).choose(i), 0);
            assert_eq!((-3).choose(i), 0);
            assert_eq!((-3241).choose(i), 0);
        }
    }
}
