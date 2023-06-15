//**********************************************************************
// This file is part of sci-rs                                         *
// Copyright 2023 Matthew R. Hennefarth                                *
//**********************************************************************

use crate::special::Combinatorics;

/// Returns first $n$ Bernoulli numbers using recurssion.
///
/// bernoulli_rec(3) -> [0.0, -0.5, 1/6, 0.0] (list of 4)
pub(crate) fn bernoulli_rec(n: usize) -> Vec<f64> {
    let result_len = n + 1;
    const CACHE_SIZE: usize = 4;
    let mut result = vec![1.0, -0.5, 1.0 / 6.0, 0.0];
    debug_assert_eq!(result.len(), CACHE_SIZE);

    if result_len <= CACHE_SIZE {
        return result[..result_len].to_vec();
    }

    result.reserve_exact(result_len - CACHE_SIZE);

    for i in CACHE_SIZE..result_len {
        if i % 2 == 1 {
            result.push(0.0);
        } else {
            let np1 = i + 1; // n+1
            let mut bn = 0.0; // B_n
            for (k, bk) in result.iter().enumerate() {
                bn += np1.choose(k) as f64 * bk;
            }
            bn = -bn / np1 as f64;
            result.push(bn);
        }
    }
    result
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_bernoulli_rec() {
        for i in 0..4 {
            assert_eq!(bernoulli_rec(i).len(), i + 1);
        }

        //println!("{:?}", bernoulli_rec(4));
        //assert!(false);
        eprintln!("{:?}", bernoulli_rec(4));
        eprintln!("{}", -1.0 / 30.0);
        assert!(false);
        //assert!(false);
        //assert_eq!(bernoulli_rec(4)[4], -1.0 / 30.0);
        // assert_eq!(*bernoulli_rec(6).last().unwrap(), 1.0 / 42.0);
        // assert_eq!(*bernoulli_rec(8).last().unwrap(), -1.0 / 30.0);
        // assert_eq!(*bernoulli_rec(10).last().unwrap(), 5.0 / 66.0);
        //
        // // Odds should be zero
        // assert_eq!(*bernoulli_rec(5).last().unwrap(), 0.0);
    }
}
