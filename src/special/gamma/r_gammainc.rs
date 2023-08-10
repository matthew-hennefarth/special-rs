//**********************************************************************
// This file is part of sci-rs                                         *
// Copyright 2023 Matthew R. Hennefarth                                *
//**********************************************************************

use num_traits::Float;

enum Sign {
    Positive,
    Negative,
}

/// log(1+x)-x
fn log1pmx<T>(x: T) -> T
where
    T: Float,
{
    // TODO: maybe replace with this C code?
    // if (fabs(x) < 0.5) {
    //     int n;
    //     double xfac = x;
    //     double term;
    //     double res = 0;
    //
    //     for(n = 2; n < MAXITER; n++) {
    //         xfac *= -x;
    //         term = xfac / n;
    //         res += term;
    //         if (fabs(term) < MACHEP * fabs(res)) {
    //         break;
    //         }
    //     }
    //     return res;
    // }
    // else {
    //     return log1p(x) - x;
    // }
    (T::one() + x).ln() - x
}

fn asymptotic_series<T>(s: T, x: T, sign: Sign) -> T
where
    T: Float,
{
    let sign = match sign {
        Sign::Positive => T::one(),
        Sign::Negative => -T::one(),
    };

    let lambda = x / s;
    let sigma = (x - s) / s;
    let two: T = T::one() + T::one();

    let eta = if lambda > T::one() {
        (-two * log1pmx(sigma)).sqrt()
    } else if lambda < T::one() {
        -(-two * log1pmx(sigma)).sqrt()
    } else {
        T::zero()
    };
    let mut res = erfc(sign * eta * (s / two).sqrt()) / two;

    s
}

pub(crate) trait RealLowerGammaIncConsts: Sized {
    const SMALL: Self;
    const LARGE: Self;
    const SMALL_RATIO: Self;
    const LARGE_RATIO: Self;
}

macro_rules! impl_reallowergammaincconsts {
    ($($T: ty)*) => ($(
        impl RealLowerGammaIncConsts for $T {
            const SMALL: Self = 20.0;
            const LARGE: Self = 200.0;
            const SMALL_RATIO: Self = 0.3;
            const LARGE_RATIO: Self = 4.5;
        }
)*)
}

impl_reallowergammaincconsts! {f32 f64}

/// Implementation of the regularized incomplete lower gamma function.
///
/// Implementation follows that from the [cephes library in scipy](cephes).
///
/// [cephes]: https://github.com/scipy/scipy/blob/main/scipy/special/cephes/igam.c#L368
pub(crate) fn r_gammainc<T>(s: T, x: T) -> T
where
    T: Float + RealLowerGammaIncConsts,
{
    if s < T::zero() || x < T::zero() {
        panic!("Invalid!");
    } else if s.is_zero() {
        if x.is_sign_positive() {
            return T::one();
        }
        return T::nan();
    } else if x.is_zero() {
        // Zero integration limit
        return T::zero();
    } else if s.is_infinite() {
        if x.is_infinite() {
            return T::nan();
        }
        return T::zero();
    } else if x.is_infinite() {
        return T::one();
    }

    let absxma_a = (x - s).abs() / s;
    if s > T::SMALL && s < T::LARGE && absxma_a < T::SMALL_RATIO {
        return asymptotic_series(s, x, Sign::Positive);
    } else if s > T::LARGE && absxma_a < (T::LARGE_RATIO / s.sqrt()) {
        return asymptotic_series(s, x, Sign::Positive);
    }
    s
}
