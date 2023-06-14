//**********************************************************************
// This file is part of Sci-rs                                         *
// Copyright 2023 Matthew R. Hennefarth                                *
//**********************************************************************

use crate::special::gamma::gamma_util::{eval_poly, lngamma_stirling, LnGammaStirlingConsts};
use crate::traits::FloatSciConst;
use num_complex::{Complex, ComplexFloat};
use num_traits::{Float, FloatConst, One, Zero};
use std::ops::{Add, AddAssign, Mul, MulAssign, Sub};

/// Coefficients should be zeta(i)/i)
pub(crate) trait LogGammaTaylorCoeffs: Sized {
    const COEFFS: [Self; 23];
}

macro_rules! impl_loggamma_taylorconsts {
    ($($T: ty)*) => ($(
        impl LogGammaTaylorCoeffs for $T {
            const COEFFS: [Self; 23] = [
                -4.3478266053040259361e-2, 4.5454556293204669442e-2,
                -4.7619070330142227991e-2, 5.000004769810169364e-2,
                -5.2631679379616660734e-2, 5.5555767627403611102e-2,
                -5.8823978658684582339e-2, 6.2500955141213040742e-2,
                -6.6668705882420468033e-2, 7.1432946295361336059e-2,
                -7.6932516411352191473e-2, 8.3353840546109004025e-2,
                -9.0954017145829042233e-2, 1.0009945751278180853e-1,
                -1.1133426586956469049e-1, 1.2550966952474304242e-1,
                -1.4404989676884611812e-1, 1.6955717699740818995e-1,
                -2.0738555102867398527e-1, 2.7058080842778454788e-1,
                -4.0068563438653142847e-1, 8.2246703342411321824e-1,
                -5.7721566490153286061e-1
            ];
        }
)*)
}

impl_loggamma_taylorconsts! {f32 f64}

/// Taylor Series of Log-Gamma around $z=1$.
/// $$
/// \ln\Gamma(z+1) = -\gamma z + \frac{\zeta(2)}{2}z^2 - \frac{\zeta(3)}{3}z^3\ldots
/// $$
/// Here $\gamma$ is the Euler-Mascheroni constant.
fn loggamma_taylor<T>(z: T) -> T
where
    T: ComplexFloat
        + Add<<T as ComplexFloat>::Real, Output = T>
        + Mul<<T as ComplexFloat>::Real, Output = T>,
    <T as ComplexFloat>::Real: LogGammaTaylorCoeffs,
{
    let z = z - T::one();
    z * eval_poly(z, &T::Real::COEFFS)
}

pub(crate) trait LogGammaConsts {
    const MIN_USE_STIRLING: Self;
    const TAYLOR_RADIUS: Self;
    const MAX_REFLECT: Self;
}

macro_rules! impl_loggammaconsts {
    ($($T: ty)*) => ($(
        impl LogGammaConsts for $T {
            const MIN_USE_STIRLING: Self = 7.0;
            const TAYLOR_RADIUS: Self = 0.2;
            const MAX_REFLECT: Self = 0.1;
        }
)*)
}

impl_loggammaconsts! {f32 f64}

/// Recurse upwards until we can use the stirling approximation
fn c_lngamma_recurrence<T>(mut z: T) -> Complex<<T as ComplexFloat>::Real>
where
    T: ComplexFloat
        + AddAssign
        + MulAssign
        + Add<<T as ComplexFloat>::Real, Output = T>
        + Mul<<T as ComplexFloat>::Real, Output = T>,
    <T as ComplexFloat>::Real: LogGammaConsts + AddAssign + FloatSciConst + LnGammaStirlingConsts,
{
    let mut product = z;
    let mut sb = false;
    let mut signflips = T::Real::zero();

    z += T::one();
    while z.re() <= T::Real::MIN_USE_STIRLING {
        product *= z;
        if product.im().is_sign_negative() && !sb {
            signflips += T::Real::one();
        }
        sb = product.im().is_sign_negative();
        z += T::one();
    }

    let partial_result = lngamma_stirling(z) - product.ln();

    return Complex::new(partial_result.re(), partial_result.im())
        - Complex::i() * signflips * T::Real::TAU();
}

/// Complex implementation of `Log(Gamma(z))`
///
/// Input can be a real-valued float, but output will always be complex.
///
/// For $z.re$ > 7 or $|z.im| > 2$ we use the stirling approximation.
///
/// For $z$ close to 1 or 2, we use a Taylor expansion around 1.
///
/// Based off of the SciPy implementation: https://github.com/scipy/scipy/blob/main/scipy/special/_loggamma.pxd
pub(crate) fn c_lngamma<T>(z: T) -> Complex<<T as ComplexFloat>::Real>
where
    T: ComplexFloat
        + AddAssign
        + MulAssign
        + Add<<T as ComplexFloat>::Real, Output = T>
        + Mul<<T as ComplexFloat>::Real, Output = T>,
    <T as ComplexFloat>::Real: LogGammaConsts
        + FloatSciConst
        + LnGammaStirlingConsts
        + LogGammaTaylorCoeffs
        + Add<Output = <T as ComplexFloat>::Real>
        + AddAssign,
    Complex<<T as ComplexFloat>::Real>: Sub<T>,
{
    if !z.is_finite() {
        return Complex::new(z.re(), z.im());
    }
    if z.re() <= T::Real::zero() && z.re() == z.re().floor() && z.im().is_zero() {
        return Complex::new(T::Real::nan(), T::Real::nan());
    }

    if z.re() > T::Real::MIN_USE_STIRLING || Float::abs(z.im()) > T::Real::MIN_USE_STIRLING {
        let result = lngamma_stirling(z);
        return Complex::new(result.re(), result.im());
    }
    if (z - T::one()).abs() <= T::Real::TAYLOR_RADIUS {
        let result = loggamma_taylor(z);
        return Complex::new(result.re(), result.im());
    }
    if (z - T::one() - T::one()).abs() <= T::Real::TAYLOR_RADIUS {
        let z = z - T::one();
        let result = z.ln() + loggamma_taylor(z);
        return Complex::new(result.re(), result.im());
    }
    if z.re() < T::Real::MAX_REFLECT {
        let two = T::Real::one() + T::Real::one();
        let four = two * two;

        let tmp =
            z.im().signum() * T::Real::TAU() * ((two * z.re() + T::Real::one()) / four).floor();

        // Winding issue I believe like ln(z) = ln(r) + i\theta
        let winding = Complex::new(Float::ln(T::Real::PI()), tmp);

        // Reflection PHase
        let phase = -(z * T::Real::PI()).sin().ln();
        let phase = Complex::new(phase.re(), phase.im());

        // Reflected
        let partial_result = -c_lngamma(T::one() - z);

        return phase + partial_result + winding;
    }

    if !z.im().is_sign_negative() {
        return c_lngamma_recurrence(z);
    }
    c_lngamma_recurrence(z.conj()).conj()
}

/// Gamma function evaluated from `lngamma`.
/// $$
/// \Gamma(z) = e^{\ln\Gamma(z)}
/// $$
pub(crate) fn c_gamma<T>(z: T) -> Complex<<T as ComplexFloat>::Real>
where
    T: ComplexFloat
        + AddAssign
        + MulAssign
        + Add<<T as ComplexFloat>::Real, Output = T>
        + Mul<<T as ComplexFloat>::Real, Output = T>,
    <T as ComplexFloat>::Real: LogGammaConsts
        + FloatSciConst
        + LnGammaStirlingConsts
        + LogGammaTaylorCoeffs
        + Add<Output = <T as ComplexFloat>::Real>
        + AddAssign,
    Complex<<T as ComplexFloat>::Real>: Sub<T>,
{
    c_lngamma(z).exp()
}

/// Reciprocal of the Gamma function evaluated from `lngamma`.
/// $$
/// \frac{1}{\Gamma(z)} = e^{-\ln\Gamma(z)}
/// $$
pub(crate) fn c_rgamma<T>(z: T) -> Complex<<T as ComplexFloat>::Real>
where
    T: ComplexFloat
        + AddAssign
        + MulAssign
        + Add<<T as ComplexFloat>::Real, Output = T>
        + Mul<<T as ComplexFloat>::Real, Output = T>,
    <T as ComplexFloat>::Real: LogGammaConsts
        + FloatSciConst
        + LnGammaStirlingConsts
        + LogGammaTaylorCoeffs
        + Add<Output = <T as ComplexFloat>::Real>
        + AddAssign,
    Complex<<T as ComplexFloat>::Real>: Sub<T>,
{
    (-c_lngamma(z)).exp()
}

pub(crate) fn c_lgamma<T>(z: T) -> <T as ComplexFloat>::Real
where
    T: ComplexFloat
        + AddAssign
        + MulAssign
        + Add<<T as ComplexFloat>::Real, Output = T>
        + Mul<<T as ComplexFloat>::Real, Output = T>,
    <T as ComplexFloat>::Real: LogGammaConsts
        + FloatSciConst
        + LnGammaStirlingConsts
        + LogGammaTaylorCoeffs
        + Add<Output = <T as ComplexFloat>::Real>
        + AddAssign,
    Complex<<T as ComplexFloat>::Real>: Sub<T>,
{
    c_lngamma(z).re()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::special::gamma::r_lgamma;

    const PRECISION: f64 = 1.0e-14;

    #[test]
    fn test_c_loggamma() {
        // Values take from Wolframalpha
        assert_almost_eq!(
            c_lngamma(Complex { re: 7.5, im: 0.0 }),
            Complex {
                re: 7.53436423675873295515836763243,
                im: 0.0
            },
            PRECISION
        );
        assert_almost_eq!(
            c_lngamma(Complex { re: 7.5, im: 1.0 }),
            Complex {
                re: 7.46329489273832466759034022127,
                im: 1.95012140717825057478945773428
            },
            PRECISION
        );
        assert_almost_eq!(
            c_lngamma(0.9),
            Complex {
                re: r_lgamma(0.9),
                im: 0.0
            },
            PRECISION
        );
        assert_almost_eq!(
            c_lngamma(Complex { re: 0.95, im: 0.05 }),
            Complex {
                re: 0.028753589066549549997223245584,
                im: -0.03307300849770222888977061702
            },
            PRECISION
        );

        // From SciPy version 1.10.1
        assert_almost_eq!(
            c_lngamma(Complex { re: 0.0, im: 7.1 }),
            Complex {
                re: -11.2137627790627281143543,
                im: 6.0195299083889901581301
            },
            PRECISION
        );
        // From SciPy version 1.10.1
        assert_almost_eq!(
            c_lngamma(Complex { re: 0.0, im: -7.22 }),
            Complex {
                re: -11.4106384227068478054434,
                im: -6.2559451620947514882687
            },
            PRECISION
        );
        // From SciPy version 1.10.1
        assert_almost_eq!(
            c_lngamma(-7.2),
            Complex {
                re: -7.2548056050797278260234,
                im: -25.1327412287183449279837
            },
            PRECISION
        );
        // From SciPy version 1.10.1
        assert_almost_eq!(
            c_lngamma(Complex { re: -7.33, im: 4.3 }),
            Complex {
                re: -19.7421179905029617884793,
                im: -15.5482125579457619579671
            },
            PRECISION
        );
        // From SciPy version 1.10.1
        assert_almost_eq!(
            c_lngamma(Complex { re: 1.5, im: 0.5 }),
            Complex {
                re: -0.2341863474703487213446,
                im: 0.0346689612753978693149
            },
            PRECISION
        );
        // From SciPy version 1.10.1
        assert_almost_eq!(
            c_lngamma(Complex { re: 1.5, im: -0.5 }),
            Complex {
                re: -0.2341863474703487213446,
                im: -0.0346689612753978693149
            },
            PRECISION
        );
        // From SciPy version 1.10.1
        assert_almost_eq!(
            c_lngamma(Complex {
                re: 140.233,
                im: -54.21
            }),
            Complex {
                re: 541.1630267673411935902550,
                im: -269.0855433781250667379936
            },
            PRECISION
        );
        // From SciPy version 1.10.1
        assert_almost_eq!(
            c_lngamma(Complex { re: 0.0, im: 1.0 }),
            Complex {
                re: -0.6509231993018549378149,
                im: -1.8724366472624294210902
            },
            PRECISION
        );

        assert_almost_eq!(
            c_lngamma(1.5),
            Complex {
                re: r_lgamma(1.5),
                im: 0.0
            },
            PRECISION
        );
    }

    #[test]
    fn test_c_lgamma() {
        for i in 1..5 {
            assert_eq!(
                c_lgamma(Complex {
                    re: 0.0,
                    im: i as f32
                }),
                c_lgamma(Complex {
                    re: 0.0,
                    im: -i as f32
                })
            );
        }
        // From Wolframalpha
        assert_almost_eq!(
            c_lgamma(Complex { re: 0.0, im: 1.0 }),
            -0.6509231993018563388852168315,
            PRECISION
        );
        assert_almost_eq!(
            c_lgamma(Complex { re: 0.0, im: 2.0 }),
            -2.56922596699087465064722769985,
            PRECISION
        );
        assert_almost_eq!(
            c_lgamma(Complex { re: 0.0, im: 3.0 }),
            -4.342756588257865882968429589295,
            PRECISION
        );
        assert_almost_eq!(
            c_lgamma(Complex { re: 0.0, im: 4.0 }),
            -6.057393954528778266207447521547,
            PRECISION
        );
        assert_almost_eq!(
            c_lgamma(Complex { re: 0.0, im: 5.0 }),
            -7.739762056986849186171316767808,
            PRECISION
        );
        assert_almost_eq!(
            c_lgamma(Complex {
                re: 0.001,
                im: 0.001
            }),
            6.56060447383755263956515723187072,
            PRECISION
        );
        assert_almost_eq!(
            c_lgamma(Complex {
                re: 1.0e-9,
                im: 1.0e-9
            }),
            20.3766922460892228365517741716250487,
            PRECISION
        );
        assert_almost_eq!(
            c_lgamma(Complex {
                re: 150.01,
                im: 1.0e-20
            }),
            600.059543872337641753611684656568302,
            PRECISION
        );
    }
}
