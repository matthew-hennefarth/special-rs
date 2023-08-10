//**********************************************************************
// This file is part of sci-rs                                         *
// Copyright 2023 Matthew R. Hennefarth                                *
//**********************************************************************

use crate::special::tools::eval_poly;
use crate::special::{frexp, ldexp};
use num_traits::Float;

pub(crate) trait RealErfConsts: Sized {
    const C: [Self; 2];
    const Y: [Self; 5];
    const P_1: [Self; 6];
    const Q_1: [Self; 6];
    const P_2: [Self; 8];
    const Q_2: [Self; 7];
    const P_3: [Self; 7];
    const Q_3: [Self; 7];
    const P_4: [Self; 7];
    const Q_4: [Self; 7];
    const P_5: [Self; 9];
    const Q_5: [Self; 9];
}

macro_rules! impl_realerfconsts {
    ($($T: ty)*) => ($(
        impl RealErfConsts for $T {
            const C: [Self; 2] = [1.125, 0.003379167095512573896158903121545171688];
            const Y: [Self; 5] = [1.044948577880859375, 0.405935764312744140625, 0.50672817230224609375, 0.5405750274658203125, 0.55825519561767578125];
            const P_1: [Self; 6] = [0.0834305892146531988966,-0.338097283075565413695, -0.0509602734406067204596, -0.00904906346158537794396, -0.000489468651464798669181, -0.200305626366151877759e-4];
            const Q_1: [Self; 6] = [1.0, 0.455817300515875172439, 0.0916537354356241792007, 0.0102722652675910031202, 0.000650511752687851548735, 0.189532519105655496778e-4];
            const P_2: [Self; 8] = [-0.0980905922162812031672,
            0.159989089922969141329,
            0.222359821619935712378,
            0.127303921703577362312,
            0.0384057530342762400273,
            0.00628431160851156719325,
            0.000441266654514391746428,
            0.266689068336295642561e-7];
            const Q_2: [Self; 7] = [
            1.0,
            2.03237474985469469291,
            1.78355454954969405222,
            0.867940326293760578231,
            0.248025606990021698392,
            0.0396649631833002269861,
            0.00279220237309449026796,];
            const P_3: [Self; 7] = [
            -0.024350047620769840217,
            0.0343522687935671451309,
            0.0505420824305544949541,
            0.0257479325917757388209,
            0.00669349844190354356118,
            0.00090807914416099524444,
            0.515917266698050027934e-4
            ];
            const Q_3: [Self; 7] = [
            1.0,
            1.71657861671930336344,
            1.26409634824280366218,
            0.512371437838969015941,
            0.120902623051120950935,
            0.0158027197831887485261,
            0.000897871370778031611439
            ];
            const P_4: [Self; 7] = [
            0.0029527671653097284033,
            0.0141853245895495604051,
            0.0104959584626432293901,
            0.00343963795976100077626,
            0.00059065441194877637899,
            0.523435380636174008685e-4,
            0.189896043050331257262e-5
            ];
            const Q_4: [Self; 7] = [
            1.0,
            1.19352160185285642574,
            0.603256964363454392857,
            0.165411142458540585835,
            0.0259729870946203166468,
            0.00221657568292893699158,
            0.804149464190309799804e-4
            ];
            const P_5: [Self; 9] = [
            0.00593438793008050214106,
            0.0280666231009089713937,
            -0.141597835204583050043,
            -0.978088201154300548842,
            -5.47351527796012049443,
            -13.8677304660245326627,
            -27.1274948720539821722,
            -29.2545152747009461519,
            -16.8865774499799676937
            ];
            const Q_5: [Self; 9] = [
            1.0,
            4.72948911186645394541,
            23.6750543147695749212,
            60.0021517335693186785,
            131.766251645149522868,
            178.167924971283482513,
            182.499390505915222699,
            104.365251479578577989,
            30.8365511891224291717,
            ];
        }
)*)
}

impl_realerfconsts! {f32 f64}

pub(crate) fn r_erf<T>(x: T, mut compliment: bool) -> T
where
    T: Float + RealErfConsts,
{
    if x.is_sign_negative() {
        return if !compliment {
            -r_erf(-x, compliment)
        } else if x < T::from(-0.5).unwrap() {
            T::from(2.0).unwrap() - r_erf(-x, compliment)
        } else {
            T::one() + r_erf(-x, false)
        };
    }

    let result = if x < T::from(0.5).unwrap() {
        // Calculating erf
        if x.is_zero() {
            T::zero()
        } else if x < T::from(1e-10).unwrap() {
            x * T::C[0] + T::C[1] * x
        } else {
            x * (T::Y[0] + eval_poly(x * x, &T::P_1) / eval_poly(x * x, &T::Q_1))
        }
    } else if if compliment {
        x < T::from(110.0).unwrap()
    } else {
        x < T::from(6.6).unwrap()
    } {
        // We calculate erfc here
        compliment = !compliment;
        let r = if x < T::from(1.5).unwrap() {
            T::Y[1]
                + eval_poly(x - T::from(0.5).unwrap(), &T::P_2)
                    / eval_poly(x - T::from(0.5).unwrap(), &T::Q_2)
        } else if x < T::from(2.5).unwrap() {
            T::Y[2]
                + eval_poly(x - T::from(1.5).unwrap(), &T::P_3)
                    / eval_poly(x - T::from(1.5).unwrap(), &T::Q_3)
        } else if x < T::from(4.5).unwrap() {
            T::Y[3]
                + eval_poly(x - T::from(3.5).unwrap(), &T::P_4)
                    / eval_poly(x - T::from(3.5).unwrap(), &T::Q_4)
        } else {
            T::Y[4] + eval_poly(x.recip(), &T::P_5) / eval_poly(x.recip(), &T::Q_5)
        };
        let (hi, expon): (T, isize) = frexp(x);
        let hi = ldexp(hi, 32).floor();
        let hi = ldexp(hi, expon - 32);
        let lo = x - hi;
        let x_sqr = x * x;
        let err_sqr = ((hi * hi - x_sqr) + T::from(2.0).unwrap() * hi * lo) + lo * lo;
        r * (-x_sqr).exp() * (-err_sqr).exp() / x
    } else {
        compliment = !compliment;
        T::zero()
    };

    if compliment {
        T::one() - result
    } else {
        result
    }
}
