//! EKF版の逐次較正アルゴリズム
//! 
//! 非直行性補正係数（D12, D13, D23）を推定値から外した
//! 
//! 参考文献
//! 1. John L. Crassdis, Kok-Lam Lai, Richard R. Harman, 
//!    "Real-Time Attitude-Independent Three-Axis Magnetometer Calibration", 
//!    Journal of Guidance, Control, and Dynamics, Vol.28, No.1, 2005.

use super::{SVector, SMatrix, NOISE_VAR, FIELD_NORM};

type Scalar = SVector<f64, 1>;  // 1x1行列とf64を足せないからスカラー型を作る
type Vector3 = SVector<f64, 3>;
type Vector6 = SVector<f64, 6>;
type Matrix1x3 = SMatrix<f64, 1, 3>;
type Matrix1x6 = SMatrix<f64, 1, 6>;
type Matrix3x3 = SMatrix<f64, 3, 3>;
type Matrix6x6 = SMatrix<f64, 6, 6>;

pub struct EKF {
    pub xhat: Vector6,  // 状態変数（先頭3要素がバイアス，その後ろ3要素がスケールファクター）
    p: Matrix6x6, // 誤差共分散行列
}

impl EKF {
    pub fn new() -> Self {
        let mut p_init = Matrix6x6::zeros();
        for i in 0..3 {
            p_init[(i, i)]     = 0.1;   // バイアス
            p_init[(i+3, i+3)] = 0.001; // スケールファクター
        }
        
        Self {
            xhat: Vector6::zeros(),
            p: p_init,
        }
    }
    
    pub fn update(&mut self, mag: Vector3) {
        let c = self.jacobian_h(mag);

        // ノイズの共分散行列(Ref1 - eq.(5c))は定数になるので計算後の値を埋め込んである
        let tmp = self.calibrate(mag);
        let sigma_square = 4.0 * tmp.transpose() * NOISE_VAR * tmp + Scalar::new(2.0 * (3.0 * NOISE_VAR * NOISE_VAR));
        
        let coef: f64 = 1.0 / (c * self.p * c.transpose() + sigma_square)[0];  // スカラーになるので逆行列計算は不要
        let gain = self.p * c.transpose() * coef;
        let z =  Scalar::new(mag.norm_squared() - FIELD_NORM * FIELD_NORM);
        self.xhat = self.xhat + gain * (z - self.h(mag));
        self.p = (Matrix6x6::identity() - gain * c) * self.p;
    }

    fn get_bias_scale_factor(&self) -> (Vector3, Matrix3x3) {
        (
            self.xhat.fixed_rows::<3>(0).into_owned(),
            Matrix3x3::from_diagonal(&self.xhat.fixed_rows::<3>(3))
        )
    }

    /// 推定したバイアスとスケールファクタを用いて較正した地磁気計測値を返す
    /// 
    /// * mag: 生の地磁気計測値 [x, y, z]
    pub fn calibrate(&self, mag: Vector3) -> Vector3 {
        let mut calib = Vector3::zeros();
        for i in 0..3 {
            calib[i] = (mag[i] - self.xhat[i]) + self.xhat[i+3] * mag[i];
        }

        calib
    }

    fn h(&self, mag: Vector3) -> Scalar {
        let (b, d) = self.get_bias_scale_factor();
        -mag.transpose() * (2.0 * d + d * d) * mag + 2.0 * mag.transpose() * (Matrix3x3::identity() + d) * b - Scalar::new(b.norm_squared())
    }

    fn jacobian_h(&self, mag: Vector3) -> Matrix1x6 {
        let (b, d) = self.get_bias_scale_factor();

        let mut s = Matrix1x3::zeros();
        let mut j = Matrix1x3::zeros();
        let mut de_dd = Matrix3x3::zeros();
        // 初期化
        for i in 0..3 {
            s[i] = mag[i] * mag[i];
            j[i] = mag[i] * b[i];
            de_dd[(i, i)] = 2.0 * (1.0 + d[(i,i)]);
        }
        
        let left = 2.0 * mag.transpose() * (Matrix3x3::identity() + d) -  2.0 * b.transpose();
        let right = -s * de_dd + 2.0 * j;
        let mut dh_dx = Matrix1x6::zeros();
        for i in 0..3 {  // 行列を結合する方法無いんだっけ
            dh_dx[(0, i)] = left[(0, i)];
            dh_dx[(0, i + 3)] = right[(0, i)];
        }

        dh_dx
    }
}