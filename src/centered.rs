//! 最小二乗法とCentering近似に基づく地磁気センサの逐次較正アルゴリズム
//! 
//! Ref.1を参考にしつつ，非直行性補正係数は要らないからバイアスと
//! スケールファクタの推定だけに変更したもの．
//! 
//! 参考文献
//! 1. John L. Crassdis, Kok-Lam Lai, Richard R. Harman, 
//!    "Real-Time Attitude-Independent Three-Axis Magnetometer Calibration", 
//!    Journal of Guidance, Control, and Dynamics, Vol.28, No.1, 2005.
//! 2. Robert Alonso, Malcolm D.Shuster, 
//!    "Complete Linear Attitude-Independent Magnetometer Calibration",
//!     The Journal of the Astronautical Sciences, Vol.50, No.4, 2002.
//! 3. John L. Crassidis, F. Landis Markley, E. Glenn Lightsey,
//!    "Global Positioning System Integer Ambiguity Resolution Without Attitude Knowledge",
//!    Journal of Guidance, Control, and Dynamics, Vol.22, No.2,1999.
//! 4. Roberto Alonso, Malcolm D. Shuster,
//!    "TWOSTEP: A Fast Robust Algorithm for Attitude-Independent Magnetometer-Bias Determination",
//!     The Journal of the Astronautical Sciences, Vol.50, No.4, 2002.

use super::{SVector, SMatrix, NOISE_VAR, FIELD_NORM};

type Vector3 = SVector<f64, 3>;
type Vector6 = SVector<f64, 6>;
type Matrix1x6 = SMatrix<f64, 1, 6>;
type Matrix3x3 = SMatrix<f64, 3, 3>;
type Matrix6x6 = SMatrix<f64, 6, 6>;


pub struct SequencialCentered {
    theta_dash: Vector6, // Ref.1 - eq.(4b)
    p: Matrix6x6,        // Ref.1 - eq.(9b)
    l_bar: Matrix1x6,    // Ref.1 - eq.(4a)
    z_bar: f64,          // Ref.1 - eq.(3)
    mu_bar: f64,         // Ref.1 - eq.(5a)
    sigma_bar: f64,      // Ref.1 - Square of eq.(7)
}

impl SequencialCentered {
    pub fn new() -> Self {
        // 初期値は適当
        let b_init = 0.01;  // ゼロだと収束が保証されないらしい（Ref.2）
        let d_init = 0.001; // スケールファクタの初期値．ゼロだとNaNになる．
        let mut theta_dash_init = Vector6::zeros();
        for i in 0..3 {
            theta_dash_init[i] = (1.0 + d_init) * b_init;
            theta_dash_init[i + 3] = d_init * (2.0 + d_init);
        }

        // ここは適当に決めてもなんとかなる
        let mut p_init = Matrix6x6::zeros();
        for i in 0..3 {
            p_init[(i, i)]     = 0.1;   // バイアス
            p_init[(i+3, i+3)] = 0.001; // スケールファクタ
        }

        // sigmaの初期値をゼロにするとゼロ除算が起こるのでそれっぽい初期値を入れておく
        let tmp = (Matrix3x3::identity() + d_init * Matrix3x3::identity()) * Vector3::new(0.0, FIELD_NORM, 0.0);
        let sigma_bar_init = (4.0 * tmp.transpose() * NOISE_VAR * tmp)[0] + 2.0 * (3.0 * NOISE_VAR * NOISE_VAR);

        let mut l_bar_init = Matrix1x6::zeros();
        l_bar_init[(0, 2)] = 2.0 * FIELD_NORM;
        l_bar_init[(0, 4)] = -FIELD_NORM * FIELD_NORM;

        // 逐次計算で変数の値は素早く収束していくので初期値は適当で良い（Ref.3）．
        // ただし，sigma_barに関しては初期値がゼロだとゼロ除算が起こるので注意．
        Self {
            theta_dash: theta_dash_init,
            p: p_init,
            l_bar: l_bar_init,
            z_bar: 0.0,
            mu_bar: -3.0 * NOISE_VAR,
            sigma_bar: sigma_bar_init,
        }
    }

    /// 推定値を更新する
    /// * mag: 地磁気計測値 [x, y, z]
    pub fn update(&mut self, mag: Vector3) {
        // ノイズの共分散行列(Ref1 - eq.(5c))は定数になるので計算後の値を埋め込んである
        let tmp = self.calibrate(mag);
        let sigma = (4.0 * tmp.transpose() * NOISE_VAR * tmp)[0] + 2.0 * (3.0 * NOISE_VAR * NOISE_VAR);

        let coef = 1.0 / (sigma + self.sigma_bar);
        self.sigma_bar = ((1.0 / self.sigma_bar) + (1.0 / sigma)).recip();

        let mut l = Matrix1x6::zeros();
        for i in 0..3 {
            l[(0, i)] = 2.0 * mag[i];
            l[(0, i+3)] = -mag[i] * mag[i];
        }
        self.l_bar = coef * (sigma * self.l_bar + self.sigma_bar * l);

        let z = mag.norm_squared() - (FIELD_NORM * FIELD_NORM);
        self.z_bar = coef * (sigma * self.z_bar + self.sigma_bar * z);

        let mu = -3.0 * NOISE_VAR;
        self.mu_bar = coef * (sigma * self.mu_bar + self.sigma_bar * mu);
        
        let l_tilde = l - self.l_bar;
        let z_tilde = z - self.z_bar;
        let mu_tilde = mu - self.mu_bar;

        let tmp = ((l_tilde * self.p * l_tilde.transpose())[0] + sigma).recip();
        let gain = Matrix6x6::identity() - self.p * l_tilde.transpose() * tmp * l_tilde;
        self.p = gain * self.p;
        self.theta_dash = gain * self.theta_dash + sigma.recip() * (z_tilde - mu_tilde) * self.p * l_tilde.transpose();
    }

    /// Get magnetometer bias and scale-factor
    /// 
    /// Return: (Bias, Scale-factor) = ([b1, b2, b3], [D11, D22, D33])
    pub fn get_bias_scale_factor(&self) -> (Vector3, Vector3) {
        let mut b = Vector3::zeros();  // バイアス
        let mut d = Vector3::zeros();  // スケールファクタ
        for i in 0..3 {
            let tmp = (1.0 + self.theta_dash[i + 3]).sqrt();
            b[i] = self.theta_dash[i] / tmp;
            d[i] = tmp - 1.0;
        }
        (b, d)
    }

    /// 推定したバイアスとスケールファクタを用いて較正した地磁気計測値を返す
    /// 
    /// * mag: 生の地磁気計測値 [x, y, z]
    pub fn calibrate(&self, mag: Vector3) -> Vector3 {
        let (b, d) = self.get_bias_scale_factor();

        let mut calib = Vector3::zeros();
        for i in 0..3 {
            calib[i] = (mag[i] - b[i]) + d[i] * mag[i];
        }

        calib
    }
}