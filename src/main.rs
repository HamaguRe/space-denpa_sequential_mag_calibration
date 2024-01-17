//! 地磁気センサの逐次較正アルゴリズムを実装する

use std::fs::{self, File};
use std::io::{Write, BufWriter};
use nalgebra::{SVector, SMatrix};

type Vector3 = SVector<f64, 3>;
type Vector4 = SVector<f64, 4>;
type Matrix4x4 = SMatrix<f64, 4, 4>;

/// 観測ノイズの共分散（各軸の干渉は無く，3軸全て同じ条件と仮定する）
const NOISE_VAR: f64 = 0.04;  // > 0

/// 地磁気ベクトルのノルム（全磁力）
const FIELD_NORM: f64 = 0.47;  // [Gauss]

mod centered;
mod ekf;
mod ukf;


fn main() {
    // 逐次較正に使用するデータを読みこむ
    let mag_samples = load_csv("./magnetic_data/data-set_1.csv");
    let mut writer_centered = BufWriter::new(File::create("./result_centered.csv").unwrap());
    let mut writer_ekf = BufWriter::new(File::create("./result_ekf.csv").unwrap());
    let mut writer_ukf = BufWriter::new(File::create("./result_ukf.csv").unwrap());

    let mut centered = centered::SequencialCentered::new();
    let mut ekf = ekf::EKF::new();
    let mut ukf = ukf::UKF::new();

    // 収束値を見る
    /*
    for _ in 0..500 {
        for mag in &mag_samples {
            // 更新
            centered.update(*mag);
            ekf.update(*mag);
            ukf.update(*mag);
        }
    }
    */

    let csv_header = b"TimeStep,b_x,b_y,b_z,D11,D22,D33,calib_x,calib_y,calib_z\n";
    writer_centered.write(csv_header).unwrap();
    writer_ekf.write(csv_header).unwrap();
    writer_ukf.write(csv_header).unwrap();
    for k in 0..mag_samples.len() {
        // 推定値更新
        centered.update(mag_samples[k]);
        ekf.update(mag_samples[k]);
        ukf.update(mag_samples[k]);

        // ------ 計算結果書き出し ------ //
        // ---- Centeredアルゴリズム
        // 計算ステップ
        writer_centered.write(format!("{}", k).as_bytes()).unwrap();
        // スケールファクタとバイアス
        let (b, d) = centered.get_bias_scale_factor();
        for i in 0..3 {
            writer_centered.write(format!(",{}", b[i]).as_bytes()).unwrap();
        }
        for i in 0..3 {  // D11, D22, D33
            writer_centered.write(format!(",{}", d[i]).as_bytes()).unwrap();
        }
        // 逐次較正した地磁気計測値
        let calibrated = centered.calibrate(mag_samples[k]);
        for i in 0..3 {
            writer_centered.write(format!(",{}", calibrated[i]).as_bytes()).unwrap();
        }
        writer_centered.write(b"\n").unwrap();

        // ---- EKF
        // 計算ステップ
        writer_ekf.write(format!("{}", k).as_bytes()).unwrap();
        // バイアス，スケールファクタ 各3要素
        for val in &ekf.xhat {
            writer_ekf.write(format!(",{}", val).as_bytes()).unwrap();
        }
        // 逐次較正した地磁気計測値
        let ekf_calibrated = ekf.calibrate(mag_samples[k]);
        for val in &ekf_calibrated {
            writer_ekf.write(format!(",{}", val).as_bytes()).unwrap();
        }
        writer_ekf.write(b"\n").unwrap();

        // ---- UKF
        // 計算ステップ
        writer_ukf.write(format!("{}", k).as_bytes()).unwrap();
        // バイアス，スケールファクタ 各3要素
        for val in &ukf.xhat {
            writer_ukf.write(format!(",{}", val).as_bytes()).unwrap();
        }
        // 逐次較正した地磁気計測値
        let ukf_calibrated = ukf.calibrate(mag_samples[k]);
        for val in &ukf_calibrated {
            writer_ukf.write(format!(",{}", val).as_bytes()).unwrap();
        }
        writer_ukf.write(b"\n").unwrap();
        // --- 計算結果書き出しここまで --- //
    }

    // --- 以下，計算結果の出力 ---
    // 最終ステップにおけるバイアスとスケールファクタ
    println!("--- Bias and scale factor in the final step ---");
    let (b, d) = centered.get_bias_scale_factor();
    println!("Centered: b = {:?}, d = {:?}", b, d);
    let mut b = [0.0; 3];
    let mut d = [0.0; 3];
    for i in 0..3 {
        b[i] = ekf.xhat[i];
        d[i] = ekf.xhat[i+3];
    }
    println!("EKF: b = {:?}, d = {:?}", b, d);
    for i in 0..3 {
        b[i] = ukf.xhat[i];
        d[i] = ukf.xhat[i+3];
    }
    println!("UKF: b = {:?}, d = {:?}", b, d);

    // 元データの中心と半径を計算
    let (center, radius) = sphere_fitting(&mag_samples);
    println!("--- Sphere fitting result(center & radius) of original data ---");
    println!("Center: [X0, Y0, Z0] = {:?}", center);
    println!("Radius: r = {}", radius);

    // ある程度収束したパラメータで較正して，スケールファクタとバイアスがちゃんと推定出来ているか確認
    let mut centered_calibrated_samples = Vec::<Vector3>::new();
    let mut ekf_calibrated_samples = Vec::<Vector3>::new();
    let mut ukf_calibrated_samples = Vec::<Vector3>::new();
    for mag in mag_samples {
        centered_calibrated_samples.push( centered.calibrate(mag) );
        ekf_calibrated_samples.push( ekf.calibrate(mag) );
        ukf_calibrated_samples.push( ukf.calibrate(mag) );
    }
    println!("--- Sphere fitting result(center & radius) of sequential calibrated data ---");
    println!("The center coordinates should be close to 0.");
    // Centeredアルゴリズム
    let (center, radius) = sphere_fitting(&centered_calibrated_samples);
    println!("Centered: r = {:.4}, [X0, Y0, Z0] = {:?}", radius, center);
    // EKF
    let (center, radius) = sphere_fitting(&ekf_calibrated_samples);
    println!("EKF     : r = {:.4}, [X0, Y0, Z0] = {:?}", radius, center);
    // UKF
    let (center, radius) = sphere_fitting(&ukf_calibrated_samples);
    println!("UKF     : r = {:.4}, [X0, Y0, Z0] = {:?}", radius, center);
}

/// 実験データを読み込んで，各サンプル時刻における計測値をまとめたベクトルとして返す．
fn load_csv(path: &str) -> Vec<Vector3> {
    let txt = fs::read_to_string(path).unwrap();

    let mut samples = Vec::<Vector3>::new();
    for line in txt.lines().skip(1) {  // ヘッダを読み飛ばす
        let mut tmp = Vector3::zeros();
        for (i, row) in line.split(',').enumerate() {
            tmp[i] = row.parse().unwrap();  // 文字列から数値型に変換
        }
        samples.push(tmp);
    }

    samples
}

/// 球面フィッティングにより球の中心座標と半径を求める
/// 
/// 参考文献
/// 1. j_rocket_boy, "球面フィッティングの導出と実装", 2018.
///    (https://www.slideshare.net/j_rocket_boy/fitting-88311197)
/// 
/// * Argument: 球をなす点群をまとめたベクトル[[x1,y1,z1], [x2,y2,z2], ...]
/// * Return  : (球の中心座標[x,y,z], 半径)
fn sphere_fitting(mag_samples: &Vec<Vector3>) -> ([f64; 3], f64) {
    let mut a = Matrix4x4::zeros();
    let mut b = Vector4::zeros();
    for mag in mag_samples {
        let x = mag[0];
        let y = mag[1];
        let z = mag[2];

        let x2 = x * x;
        let y2 = y * y;
        let z2 = z * z;

        a[(0, 0)] += x2;
        a[(0, 1)] += x * y;
        a[(0, 2)] += x * z;
        a[(0, 3)] += x;
        a[(1, 1)] += y2;
        a[(1, 2)] += y * z;
        a[(1, 3)] += y;
        a[(2, 2)] += z2;
        a[(2, 3)] += z;
        a[(3, 3)] += 1.0;

        let tmp = x2 + y2 + z2;
        b[0] += x * tmp;
        b[1] += y * tmp;
        b[2] += z * tmp;
        b[3] += tmp;
    }

    // 下三角行列の式は上三角行列と同じ
    for i in 1..4 {
        for j in 0..i {
            a[(i, j)] = a[(j, i)];
        }
    }
    
    b = -b;

    // Ax=Bを解く
    let x = a.lu().solve(&b).unwrap();

    let x0 = -0.5 * x[0];
    let y0 = -0.5 * x[1];
    let z0 = -0.5 * x[2];
    let r = (x0*x0 + y0*y0 + z0*z0 - x[3]).sqrt();

    ([x0, y0, z0], r)
}
