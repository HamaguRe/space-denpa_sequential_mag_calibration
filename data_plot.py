# result_*.csvを読み込んでグラフ作成

import csv
import matplotlib.pyplot as plt


# CSVファイルを読んでヘッダと値を対応させた辞書型を返す
#
# --- CSV File ---
# TimeStep,x,y
# 0,0.01,-0.1
# 1,0.02,-0.2
# 2,0.03,-0.3
# ---- Return -----
# {'TimeStep': [0,1,2], 'x': [0.01,0.02,0.03], 'y': [-0.1,-0.2,-0.3]}
def load_csv(file_path):
    data_dict = {}
    with open(file_path, 'r', newline='', encoding='utf-8') as f:
        reader = csv.DictReader(f)
        for row in reader:
            for key, value in row.items():
                if key not in data_dict:
                    data_dict[key] = []
                data_dict[key].append(float(value))

    return data_dict


# 較正前のデータを読み込む（比較表示用）
mag_raw = load_csv("./magnetic_data/data-set_1.csv")

# 較正結果のデータを読み込む
result_centered = load_csv("./result_centered.csv")
result_ekf = load_csv("./result_ekf.csv")
result_ukf = load_csv("./result_ukf.csv")

# スケールファクタとバイアスの時間遷移を表示するウィンドウ
fig1 = plt.figure(figsize = (17, 9))
fig1.suptitle("Time evolution of estimated parameters", fontsize=16)

centered_ax1 = fig1.add_subplot(2, 3, 1, ylabel="Bias [Gauss]")
centered_ax2 = fig1.add_subplot(2, 3, 4, xlabel="Time step", ylabel="Scale factor")
ekf_ax1 = fig1.add_subplot(2, 3, 2, ylabel="Bias [Gauss]")
ekf_ax2 = fig1.add_subplot(2, 3, 5, xlabel="Time step", ylabel="Scale factor")
ukf_ax1 = fig1.add_subplot(2, 3, 3, ylabel="Bias [Gauss]")
ukf_ax2 = fig1.add_subplot(2, 3, 6, xlabel="Time step", ylabel="Scale factor")
centered_ax1.set_title("Centered")
ekf_ax1.set_title("EKF")
ukf_ax1.set_title("UKF")

# --- 描画 --- 
k = result_centered["TimeStep"]  # 計算ステップ

# バイアス
labels = ['x', 'y', 'z']
for lb in labels:
    centered_ax1.plot(k, result_centered["b_" + lb], label=lb)
    ekf_ax1.plot(k, result_ekf["b_" + lb], label=lb)
    ukf_ax1.plot(k, result_ukf["b_" + lb], label=lb)

centered_ax1.set_xlim(0, k[-1])
ekf_ax1.set_xlim(0, k[-1])
ukf_ax1.set_xlim(0, k[-1])
lim_bias = [-0.6, 0.6]
centered_ax1.set_ylim(lim_bias[0], lim_bias[1])
ekf_ax1.set_ylim(lim_bias[0], lim_bias[1])
ukf_ax1.set_ylim(lim_bias[0], lim_bias[1])

# スケールファクタ
labels = ["D11", "D22", "D33"]
for lb in labels:
    centered_ax2.plot(k, result_centered[lb], label=lb)
    ekf_ax2.plot(k, result_ekf[lb], label=lb)
    ukf_ax2.plot(k, result_ukf[lb], label=lb)

centered_ax2.set_xlim(0, k[-1])
ekf_ax2.set_xlim(0, k[-1])
ukf_ax2.set_xlim(0, k[-1])

centered_ax1.legend(loc="lower right")
centered_ax2.legend(loc="lower right")
ekf_ax1.legend()
ekf_ax2.legend()
ukf_ax1.legend()
ukf_ax2.legend()

# 3Dプロットしない場合はコメントアウト（デバッグ時など）
#plt.show(); quit()

# --- 元データと較正したデータを3Dプロットで比較するウィンドウ
fig2 = plt.figure(figsize = (12, 6))
#fig2.suptitle("Comparison of raw data and calibrated data", fontsize=16)

centered_ax3 = fig2.add_subplot(1, 3, 1, projection="3d")
ekf_ax3 = fig2.add_subplot(1, 3, 2, projection="3d")
ukf_ax3 = fig2.add_subplot(1, 3, 3, projection="3d")
centered_ax3.set_title("Centered")
ekf_ax3.set_title("EKF")
ukf_ax3.set_title("UKF")

# 座標系の原点
centered_ax3.scatter([0],[0],[0],marker="o", color="r", label="Origin of coordinate system")
ekf_ax3.scatter([0],[0],[0],marker="o", color="r", label="Origin of coordinate system")
ukf_ax3.scatter([0],[0],[0],marker="o", color="r", label="Origin of coordinate system")

centered_ax3.scatter(mag_raw["x"], mag_raw["y"], mag_raw["z"], s=5, label="Raw data")
centered_ax3.scatter(result_centered["calib_x"], result_centered["calib_y"], result_centered["calib_z"], s=5, color="darkorange", label="Sequential calibrated data")
ekf_ax3.scatter(mag_raw["x"], mag_raw["y"], mag_raw["z"], s=5, label="Raw data")
ekf_ax3.scatter(result_ekf["calib_x"], result_ekf["calib_y"], result_ekf["calib_z"], s=5, color="darkorange", label="Sequential calibrated data")
ukf_ax3.scatter(mag_raw["x"], mag_raw["y"], mag_raw["z"], s=5, label="Raw data")
ukf_ax3.scatter(result_ukf["calib_x"], result_ukf["calib_y"], result_ukf["calib_z"], s=5, color="darkorange", label="Sequential calibrated data")

lim = 1.0
centered_ax3.set_xlim(-lim, lim)
centered_ax3.set_ylim(-lim, lim)
centered_ax3.set_zlim(-lim, lim)
centered_ax3.set_xlabel("X-axis")
centered_ax3.set_ylabel("Y-axis")
centered_ax3.set_zlabel("Z-axis")
ekf_ax3.set_xlim(-lim, lim)
ekf_ax3.set_ylim(-lim, lim)
ekf_ax3.set_zlim(-lim, lim)
ekf_ax3.set_xlabel("X-axis")
ekf_ax3.set_ylabel("Y-axis")
ekf_ax3.set_zlabel("Z-axis")
ukf_ax3.set_xlim(-lim, lim)
ukf_ax3.set_ylim(-lim, lim)
ukf_ax3.set_zlim(-lim, lim)
ukf_ax3.set_xlabel("X-axis")
ukf_ax3.set_ylabel("Y-axis")
ukf_ax3.set_zlabel("Z-axis")

ukf_ax3.legend(loc="upper right")

plt.show()