import numpy as np
import random

# データの読み込み#####################################################
expfile = open('expData/VLE_ACN-BZ_45.txt', 'r')
datalist = expfile.readlines()

system12 = datalist[0].rstrip('\n') # system12:系の名前
tdc = float(datalist[1].rstrip('\n')) # tdc:系の温度[C]
NVLE = int(datalist[2].rstrip('\n')) # exp_point_num:実験のデータ点の数

x1_list = []
g1_exp_list = []
g2_exp_list = []
for data in datalist[3::]:
  x1_list.append(float(data[0:6]))
  g1_exp_list.append(float(data[7:14]))
  g2_exp_list.append(float(data[15:22]))

expfile.close()
####################################################################


# 評価関数(目的関数)
def f(g_sum, NVLE):
  z = g_sum * 100 / (2 * NVLE) 
  return z

# 目的関数に必要な計算
# c21i, c12i, alpの最適解を探索する
############for i in range(NVLE)で回す範囲######################################################################
def g_calc(c21ix, c12ix, alpx, x1):
  c21p = alpx * c21ix
  c12p = alpx * c12ix
  x2 = 1.0 - x1
  c1 = c21p * x1 + c21ix * x2
  c2 = c12p * x2 + c12ix * x2
  asa1 = c1 * x1 / (c1 * x1 + c2 * x2) # 平均表面積
  asa2 = 1.0 - asa1 # 平均表面積
  g1_calc = np.exp((c21p + (c21ix - c21p) * x2 ** 2) * asa2 ** 2 + (c12ix - c12p) * x2 ** 2 * asa1 ** 2)
  g2_calc = np.exp((c12p + (c12ix - c12p) * x1 ** 2) * asa1 ** 2 + (c21ix - c21p) * x1 ** 2 * asa2 ** 2)
  return g1_calc, g2_calc

def g1_g2(g1_exp, g2_exp, g1_calc, g2_calc):
  g1g2 = abs(g1_exp - g1_calc) / g1_exp + abs(g2_exp - g2_calc) / g2_exp
  return g1g2
###############################################################################################################

# 粒子の位置の更新を行う関数
def update_position(c21ix, c12ix, alpx, c21iv, c12iv, alpv):
  new_c21ix = c21ix + c21iv
  new_c12ix = c12ix + c12iv
  new_alpx = alpx + alpv
  return new_c21ix, new_c12ix, new_alpx

# 粒子の速度の更新を行う関数
def update_velocity(c21ix, c12ix, alpx, c21iv, c12iv, alpv, p, g, I, rho_max=0.5):
  # パラメーターrhoはランダムに与える
  rho1 = random.uniform(0, rho_max)
  rho2 = random.uniform(0, rho_max)
  # 粒子速度の更新を行う
  new_c21iv = I * c21iv + rho1 * (p["c21ix"] - c21ix) + rho2 * (g["c21ix"] - c21ix)
  new_c12iv = I * c12iv + rho1 * (p["c12ix"] - c12ix) + rho2 * (g["c12ix"] - c12ix)
  new_alpv = I * alpv + rho1 * (p["alpx"] - alpx) + rho2 * (g["alpx"] - alpx)
  return new_c21iv, new_c12iv, new_alpv

N = 250  # 粒子の数
c21ix_min, c21ix_max = 0, 2
c12ix_min, c12ix_max = 0, 2
alpx_min, alpx_max = 0, 2

# 粒子位置, 速度, パーソナルベスト, グローバルベストの初期化を行う
ps = [{"c21ix": random.uniform(c21ix_min, c21ix_max), "c12ix": random.uniform(c12ix_min, c12ix_max), "alpx": random.uniform(alpx_min, alpx_max)} for i in range(N)]
vs = [{"c21ix": 0.0, "c12ix": 0.0, "alpx": 0.0} for i in range(N)]
personal_best_positions = ps

##################### for文で計算###########################################
personal_best_scores = []
for p in ps:
  g_sum = 0
  for n in range(NVLE):
    x1 = x1_list[n]
    g1_exp = g1_exp_list[n]
    g2_exp = g2_exp_list[n]
    g1_calc, g2_calc = g_calc(p["c21ix"], p["c12ix"], p["alpx"], x1)
    g1g2 = g1_g2(g1_exp, g2_exp, g1_calc, g2_calc)
    g_sum += g1g2
  personal_best_scores.append(f(g_sum, NVLE))
###########################################################################

best_particle = np.argmin(personal_best_scores)
global_best_position = personal_best_positions[best_particle]
# print(personal_best_positions)
# print(personal_best_scores)
# print(best_particle)
# print(global_best_position)

T = 1000  # 制限時間(ループの回数)
for t in range(T):
    for i in range(N):
        c21ix, c12ix, alpx = ps[i]["c21ix"], ps[i]["c12ix"], ps[i]["alpx"]
        c21iv, c12iv, alpv = vs[i]["c21ix"], vs[i]["c12ix"], vs[i]["alpx"]
        p = personal_best_positions[i]
        g = global_best_position
        # 粒子の位置の更新を行う
        new_c21ix, new_c12ix, new_alpx = update_position(c21ix, c12ix, alpx, c21iv, c12iv, alpv)
        ps[i] = {"c21ix": new_c21ix, "c12ix": new_c12ix, "alpx": new_alpx}
        # 粒子の速度の更新を行う
        I = 0.7 - (0.7 - 0.5) / T * t
        new_c21iv, new_c12iv, new_alpv = update_velocity(new_c21ix, new_c12ix, new_alpx, c21iv, c12iv, alpv, p, g, I)
        vs[i] = {"c21ix": new_c21iv, "c12ix": new_c12iv, "alpx": new_alpv}
        # 評価値を求め, パーソナルベストの更新を行う
        new_g_sum = 0
        for n in range(NVLE):
          x1 = x1_list[n]
          g1_exp = g1_exp_list[n]
          g2_exp = g2_exp_list[n]
          new_g1_calc, new_g2_calc = g_calc(new_c21ix, new_c12ix, new_alpx, x1)
          new_g1g2 = g1_g2(g1_exp, g2_exp, new_g1_calc, new_g2_calc)
          new_g_sum += new_g1g2
        score = f(new_g_sum, NVLE)
        if score < personal_best_scores[i]:
            personal_best_scores[i] = score
            personal_best_positions[i] = {"c21ix": new_c21ix, "c12ix": new_c12ix, "alpx": new_alpx}    
    # グローバルベストの更新を行う
    best_particle = np.argmin(personal_best_scores)
    global_best_position = personal_best_positions[best_particle]
    # print(min(personal_best_scores))
   
# 最適解
c21px = global_best_position["c21ix"] * global_best_position["alpx"]
c12px = global_best_position["c12ix"] * global_best_position["alpx"]
dE = (global_best_position["c21ix"] + global_best_position["c12ix"]) / 2

print(system12)
print(tdc)
print("c21px   c12px   c21ix   c12ix   dE")
print(str('{:.5f}'.format(c21px)) + '\t' + str('{:.5f}'.format(c12px)) + '\t' + str('{:.5f}'.format(global_best_position["c21ix"])) + '\t' + str('{:.5f}'.format(global_best_position["c12ix"])) + '\t' + str('{:.5f}'.format(dE)))
print(min(personal_best_scores))