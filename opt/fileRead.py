from os import system
import numpy as np

expfile = open('expData/VLE_ACN-BZ_45.txt', 'r')
datalist = expfile.readlines()

system12 = datalist[0].rstrip('\n') # system12:系の名前
tdc = float(datalist[1].rstrip('\n')) # tdc:系の温度[C]
NVLE = int(datalist[2].rstrip('\n')) # NVLE:実験のデータ点の数

x1_list = []
g1_exp_list = []
g2_exp_list = []
for data in datalist[3::]:
    x1_list.append(float(data[0:6]))
    g1_exp_list.append(float(data[7:14]))
    g2_exp_list.append(float(data[15:22]))

# datalist_np = np.array(datalist)
print(g2_exp_list)
print(len(g1_exp_list))
expfile.close()