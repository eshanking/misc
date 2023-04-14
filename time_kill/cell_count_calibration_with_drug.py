from fears.utils import AutoRate
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import scipy

ab_plate_path = 'calibration_04132023/EK_AB_20230414_111621.xlsx'
od_plate_path = 'calibration_04132023/EK_single_OD600_20230413_120140.xlsx'
p_ab = AutoRate.Plate(ab_plate_path)
p_od = AutoRate.Plate(od_plate_path,mode='single_measurement')
od_data = p_od.od_data_to_dict(p_od.data)

#%% estimate od background

row_list = ['A','B','C','D','E','F','G','H']
bg_cols = [1,12]
bg_cols = [str(c) for c in bg_cols]

bg_est = 0
indx = 0
for row in row_list:
    for col in bg_cols:
        key = row + col
        bg_est += od_data[key]
        indx+=1

bg_est = bg_est/indx

col_list = np.arange(12) + 1
col_list = [str(c) for c in col_list]

for row in row_list:
    for col in col_list:
        key = row+col
        od_data[key] = od_data[key] - bg_est

#%% estimate cell count per column

def od_to_cells():
    

col_list = np.arange(10) + 2
row_list = ['B','C','D','E','F','G']
cell_count = []

for col in col_list:
    od_avg = 0
    for row in row_list:
        key = row+col
        od_avg += od_data[key]
    od_avg = od_avg/6

