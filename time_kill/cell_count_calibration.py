from fears.utils import AutoRate
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

plate_path = 'cell_count_calibration.xlsx'
p = AutoRate.Plate(plate_path)

time_vect = np.array(p.data['Time [s]'])

row_list = ['B','C','D','E','F','G','H']
col_list = [2,3,4,5,6,7,8,9,10,11]
col_list = [str(c) for c in col_list]

row_pairs = [('B','C'),('D','E'),('F','G')]
data_avg = {}

cmap = mpl.colormaps['viridis']

condition_indx = 0

fig,ax_list = plt.subplots(nrows=3,figsize=(6,6))

for rp in row_pairs:
    ax = ax_list[condition_indx]
    col_indx = 0
    for col in col_list:
        key1 = rp[0] + col
        key2 = rp[1] + col

        ts_avg = np.array(p.data[key1] + p.data[key2])/2
        ax.plot(ts_avg,color=cmap(col_indx/10))

        data_avg[str(condition_indx)+col] = ts_avg

        col_indx += 1
    ax.set_xlim(0,400)
    condition_indx+=1

#%% Time to 200000

thresh = 200000

condition = '0'

res = []
for col in col_list:
    key = condition + col
    ts = data_avg[key]
    indx = np.argwhere(ts>=thresh)[0][0]
    res.append(time_vect[indx]/60)

dilution = []
for i in range(10):
    dilution.append(5**(-i))

fig,ax = plt.subplots()

ax.scatter(res,dilution,marker='*')

ax.set_yscale('log')

ax.set_xlabel('Time to RFU 200000 (min)')
ax.set_ylabel('Dilution factor');

dil_log = np.log10(dilution[1:])
res_fit = np.array(res[1:])*60

lin_fit = np.polyfit(res_fit,dil_log,1)

time_to_thresh_fit = np.linspace(0,600*60)
dil_fit = lin_fit[0]*time_to_thresh_fit + lin_fit[1]

ax.plot(time_to_thresh_fit/60,10**dil_fit)
# %%
