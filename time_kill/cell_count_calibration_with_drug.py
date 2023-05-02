from fears.utils import AutoRate
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy import stats

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

def od_to_cells(od):
    res = [297761.03865714, 324941.3815491, 17609.09483884]
    return res[2] + res[1]*od + res[0]*od**2

col_list = np.arange(10) + 2
row_list = ['B','C','D','E','F','G']
cell_count = []
od = []

for col in col_list:
    od_avg = 0
    for row in row_list:
        key = row + str(col)
        od_avg += od_data[key]
    od_avg = od_avg/6
    cell_count.append(od_to_cells(od_avg))
    od.append(od_avg)

cell_count_est = cell_count

cell_count_col_1 = cell_count[0]
cell_count = []

for i in range(10):
    cell_count.append(cell_count_col_1/(2**i))

# %%
time_vect = np.array(p_ab.data['Time [s]'])

fig,ax_list = plt.subplots(nrows=2)
cmap = mpl.colormaps['viridis']

col_indx = 0

for col in col_list[1:]:
    ax = ax_list[0] 
    ts_avg = np.zeros(len(time_vect))
    for row in ['B','C','D']:
        key = row + str(col)
        ts_avg += np.array(p_ab.data[key]).astype('float64')
    ts_avg = ts_avg/3
    ax.plot(time_vect,ts_avg,color=cmap(col_indx/8))
    # col_indx += 1
    ax = ax_list[1] 
    ts_avg = np.zeros(len(time_vect))
    for row in ['E','F','G']:
        key = row + str(col)
        ts_avg += np.array(p_ab.data[key]).astype('float64')
    ts_avg = ts_avg/3
    ax.plot(time_vect,ts_avg,color=cmap(col_indx/8))
    col_indx += 1

for ax in ax_list:
    ax.set_xlim(0,3600)
    ax.set_ylim(0,400000)
# %%
fig,ax_list = plt.subplots(nrows=4,figsize=(4,8))
cmap = mpl.colormaps['viridis']

dilution = []
for i in range(len(col_list)):
    dilution.append(str(2**i) + 'x')

col_indx = 0

for col in col_list[1:]:
    ax = ax_list[0] 
    ts_avg = np.zeros(len(time_vect))
    for row in ['B','C','D']:
        key = row + str(col)
        ts_avg += np.array(p_ab.data[key]).astype('float64')
    ts_avg = ts_avg/3
    ax.plot(time_vect/60,ts_avg,color=cmap(col_indx/8))

    row_indx = 1
    for row in ['E','F','G']:
        ax = ax_list[row_indx]
        key = row + str(col)
        ts = np.array(p_ab.data[key]).astype('float64')
        ax.plot(time_vect/60,ts,color=cmap(col_indx/8),label=dilution[col_indx])
        row_indx+=1
    col_indx+=1

ax_list[0].set_title('No drug')
ax_list[1].set_title('1x MIC')
ax_list[2].set_title('10x MIC')
ax_list[3].set_title('100x MIC')

for ax in ax_list:
    ax.set_xlim(0,60)
    ax.set_ylim(0,500000)
    ax.ticklabel_format(axis='y',style='sci',scilimits=(0,0))
    ax.set_xlabel('Time (min)')
    ax.set_ylabel('RFU (au)')

ax_list[3].legend(ncol=4,frameon=False,loc=(-0.05,-1.35))

fig.tight_layout()
# %%
# fig,ax_list = plt.subplots(nrows=4,figsize=(4,8))
cmap = mpl.colormaps['viridis']

dilution = []
for i in range(len(col_list)):
    dilution.append(1/(2**i))

hr1_time_indx = np.argwhere(time_vect>=3600)[0][0]

col_indx = 0

no_drug_fluor = []

for col in col_list[1:]:
    ts_avg = np.zeros(len(time_vect))
    for row in ['B','C','D']:
        key = row + str(col)
        ts_avg += np.array(p_ab.data[key]).astype('float64')
    ts_avg = ts_avg/3
    no_drug_fluor.append(ts_avg[hr1_time_indx])

#1xMIC

MICx1_fluor = []

for col in col_list[1:]:

    key = 'E' + str(col)
    ts = np.array(p_ab.data[key]).astype('float64')

    MICx1_fluor.append(ts[hr1_time_indx])

#10xMIC

MICx10_fluor = []

for col in col_list[1:]:

    key = 'F' + str(col)
    ts = np.array(p_ab.data[key]).astype('float64')

    MICx10_fluor.append(ts[hr1_time_indx])

#100xMIC

MICx100_fluor = []

for col in col_list[1:]:

    key = 'G' + str(col)
    ts = np.array(p_ab.data[key]).astype('float64')

    MICx100_fluor.append(ts[hr1_time_indx])

fig,ax = plt.subplots()

dilution_log = np.arange(10)
ax.plot(dilution_log[1:],no_drug_fluor,color=cmap(0),linewidth=2,label='no drug')
ax.plot(dilution_log[1:],MICx1_fluor,color=cmap(0.33),linewidth=2,label='1xMIC')
ax.plot(dilution_log[1:],MICx10_fluor,color=cmap(0.66),linewidth=2,label='10xMIC')
ax.plot(dilution_log[1:],MICx100_fluor,color=cmap(0.99),linewidth=2,label='100xMIC')

dilution_xlabel = ['$2^{' + str(d) + '}$' for d in dilution_log]
ax.set_xticks(np.arange(1,10))
ax.set_xticklabels(dilution_xlabel[1:])

ax.ticklabel_format(axis='y',style='sci',scilimits=(0,0))
ax.tick_params(axis='both',labelsize=12)

ax.set_ylabel('RFU (a.u.)',fontsize=14)
ax.set_xlabel('Dilution',fontsize=14)

ax.legend(frameon=False,fontsize=12)
fig.savefig('cell_count_vs_dilution.png',bbox_inches='tight')
# %% linear fit

ydata = np.mean((no_drug_fluor,MICx100_fluor,MICx10_fluor,MICx1_fluor),axis=0)

res = stats.linregress(dilution_log[1:],ydata)

fig,ax = plt.subplots()

ax.scatter(dilution_log[1:],no_drug_fluor,color=cmap(0),label='no drug')
ax.scatter(dilution_log[1:],MICx1_fluor,color=cmap(0.33),label='1xMIC')
ax.scatter(dilution_log[1:],MICx10_fluor,color=cmap(0.66),label='10xMIC')
ax.scatter(dilution_log[1:],MICx100_fluor,color=cmap(0.99),label='100xMIC')
ax.plot(dilution_log[1:], res.intercept + res.slope*dilution_log[1:], 'r', label='linear fit')

dilution_xlabel = ['$2^{' + str(d) + '}$' for d in dilution_log]
ax.set_xticks(np.arange(1,10))
ax.set_xticklabels(dilution_xlabel[1:])

ax.ticklabel_format(axis='y',style='sci',scilimits=(0,0))
ax.tick_params(axis='both',labelsize=12)

ax.set_ylabel('RFU (a.u.)',fontsize=14)
ax.set_xlabel('Dilution',fontsize=14)

ax.annotate('$R^{2}$ = ' + str(round(res.rvalue**2,3)),(1,50000),fontsize=12)

fig.savefig('cell_count_lin_regress.png',bbox_inches='tight')
# %% RFU to cell count

res = stats.linregress(ydata,dilution_log[1:])

def rfu_to_cell_count(rfu):
    dilution_factor = 8.66 + (-1.812*10**-5)*rfu
    dilution_factor=2**(-dilution_factor)
    return dilution_factor*90000
# %%
