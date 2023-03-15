#%%
import pickle
from fears.utils import AutoRate
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy.signal import butter, lfilter, freqz

data_file = 'EK_AB_20230308_125147.xlsx'

row_list = ['A','B','C','D','E','F','G','H']
col_list = [2,3,4,5,6,7,8,9,10,11]
col_list = [str(col) for col in col_list]

sample_times = [1,2,3,4,5,6,7,8,9,10,11]
sample_times = [30*st for st in sample_times]

dc = [0,10**-2,10**-1,1,10,10**2,10**3,10**3]
mic = 0.1

dc = [c*mic for c in dc]

# load normalizer object

filehandler = open('normalizer.p', 'rb') 
norm = pickle.load(filehandler)
#%%
def rolling_average(x,N):
    
    indx = 0
    res = []
    while indx+N < len(x):
        xt = x[indx:indx+N]
        res.append(np.nanmean(xt))
        indx+=1
    
    for i in range(len(x[indx:])):
        res.append(np.nanmean(x[indx+i:]))

    res = np.array(res)
    x = np.array(x)
    res[x == 0] = 0

    return res

def butter_lowpass(cutoff, fs, order=5):
    return butter(order, cutoff, fs=fs, btype='low', analog=False)

def butter_lowpass_filter(data, cutoff, fs, order=5):
    b, a = butter_lowpass(cutoff, fs, order=order)
    y = lfilter(b, a, data)
    return y

order = 1
fs = 0.01      # sample rate, Hz
cutoff = 0.0001  # desired cutoff frequency of the filter, Hz

# %%
p = AutoRate.Plate(data_file)
time_vect = np.array(p.data['Time [s]'])

fig,ax_list = plt.subplots(nrows=4,ncols=2,figsize=(8,10))
ax_list_t = ax_list.reshape(-1)
cmap = mpl.colormaps['viridis']

norm_data = {}

row_indx = 0
for row in row_list:
    ax = ax_list_t[row_indx]
    col_indx = 0
    max_fluor = 0
    for col in col_list:
        key = row+col
        ts = np.array(p.data[key])
        conc = dc[row_indx]
        if conc>0:
            conc = np.log10(conc)
            for indx in range(len(ts)):
                if time_vect[indx] >= np.max(np.array(norm.data['Time [s]'])):
                    time_t = np.max(np.array(norm.data['Time [s]']))
                else:
                    time_t = time_vect[indx]
                norm_factor = norm.get_ctx_norm_factor(conc,time_t)
                ts[indx] = ts[indx]/norm_factor
        ts = butter_lowpass_filter(ts, cutoff, fs, order)
        norm_data[key] = ts

        if np.max(ts) > max_fluor:
            max_fluor = np.max(ts)

        ax.plot(time_vect,ts,color=cmap(col_indx/10))
        # ax.set_xlim(0,20000)
        col_indx+=1
    xl = ax.get_xlim()
    max_fluor = max_fluor/4
    ax.plot([xl[0],xl[1]],[max_fluor,max_fluor],color='black',linewidth=2)

    row_indx+=1

# %%
