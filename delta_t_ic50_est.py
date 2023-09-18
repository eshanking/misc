from fears.utils import AutoRate
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import butter, lfilter

def rolling_average(x,N):
    """Rolling average function

    Args:
        x (array-like): array to be smoothed 
        N (int): window size

    Returns:
        list: smoothed array
    """
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

data_file = 'EK_OD600_24hr_20230916_180536.xlsx'

plate = AutoRate.Plate(data_file)

data = plate.data
# data = plate.od_data_to_dict(data)

cols = list(np.arange(12)+1)
cols = [str(c) for c in cols]
rows = ['A','B','C','D','E','F','G','H']

fig,ax_list = plt.subplots(nrows = 8,ncols=6,figsize=(10,10),sharex=True,sharey=True)

t = np.array(data['Time [s]'])

bg_time = 1000
bg_indx = np.argwhere(t >= bg_time)[0][0]

order = 1
fs = 1/(332)    # sample rate, Hz
cutoff = fs/10  # desired cutoff frequency of the filter, Hz

row_indx = 0
for r in rows:
   
    col_indx = 0
    for c in cols[0:6]:
        ax = ax_list[row_indx,col_indx]
        key = r + c
        ts = np.array(data[key])
        # ts = rolling_average(ts,5)
        ts = butter_lowpass_filter(ts,cutoff,fs,order)
        bg = np.mean(ts[0:bg_indx])
        # ts = ts - bg
        # ts = ts - np.min(ts)
        # if col_indx < 6:
        ax.plot(t,ts,color='black')
        ax.set_title(key)
        # else:
        #     ax.plot(t,ts,'--',color='black')
        col_indx+=1
    row_indx += 1
    # ax.set_yscale('log')
    ax.set_ylim(0,0.75)

fig.tight_layout()
