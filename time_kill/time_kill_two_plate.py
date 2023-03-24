#%%
from fears.utils import AutoRate
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy.signal import butter, lfilter
import scipy

data_file = 'EK_AB_20230324_090141.xlsx'

row_list = ['B','C','D','E','F','G']
col_list = [2,3,4,5,6,7,8,9,10,11]
col_list = [str(col) for col in col_list]

sample_times = [0,30,90,210,390]
# 11:12:44
st_data_cleaning = [32*60-13,93*60+7,213*60+-18,397*60-26]

dc = [0,10**-1,10**0,10,10**2,10**3]
mic = 0.1

dc = [c*mic for c in dc]

#%% Helper functions
def est_time_kill(xdata,ydata,debug=False):
    # alpha,A,l
    
    indx = np.argwhere(np.abs(ydata)==np.max(np.abs(ydata)))[0][0]
    A_est = ydata[indx]

    alpha_est = 0
    p0 = [alpha_est,A_est,0]
    # bounds = [[-np.inf,A_est-0.2*(np.abs(A_est)),-np.inf],
    #     [np.inf,A_est+0.2*(np.abs(A_est)),np.inf]]

    popt,pcov = scipy.optimize.curve_fit(time_kill_model,xdata,ydata,p0=p0,maxfev=10000)

    if debug:

        y_opt = time_kill_model(xdata,popt[0],popt[1],popt[2])
        
        fig,ax = plt.subplots()
        ax.plot(xdata,y_opt,label='est')
        ax.scatter(xdata,ydata,marker="*",label='data')
        ax.set_title(popt[0])
        ax.legend()
        # ax.set_ylim(-1000,50000)

    return popt

def net_growth(t,N0,K,Kss,alpha):

    g = []

    for tt in t:
        g.append(N0 + (1/2.303)*((K-Kss)*tt + (Kss/alpha)*(1-np.exp(-alpha*tt))))

    return g

def est_net_growth(xdata,ydata,K=None,N0=None,debug=False):
    
    if K is None: # want to estimate K
        p0 = [ydata[0],0.01,0,0.1] # N0, K, Kss, alpha
        bounds = [[ydata[0]-0.1*ydata[0],0,-10**-5,0.1],
                  [ydata[0]+0.1*ydata[0],np.inf,np.inf,20]]
        popt,pcov = scipy.optimize.curve_fit(net_growth,
                                    xdata,ydata,p0 = p0,bounds=bounds,maxfev=10000)
    else:
        p0 = [0,10**-5] # Kss, alpha
        popt,pcov = scipy.optimize.curve_fit(lambda t,Kss,alpha: net_growth(t,N0,K,Kss,alpha),
                                    xdata,ydata,p0 = p0)
    
    if debug:

        if K is None:
            Kss = popt[2]
            y_opt = net_growth(xdata,popt[0],popt[1],popt[2],popt[3])
        else:
            Kss = popt[0]
            y_opt = net_growth(xdata,N0,K,popt[0],popt[1])
        
        fig,ax = plt.subplots()
        ax.plot(xdata,y_opt,label='est')
        ax.scatter(xdata,ydata,marker="*",label='data')
        ax.set_title(round(Kss,3))
        ax.legend()

    return popt
    

def time_kill_model(t,alpha,A,l):
    res = []
    for tt in t:
        res.append(A*(1-np.exp(-alpha*(tt-l)/A)))
    return res

# def time_kill_model(t,kss,alpha,y0):
#     # Growth rate constant of control arbitrarily set to zero
#     res = []
#     for tt in t:
#         res.append(y0 - kss*tt + (kss/alpha)*(1-np.exp(-1*alpha*tt)))
#     return res

def butter_lowpass(cutoff, fs, order=5):
    return butter(order, cutoff, fs=fs, btype='low', analog=False)

def butter_lowpass_filter(data, cutoff, fs, order=5):
    b, a = butter_lowpass(cutoff, fs, order=order)
    y = lfilter(b, a, data)
    return y

order = 1
fs = 0.01      # sample rate, Hz
cutoff = 0.0001  # desired cutoff frequency of the filter, Hz

def est_pharm_curve(xdata,ydata,debug=False):


    # gmax,gmin,mic,k

    g_max_est = ydata[0]
    gmin_est = np.min(ydata)

    mic_est = np.argwhere(np.array(ydata)<=0)[0][0]
    mic_est = xdata[mic_est]

    print(mic_est)

    p0 = [g_max_est,gmin_est,mic_est,0.1]
    bounds = [[g_max_est-0.1*g_max_est,gmin_est+0.1*gmin_est,mic_est/10,0.1],
              [g_max_est+0.1*g_max_est,gmin_est-0.1*gmin_est,mic_est*10,10]]

    # p0 = [1,g_drugless_est,-1,0]
    # bounds = ([-np.inf,g_drugless_est-g_drugless_est*0.05,-2,0],
    #           [np.inf,g_drugless_est+g_drugless_est*0.05,-0.1,np.inf])

    popt,pcov = scipy.optimize.curve_fit(pharmacodynamic_curve,xdata,ydata,p0=p0,
                                         bounds=bounds)

    if debug:

        # xmin = np.log10(xdata[1])
        # xmax = np.log10(xdata[-1])
        dc_fit = np.logspace(-3,3)
        dc_fit[0] = 0
        y_opt = pharmacodynamic_curve(dc_fit,popt[0],popt[1],popt[2],popt[3])
        fig,ax = plt.subplots()
        ax.plot(dc_fit,y_opt)
        ax.plot(xdata,ydata,'o')

        ax.set_xscale('log')

    return popt

def pharmacodynamic_curve(c, gmax, gmin, mic, k):
    """pharmacodynamic model adapted from Foerster et al.

    Foerster, S., Unemo, M., Hathaway, L.J. et al. Time-kill curve analysis and
    pharmacodynamic modelling for in vitro evaluation of antimicrobials against Neisseria
    gonorrhoeae . BMC Microbiol 16, 216 (2016). 
    https://doi.org/10.1186/s12866-016-0838-9

    Args:
        c (float): drug concentration
        gmax (float): max growth rate
        gmin (float): min growth rate
        mic (float): estimated minimum inhibitory concentration
        k (float): hill coefficient
    """
    g = []
    for c_t in c:
        if c_t == 0:
            g.append(gmax)
        else:
            g.append(gmax - (((gmax-gmin)*(c_t/mic)**k)/((c_t/mic)**k-(gmin/gmax))))
    
    return g

#%%

col_pairs = [(2,3), (4,5), (6,7), (8,9), (10,11)]
data_avg = {}

p = AutoRate.Plate(data_file)
time_vect = np.array(p.data['Time [s]'])

data_cleaning_indx = []
for tt in st_data_cleaning:
    indx = np.argwhere(time_vect>=tt)[0][0]
    p.data = p.data.drop(p.data.index[indx-1:indx+2])
    data_cleaning_indx.append(indx)
    time_vect = np.array(p.data['Time [s]'])

p.data = p.data.drop(p.data.index[-1]) # remove last datapoint because of incomplete data
time_vect = np.array(p.data['Time [s]'])

fig,ax_list = plt.subplots(nrows=3,ncols=2,figsize=(8,10))
ax_list_t = ax_list.reshape(-1)
cmap = mpl.colormaps['viridis']

row_indx = 0
for row in row_list:
    ax = ax_list_t[row_indx]
    sample_indx = 0

    for cp in col_pairs:
        key1 = row + str(cp[0])
        key2 = row + str(cp[1])

        ts1 = np.array(p.data[key1])
        ts1[ts1=='OVER'] = np.nan

        ts2 = np.array(p.data[key2])
        ts2[ts2=='OVER'] = np.nan
        
        ts_avg = np.nanmean((ts1,ts2),axis=0)

        data_avg[row+str(sample_indx)] = ts_avg
        si = sample_times[sample_indx]*60

        ts1 = butter_lowpass_filter(ts1, cutoff, fs, order)
        ts2 = butter_lowpass_filter(ts2, cutoff, fs, order)
        # ax.plot(time_vect-si,ts_avg,color=cmap(sample_indx/5))
        ax.plot(time_vect-si,ts1,color=cmap(sample_indx/5))
        ax.plot(time_vect-si,ts2,'--',color=cmap(sample_indx/5))
        ax.set_ylim(0,200000)
        ax.set_xlim(0,60000)
        sample_indx+=1

    row_indx+=1
# %% time to 200,000 

thresh = 200000
res = {}

samples = ['0','1','2','3','4']

fig,ax_list = plt.subplots(nrows=3,ncols=2,figsize=(8,10))
ax_list_t = ax_list.reshape(-1)

row_indx = 0
for row in row_list:

    time_to_thresh = []

    sample_indx = 0

    for s in samples:
        key = row + s

        ts = data_avg[key]
        
        indx = np.argwhere(ts>=thresh)[0][0]

        time_t = time_vect[indx] - sample_times[sample_indx]*60

        time_to_thresh.append(time_t)
        sample_indx +=1

    ax = ax_list_t[row_indx]
    ax.plot(sample_times,time_to_thresh)
    res[str(row_indx)] = time_to_thresh

    row_indx+=1

# %%

growth_rates = []

fig,ax = plt.subplots()
# ax_list_t = ax_list.reshape(-1)

y0 = res['0'][0]
cell_count = []

indx = 0

for key in res.keys():

    ydata = np.array(res[key])

    ydata[0] = y0

    # convert to cell count
    ydata = (-1.7*10**-4)*ydata - 2.96*10**-1
    # ydata = 10**ydata

    cell_count.append(ydata + 9)

    st = sample_times

    if key == '4':
        ydata = ydata[0:-1]
        st = st[0:-1]

    if key == '3':
        ydata = ydata[0:-1]
        st = st[0:-1]

    ax.plot(st,ydata,color=cmap(indx/5),label=round(dc[indx],2))

    ydata = ydata - ydata[0]

    indx+=1

    popt = est_time_kill(st,ydata,debug=False)

    r = popt[0]
    A = popt[1]


    growth_rates.append(r)

ax.set_xlabel('Time (min)',fontsize=14)
ax.set_ylabel('Log proportion of carrying capacity',fontsize=14)
ax.legend(frameon=False)
# %%
fig,ax = plt.subplots()
dc_plot = [-2,-1,0,1,2,3]
ax.scatter(dc_plot,growth_rates)
# ax.set_xscale('log')
# %%

popt = est_pharm_curve(dc,growth_rates,debug=True)

fig,ax = plt.subplots()

dc_fit = np.linspace(-3,3,num=100)
dc_fit = 10**dc_fit
dc_fit[0] = 0
g_fit = pharmacodynamic_curve(dc_fit,popt[0],popt[1],popt[2],1)
dc_fit_plot = np.linspace(-3,3,num=100)

ax.plot(dc_fit_plot,g_fit)
# %%

cell_count_trunc = []

for i in range(len(cell_count)):
    
    if i in [0,1,3,4]:
        cell_count_trunc.append(cell_count[i][0:-1])
    elif i == 5:
        cell_count_trunc.append(cell_count[i][0:-2])
    else:
        cell_count_trunc.append(cell_count[i])

net_killing_rate = []
alpha_est = []

# estimate K
st = sample_times[0:len(cell_count_trunc[0])]
popt = est_net_growth(st,cell_count_trunc[0],debug=True)
K = popt[1]
N0 = popt[0]

net_killing_rate.append(popt[2])
alpha_est.append(popt[-1])

for cc in cell_count_trunc[1:]:
    st = sample_times[0:len(cc)]
    popt = est_net_growth(st,cc,K=K,N0=N0,debug=True)
    net_killing_rate.append(popt[0])
    alpha_est.append(popt[-1])


# %%
