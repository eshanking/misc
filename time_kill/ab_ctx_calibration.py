#%%
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.animation as animation
from fears.utils import AutoRate
import scipy

def get_polyfit_result(x,res):
    y = 0
    n = len(res)-1
    indx = 0
    for r in res:
        y += r*x**(n-indx)
        indx+=1
    return y

data_file = 'EK_AB_20230310_101015.xlsx'

p = AutoRate.Plate(data_file)

final_plate = np.zeros((8,10))

row_list = ['A','B','C','D','E','F','G','H']

for row_indx in range(8):
    for col in range(1,11):
        key = row_list[row_indx] + str(col+1)
        ts = np.array(p.data[key])
        final_plate[row_indx,col-1] = ts[-1]

# final_plate = final_plate - np.min(final_plate)
final_plate = final_plate/np.min(final_plate)

fig,ax = plt.subplots()
im = ax.imshow(final_plate)
fig.colorbar(im,ax=ax)
#%%

fig,ax = plt.subplots()

# dc = [1000,100,10,1,0.1,0.01,0.001,0]
dc = [0,0.001,0.01,0.1,1,10,100,1000]

bp = ax.boxplot(final_plate.T,labels=dc);
ax.set_ylabel('Normalized fluorescence',fontsize=14)
ax.set_xlabel('Drug concentration (ug/mL)',fontsize=14)

dc_log = [-4,-3,-2,-1,0,1,2,3]

mean_fluor = np.mean(final_plate,axis=1)

res = np.polyfit(dc_log,mean_fluor,3)

dc_fit = np.linspace(-4,3,num=100)

yfit = get_polyfit_result(dc_fit,res)

xt = ax.get_xticks()
dc_fit_plot = np.linspace(xt[0],xt[-1],num=100)

ax.plot(dc_fit_plot,yfit,color='red')
# %%
cmap = mpl.colormaps['viridis']
fig,ax_list = plt.subplots(nrows=10,figsize=(3,8))

for row_indx in range(8):
    for col in range(1,11):
        key = row_list[row_indx] + str(col+1)
        ts = np.array(p.data[key])
        ax_list[col-1].plot(ts,color=cmap(row_indx/8))
# %%

spline_fit = scipy.interpolate.CubicSpline(dc_log,mean_fluor)

fig,ax = plt.subplots()
bp = ax.boxplot(final_plate.T,labels=dc);

xt = ax.get_xticks()
dc_fit_plot = np.linspace(xt[0],xt[-1],num=100)

ax.plot(dc_fit_plot,spline_fit(dc_fit),color='red',linewidth=2)

# %%

def get_ctx_norm_factor(plate,plate_dc,conc,t,debug=False):

    fluor_data = np.zeros((8,10))

    row_list = ['A','B','C','D','E','F','G','H']

    for row_indx in range(8):
        for col in range(1,11):
            key = row_list[row_indx] + str(col+1)
            ts = np.array(plate.data[key])
            fluor_data[row_indx,col-1] = ts[t]

    fluor_data = fluor_data/np.min(fluor_data)

    mean_fluor = np.mean(fluor_data,axis=1)

    spline_fit = scipy.interpolate.CubicSpline(plate_dc,mean_fluor)

    if debug:
        fig,ax = plt.subplots()
        bp = ax.boxplot(fluor_data.T,labels=dc);

        xt = ax.get_xticks()
        dc_fit_plot = np.linspace(xt[0],xt[-1],num=100)

        ax.plot(dc_fit_plot,spline_fit(dc_fit),color='red',linewidth=2)

    return spline_fit(conc)

# %% animation


fig, ax = plt.subplots()
xdata, ydata = [], []
ln, = ax.plot([], [], 'ro')

def init():
    ax.set_xlim(dc_fit_plot[0],dc_fit_plot[-1])
    ax.set_ylim(1,2)
    ax.set_xlabel('Log drug concentration')
    ax.set_ylabel('Normalized fluorescence')
    return ln,

def update(frame):
    fluor_data = np.zeros((8,10))
    for row_indx in range(8):
        for col in range(1,11):
            key = row_list[row_indx] + str(col+1)
            ts = np.array(p.data[key])
            fluor_data[row_indx,col-1] = ts[frame]

    fluor_data = fluor_data/np.min(fluor_data)

    mean_fluor = np.mean(fluor_data,axis=1)

    spline_fit = scipy.interpolate.CubicSpline(dc_log,mean_fluor)
    ln.set_data(dc_fit_plot,spline_fit(dc_fit))
    ax.set_title(str(frame))
    ax.set_xticks(ax.get_xticks(),['ND',-3,-2,-1,0,1,2,3])
    return ln,

ani = animation.FuncAnimation(fig, update, frames=670,
                    init_func=init, blit=True)

writervideo = animation.FFMpegWriter(fps=60) 
ani.save('spline_through_time.mov', writer=writervideo)
# %%
