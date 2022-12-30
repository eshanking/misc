#%%
import os
import pandas as pd
import datetime
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib import colors
# import matplotlib.cm as mplcm
import scipy.optimize as sciopt
import scipy.interpolate as sciinter
import scipy.stats as stats
import re
from fears.population import Population
from fears.utils import plotter
import fears.utils.AutoRate as ar
from pathlib import Path
#%%

# def est_rate_from_diff(y,t):


def rolling_regression(xdata,ydata):

    # compute diff

    r = []
    for i in range(len(ydata)-1):
        dy = ydata[i+1] - ydata[i]
        dx = xdata[i+1] - xdata[i]
        r.append(dy/dx)

    cc_est = np.mean(ydata[-2:])

    if cc_est < 0.3:
        return 10**-6
    else:
        return np.max(r)

def logistic_pharm_curve(x,IC50,g_drugless,hill_coeff):
    """Defines the logistic dose-response curve. Use if the input is a vector of drug concentration curves

    Args:
        x (numpy array): drug concentration vector
        IC50 (float)): IC50
        g_drugless (float): drugless growth rate
        hill_coeff (float): Hill coefficient

    Returns:
        numpy array: array of growth rates
    """
    g = []

    for x_t in x:
        if x_t == 0:
            g.append(g_drugless)
        else:
            g.append(g_drugless/(1+np.exp((IC50-np.log10(x_t))/hill_coeff)))

    return g

def logistic_growth_curve(t,r,p0,k):
    """Logistic growth equation

    Args:
        t (float): time
        r (float): growth rate
        p0 (float): starting population size
        k (float): carrying capacity

    Returns:
        float: population size at time t
    """
    p = k/(1+((k-p0)/p0)*np.exp(-r*t))

    return p

def logistic_growth_with_lag(t,r,p0,k,l):

    p = p0 + k/(1+ np.exp((4*r*(l-t)/k) + 2))

    return p

def est_ic50(gr,dc):
    indx = np.argwhere(gr)[0][0]
    ic50 = dc[indx]
    return np.log10(ic50)

def fit_hill_curve(xdata,ydata,hc=None,debug=False,interpolate=False):
    """Fits dose-response curve to growth rate data

    Args:
        xdata (list or numpy array): drug concentration curve from plate experiment
        ydata (list or numpy array): growth rate versus drug concetration for a given replicate

    Returns:
        list: List of optimized paramters: IC50, drugless growth rate, and Hill coefficient
    """
    if interpolate:
        # interpolate data
        xd_t = xdata
        yd_t = ydata
        f = sciinter.interp1d(xdata,ydata)

        if min(xdata) == 0:
            xmin = np.log10(xdata[1]) # if xdata starts at zero, set the new xmin to be the log of the next smallest value
        else:
            xmin = np.log10(min(xdata))
        xmax = np.log10(max(xdata))

        xdata = np.logspace(xmin,xmax)
        if not xdata[0] == 0:
            xdata = np.insert(xdata,0,0) # add zero back to xdata if removed before (because of log(0) error)

        ydata = f(xdata) # interpolate new ydata points
    else:
        xd_t = xdata
        yd_t = ydata

    if hc is None: # want to estimate hill coefficient as well

        ic50_est = est_ic50(yd_t,xd_t)

        p0 = [ic50_est,ydata[-1],-0.1]

        # if ydata[0] == 0:
        #     g_drugless_bound = [0,1]
        # else:
        #     # want the estimated drugless growth rate to be very close to the value given in ydata
        #     g_drugless_bound = [ydata[0]-0.0001*ydata[0],ydata[0]+0.0001*ydata[0]]
        
        gd_est = ydata[-1] #drugless growth rate estimate
        gd_lb = gd_est - 0.5*gd_est
        gd_ub = gd_est + 0.5*gd_est

        bounds = ([ic50_est-0.5,gd_lb,-0.1],[ic50_est+0.5,gd_ub,-0.001]) # these aren't magic numbers these are just starting parameters that happen to work

        popt, pcov = sciopt.curve_fit(logistic_pharm_curve,
                                            xdata,ydata,p0=p0,bounds=bounds)
    
    # else: # we already know the hill coefficient, estimate everything else
    #     p0 = [0,ydata[0]]

    #     # print(p0)

    #     if ydata[0] == 0:
    #         g_drugless_bound = [0,1]
    #     else:
    #         # want the estimated drugless growth rate to be very close to the value given in ydata
    #         g_drugless_bound = [ydata[0]-0.0001*ydata[0],ydata[0]+0.0001*ydata[0]]

    #     fitfun = partial(self.logistic_pharm_curve_vectorized,hill_coeff=hc)

    #     bounds = ([-5,g_drugless_bound[0]],[4,g_drugless_bound[1]])
    #     popt, pcov = scipy.optimize.curve_fit(fitfun,
    #                                         xdata,ydata,p0=p0,bounds=bounds)            

    d = {'ic50':popt[0],
        'g_drugless':popt[1],
        'hc':popt[2]}


    if debug:
        # est = [logistic_pharm_curve(x,popt[0],popt[1],popt[2]) for x in xdata]
        est = logistic_pharm_curve(xdata,popt[0],popt[1],popt[2])
        fig,ax = plt.subplots()

        xd_plot = [np.log10(x) for x in xd_t if x > 0]

        xd_plot = xd_plot + [-3]

        ax.scatter(xd_plot,yd_t,marker='*')
        ax.scatter(popt[0],popt[1]/2,color='r',marker='o')
        # ax.plot(xdata,ydata)
        ax.plot(xd_plot,est,color='black')
        # ax.set_xscale('log')

        title_t = 'IC50 = ' + str(round(popt[0],2)) + ' IC50_est = ' + str(ic50_est) + '\n hc = ' + str(round(popt[2],3)) + \
                ' g = ' + str(round(popt[1],2)) + ' y0 = ' + str(round(ydata[-1],2))

        ax.set_title(title_t)

    return d

def atoi(text):
    return int(text) if text.isdigit() else text

def natural_keys(text):
    '''
    alist.sort(key=natural_keys) sorts in human order
    http://nedbatchelder.com/blog/200712/human_sorting.html
    (See Toothy's implementation in the comments)
    '''
    return [ atoi(c) for c in re.split(r'(\d+)', text) ]

def linear_growth_curve(t,r,p0):
    return p0+r*t

def est_logistic_params(growth_curve,t,debug=False,sigma=None,mode='logistic',
                        normalize=False):
    """Estimates growth rate from OD growth curve

    Args:
        growth_curve (list or numpy array): vector of OD data
        t (list or numpy array, optional): Time vector. If None, algorithm assumes each time step is 1 s. Defaults to None.

    Returns:
        dict: Dict of logistic growth curve paramters
    """
    # interpolate negative time

    # t_ext = t[0]

    # t = [tt+t[1] for tt in t]

    # t_ext = [t_ext] + t

    # t = t_ext

    # growth_curve = [growth_curve[0]] + growth_curve

    # normalize
    if normalize:
        norm_factor = max(growth_curve)
        growth_curve = [g/norm_factor for g in growth_curve]
    else:
        norm_factor = 1

    # estimate cc

    cc_est = np.mean(growth_curve[-2:])

    gr_est = rolling_regression(t,growth_curve)

    if cc_est <= 0:
        
        d = {'gr':0,
        'OD_0':0,
        'OD_max':0,
        'lambda':0} 
        
        return d,None
    
    bounds = ([gr_est-0.2*gr_est,growth_curve[0]-0.1,cc_est-cc_est*0.05,-1000],

              [gr_est+0.2*gr_est,growth_curve[0]+0.1,cc_est+cc_est*0.05,max(t)])

    # print(bounds)

    # p0 = [10**-3,growth_curve[0],cc_est] # starting parameters

    # popt, pcov = sciopt.curve_fit(logistic_growth_curve,
    #                                     t,growth_curve,p0=p0,sigma=sigma,
    #                                     bounds=bounds)
    p0 = [gr_est,growth_curve[0],cc_est,1000]

    popt,pcov = sciopt.curve_fit(logistic_growth_with_lag,t,growth_curve,p0=p0,
                                 bounds=bounds)
    
    rate_indx = 0 # index of the optimized data points referring to the growth rate
    p0_indx = 1 # index of the optimized data points referring to initial population size
    carrying_cap_indx = 2 # index of the optmized data points referring to carrying capacity
    lamba_indx = 3

    r = popt[rate_indx]
    p0 = popt[p0_indx]*norm_factor
    cc = popt[carrying_cap_indx]*norm_factor
    l = popt[lamba_indx]

    min_carrying_cap = 0.4

    if r < 0: # if the growth rate is negative
        r = 0
    if cc < p0: # if the carrying capacity is less than the initial population size
        r = 0
    if cc < min_carrying_cap: # if the carrying cap is less the minimum threshold
        r = 0
    if norm_factor < 0.4:
        r = 0
            
    d = {'gr':r,
        'OD_0':p0,
        'OD_max':cc,
        'lambda':l}   

    if debug:
        if r > 0:
            fig,ax = plt.subplots()
            t_plot = np.array(t)/3600

            diff = []

            for i in range(len(growth_curve)-1):
                diff.append(growth_curve[i+1]-growth_curve[i])

            ax.scatter(t_plot,growth_curve)
            ax.scatter(t_plot[1:],diff)

            est = [logistic_growth_with_lag(tt,popt[0],popt[1],popt[2],popt[3]) for tt in t]
            
            ax.plot(t_plot,est,color='red')

            p0 = round(popt[1]*10**5)/10**5
            k = round(popt[2]*10**5)/10**5
            r_t = round(popt[0]*3600,2)
            title = 'rate = ' + str(round(r_t,2)) + ' K = ' + str(k) + ' p0 = ' + str(round(p0,2))

            ax.set_title(title)
            ax.set_xlabel('Time (hr)')
            if normalize:
                ax.set_ylabel('Normalized OD')
            else:
                ax.set_ylabel('OD')
            ax.set_ylim(0,1)   

    return d,pcov

def get_background(data_dict,bg_keys):
    avg = 0
    for key in bg_keys:
        avg+=data_dict[key]
    avg = avg/len(bg_keys)
    return avg

def od_data_to_dict(df):
    """Takes an OD data file and returns a dict with each data well as a key

    OD data files refers to the type of data that is a single OD reading in time (i.e
    not timeseries data).

    Args:
        df (pandas dataframe): Parsed OD data file

    Returns:
        dict: Dict of key-value pairs where each key is a well location and each value
        is the OD reading.
    """
    rownum = 0
    d = {}

    for k0 in df['Rows']:
        for k1 in df.columns[1:]:
            if k0.isnumeric():
                key = k1+k0
            else:
                key = k0+k1
            d[key] = df[k1][rownum]
        rownum+=1

    return d                    

def parse_od_data_file(df):
    """Loads the raw OD data files and strips metadata

    OD data files refers to the type of data that is a single OD reading in time (i.e
    not timeseries data).

    Args:
        data_path (str): path to data file

    Returns:
        pandas dataframe: Formatted dataframe with metadata stripped
    """

    # get the first column as an array
    col_0 = df.columns[0]
    col_0_array = np.array(df[col_0])

    data_start_indx = np.argwhere(col_0_array=='<>')
    data_start_indx = data_start_indx[0][0]

    # the data end index should be the first NaN after the data start index
    col_0_array_bool = [pd.isna(x) for x in col_0_array]
    data_end_indx = np.argwhere(col_0_array_bool[data_start_indx:])
    data_end_indx = data_end_indx[0][0] + data_start_indx - 1

    df_filt = df.loc[data_start_indx:data_end_indx,:]

    # fix the columns
    i = 0
    columns = list(df_filt.iloc[0])
    columns_t = []
    for c in columns:
        if type(c) is not str:
            columns_t.append(str(int(c)))
        else:
            columns_t.append(c)
        i+=1
    
    columns_t[0] = 'Rows'

    df_filt.columns = columns_t

    df_filt = df_filt.drop(df_filt.index[0])
    df_filt = df_filt.reset_index(drop=True)

    return df_filt

def get_plate_paths(folder_path):
    """Gets plate data paths
    Returns:
        list: list of plate data paths
    """
    plate_files = os.listdir(path=folder_path)

    #Need to make sure we are only attempting to load .csv or .xlsx data
    plate_files = [i for i in plate_files]

    plate_files.sort()

    plate_data_paths = []

    for pf in plate_files:
        if pf != '.DS_Store':
            plate_path = folder_path + os.sep + pf
            plate_data_paths.append(plate_path)

    plate_data_paths.sort(key=natural_keys)
    return plate_data_paths

def get_data_file_paths(plate_path):
    files = os.listdir(path=plate_path)

    #Need to make sure we are only attempting to load .csv or .xlsx data
    files = [i for i in files if ('.csv' in i) or ('.xlsx' in i)]

    files.sort()

    file_data_paths = []

    for pf in files:
        if pf != '.DS_Store':
            file_path = plate_path + os.sep + pf
            file_data_paths.append(file_path)

    return file_data_paths

def get_start_time(df,col=4):

    # first start time is shaking, so we take the second (start of scan)
    f = df[df == 'Start Time'].stack().index.tolist()[1]

    row = f[0]
    date_time = df.iloc[row,col]

    yr = int(date_time[0:4])
    mon = int(date_time[5:7])
    day = int(date_time[8:10])

    hr = int(date_time[11:13])
    min = int(date_time[14:16])
    sec = int(date_time[17:19])

    dt = datetime.datetime(yr,mon,day,hour=hr,minute=min,second=sec)


    return dt

#%%
# def make_fig():
    # bg_keys = ['A12','B12','C12','D12','E12','F12','G12','H12']
drug_conc = [10000,2000,400,80,16,3.2,0.64,0.128,0.0256,0.00512,0,'control']
folder_path = '/Users/kinge2/repos/seascapes_figures/data/08312022'

plate_paths = get_plate_paths(folder_path)

# fig,ax_list = plt.subplots(ncols=4,nrows=4,figsize=(15,12))

count = 0

gr_lib = {}

rate_est_lib = {}

all_timeseries = {}
all_log_params = {}

# for pp in plate_paths:
for pp in [plate_paths[0]]:

    row = int(np.floor(count/4))
    col = int(np.mod(count,4))

    # ax = ax_list[row,col]

    # fig,ax = plt.subplots()

    data_paths0 = get_data_file_paths(pp)

    timeseries_dict = {}
    timeseries_dict['Time'] = []
    logistic_params_dict = {}

    for p in data_paths0:

        df = pd.read_excel(p)
        t = get_start_time(df)
        df = parse_od_data_file(df)

        data_dict = od_data_to_dict(df)
        data_dict['Time'] = t

        # bg = get_background(data_dict,bg_keys)
        for key in data_dict:
            if key != 'Time':
                if key in timeseries_dict.keys():
                    od = data_dict[key]
                    timeseries_dict[key].append(data_dict[key])
                else:
                    od = data_dict[key]
                    timeseries_dict[key] = [od]
        timeseries_dict['Time'].append(t)

    # sort out time
    t_vect = timeseries_dict['Time']
    t0 = t_vect[0]
    t_vect = [(t-t0).total_seconds() for t in t_vect]

    # get summary data

    replicates = ['A','B','C','D','E','F','G','H']
    conditions = np.linspace(1,12,num=12)
    conditions = [str(int(c)) for c in conditions]

    data_avg = {}
    data_std = {}
    
    gr_avg = []
    gr_std = []

    rate_est_dict = {}

    for c in conditions:
        rate_est = []
        for r in replicates:

            if not (count == 6 and (r == 'A' or r == 'H')): 

                key = r+c

                # control_key = r+'12'

                # control_vect = np.array(timeseries_dict[control_key])

                od_vect = np.array(timeseries_dict[key])

                # od_vect = od_vect-control_vect

                od_vect= np.array(od_vect)

                od_diff = []
                r_est = []
                
                cc_est = np.mean(od_vect[-2:])

                if cc_est < 0.4:
                    r_est.append(0)
                else:

                    for i in range(len(od_vect)-1):
                        diff = (od_vect[i+1]-od_vect[i])/(t_vect[i+1]-t_vect[i])
                        od_diff.append(diff)

                        if od_vect[i]/cc_est < 1 and diff>0:

                            r_est.append(diff/(od_vect[i]*(1-od_vect[i]/cc_est)))


                # plt.plot(od_diff)

                timeseries_dict[key] = od_vect

                d,pcov = est_logistic_params(od_vect,t_vect,debug=False,mode='logistic',
                        normalize=False)

                logistic_params_dict[key] = d

                if len(r_est) == 0:
                    r_est.append(0)
                # else:
                plt.plot(od_vect)

                rate_est = rate_est + r_est
            else:
                key = r+c
                del timeseries_dict[key] # remove these two data points as outliers
            # r = rolling_regression(t_vect,od_vect)
            # rate_est.append(r)

        rate_est_dict[c] = rate_est

        # num_zeros = np.arghwere(rate_est)
        print(rate_est)
        gr_avg.append(np.mean(rate_est))
        gr_std.append(np.std(rate_est))


    grl_t = {'avg':gr_avg,
            'err':gr_std}
    
    gr_lib[str(count)] = grl_t
    rate_est_lib[str(count)] = rate_est_dict

    all_timeseries[str(count)] = timeseries_dict
    all_log_params[str(count)] = logistic_params_dict

    count+=1

#%% compute seascape

seascape_lib = {}
for key in gr_lib:
    
    gr_t = gr_lib[key]['avg'][0:-1]
    gr_t = [gr*3600 for gr in gr_t]

    dc = [d for d in drug_conc if type(d) != str]

    d = fit_hill_curve(dc,gr_t,debug=False,interpolate=False)
    seascape_lib[key] = d

#%% plot all hill curve fits

cc = plotter.gen_color_cycler()

fig1,ax_list = plt.subplots(ncols=4,nrows=4,figsize=(15,12))
# ax2.set_prop_cycle(cc)

dc_fit = np.logspace(-3.5,4)
dc_plot = [np.log10(dc) for dc in dc_fit]

dc_log = drug_conc
dc_log = [dc for dc in dc_log if type(dc) != str]
dc_log = [np.log10(dc) for dc in dc_log if dc > 0]
dc_log = dc_log + [-3]

count = 0
for key in gr_lib:
    row = int(np.floor(count/4))
    col = int(np.mod(count,4))

    ax = ax_list[row,col]
    gr_t = [g*3600 for g in gr_lib[key]['avg'][0:-1]]

    gr_err = [e*3600 for e in gr_lib[key]['err'][0:-1]]

    ic50 = seascape_lib[key]['ic50']
    g_drugless = seascape_lib[key]['g_drugless']
    hc = seascape_lib[key]['hc']

    g_est = logistic_pharm_curve(dc_fit,ic50,g_drugless,hc)

    ax.errorbar(dc_log,gr_t,yerr=gr_err,fmt='o')

    ax.plot(dc_plot,g_est)
    ax.scatter(ic50,g_drugless/2,color='r',marker='*')
    
    ax.set_ylim(0,1.1)

    # plt.show()
    # xtl = ax.get_xticklabels()
    # xt = ax.get_xticks()
    # # xtl_new = xtl
    # xtl[1].set_text('nd')
    # ax.set_xticks(xt)
    # ax.set_xticklabels(xtl)

    ax.set_xlim(-3.5,4)
    ax.set_title(key)

    if col == 0:
        ax.set_ylabel('Growth rate ($hr^{-1}$)',fontsize=12)
    if row == 3:
        ax.set_xlabel('Log drug concentration (ug/mL)',fontsize=12)
    ax.tick_params(axis='both', labelsize=12)
    count +=1

fig1.tight_layout()

#%% plot seascape

seascape_lib = {}
for key in gr_lib:
    
    gr_t = gr_lib[key]['avg'][0:-1]
    gr_t = [gr*3600 for gr in gr_t]

    dc = [d for d in drug_conc if type(d) != str]

    d = fit_hill_curve(dc,gr_t,debug=False,interpolate=False)
    seascape_lib[key] = d

fig3,ax3 = plt.subplots(figsize=(10,6))
ax3.set_prop_cycle(cc)
dc_fit = np.logspace(-3.5,4,num=100)

for key in seascape_lib:
    
    ic50 = seascape_lib[key]['ic50']
    g_drugless = seascape_lib[key]['g_drugless']
    hc = seascape_lib[key]['hc']

    g_est = logistic_pharm_curve(dc_fit,ic50,g_drugless,hc)

    ax3.plot(dc_fit,g_est,linewidth=2,label=int(key))

ax3.set_xscale('log')

ax3.set_ylabel('Growth rate ($hr^{-1}$)',fontsize=15)
ax3.set_xlabel('Drug Concentration (ug/mL)',fontsize=15)
ax3.spines['top'].set_visible(False)
ax3.spines['right'].set_visible(False)

ax3.tick_params(axis='both', labelsize=12)
ax3.legend(loc=(1,0),frameon=False)
# ax3.tick_params(axis='both', which='minor', labelsize=8)
#%%
fig4,ax4 = plt.subplots(figsize=(5,5))
ax4.set_prop_cycle(cc)
ic50_list = []
gr_list = []
p = Population()

cmap = cm.get_cmap('magma',6)
cmap = cmap.colors

mut_v_gr = [[],[],[],[],[]]

mut_list = []

for key in seascape_lib.keys():

    # if key != '3':
    ic50 = seascape_lib[key]['ic50']
    g_drugless = seascape_lib[key]['g_drugless']

    ic50_list.append(ic50)
    gr_list.append(g_drugless)

    # key_bin = int(key)

    key_bin = p.int_to_binary(int(key))

    num = 0
    for s in key_bin:
        num+=int(s)

    mut_list.append(num)

    mut_v_gr[num].append(g_drugless)

    ax4.scatter(ic50,g_drugless,marker='o',s=400,facecolor=cmap[num+1],
                edgecolors='w',label=int(num))
    # ax4.annotate(key,(ic50-0.15,g_drugless-0.001),fontsize=12)
    ax4.annotate(key,(ic50,g_drugless),fontsize=12,ha='center',va='center')

# ax4.set_ylim(0.06,0.115)
# ax4.set_xlim(-3,4)
ax4.set_ylabel('Drug-free growth rate ($hr^{-1}$)',fontsize=14)
ax4.set_xlabel('Log IC50 (ug/mL)',fontsize=14)
ax4.tick_params(axis='both', labelsize=13)

handles, labels = ax4.get_legend_handles_labels()

unique_labels = sorted(set(labels))
labels = np.array(labels)
unique_handles = []

for lab in unique_labels:
    indx = np.argwhere(labels==lab)
    indx = indx[0][0]
    unique_handles.append(handles[indx])

ax4.legend(unique_handles,unique_labels,loc = (1,0),frameon=False,
             fontsize=12,title='$n_{mut}$')

gr_list = np.array(gr_list)
ic50_list = np.array(ic50_list)

# gr_list = gr_list/np.max(gr_list)
# ic50_list = ic50_list/np.max(ic50_list)

# tradeoff_stats = stats.pearsonr(ic50_list,gr_list)

# gr_list = np.delete(gr_list,3)
# ic50_list = np.delete(ic50_list,3)

# gr_list = gr_list/np.max(gr_list)
# ic50_list = ic50_list/np.max(ic50_list)

# # gr_list = np.array(gr_list)/np.max(gr_list)
# tradeoff_stats_no3 = stats.pearsonr(ic50_list,gr_list)

#%% Mutations vs growth rate
fig8, ax_list = plt.subplots(ncols=2,figsize=(8,4))

ax = ax_list[0]

wt_gr = seascape_lib['0']['g_drugless']
gr_list_norm= gr_list/wt_gr

ax.scatter(mut_list,gr_list_norm)

ax.set_xticks([0,1,2,3,4])

res = stats.linregress(mut_list,gr_list_norm)

x = np.arange(5)
y = res.slope*x + res.intercept

ax.plot(x,y,color='orange',linewidth=2)

ax.set_xlabel('# mutations',fontsize=14)
ax.set_ylabel('Normalized growth rate',fontsize=14)
ax.tick_params(axis='both', labelsize=13)

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

ann = '$r^2$ = ' + str(round(res.rvalue**2,2)) + '\n$p$ = ' + str(round(res.pvalue,3))
ax.annotate(ann,(2.5,1.2),fontsize=12)

# Mutations vs IC50
ax = ax_list[1]

wt_ic50 = seascape_lib['0']['ic50']

ic50_list_norm = [10**x for x in ic50_list]

ic50_list_norm= ic50_list_norm/(10**wt_ic50)

ic50_list_norm = [np.log10(x) for x in ic50_list_norm]

ax.scatter(mut_list,ic50_list_norm)

res = stats.linregress(mut_list,ic50_list_norm)

x = np.arange(5)
y = res.slope*x + res.intercept

ax.plot(x,y,color='orange',linewidth=2)

ax.set_xticks([0,1,2,3,4])

ax.set_xlabel('# mutations',fontsize=14)
ax.set_ylabel('Normalized IC50',fontsize=14)
ax.tick_params(axis='both', labelsize=13)

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

ann = '$r^2$ = ' + str(round(res.rvalue**2,2)) + '\n$p$ = ' + str(round(res.pvalue,3))
ax.annotate(ann,(2.5,0.2),fontsize=12)

fig8.tight_layout()

#%% Weinreich MIC comparison

mic_list = [0.088,1.4,0.063,32,0.13,3.6*10**2,0.18,3.6*10**2,0.088,23,1.4,
            3.6*10**2,1.4,2.1*10**3,0.8,2.9*10**3]

mic_list = [np.log10(m) for m in mic_list]

fig5,ax5 = plt.subplots(figsize=(5,5))
count = 0

ic50_list = []

for key in seascape_lib.keys():
    ic50 = seascape_lib[key]['ic50']
    ic50_list.append(ic50)
    mic = mic_list[count]
    ax5.scatter(mic,ic50,marker='o',s=400,facecolor='w',edgecolors='b')
    ax5.annotate(key,(mic,ic50),fontsize=12,ha='center',va='center')
    print(key)
    count+=1

ax5.set_xlabel('Weinreich MIC (log(ug/mL))',fontsize=14)
ax5.set_ylabel('IC50 (log(ug/mL))',fontsize=14)

yl = ax5.get_ylim()
xl = ax5.get_xlim()

lim = (min(yl[0],xl[0]),max(yl[1],xl[1]))
ax5.set_xlim(lim)
ax5.set_ylim(lim)

ax5.plot(lim,lim,'--',color='black')

ax5.tick_params(axis='both', labelsize=13)

weinreich_stats = stats.pearsonr(mic_list,ic50_list)


#%% Spot check plate 6

# fig7,ax_list = plt.subplots(3,4,figsize=(10,8))

# count = 0

# for c in conditions:
#     row = int(np.floor(count/4))
#     col = int(np.mod(count,4))
#     ax = ax_list[row,col]
#     for r in replicates:
#         key = r + c
#         od = timeseries_dict[key]
#         ax.plot(t_vect,od,label=r)
    
#     ax.set_ylim(0,1)

#     if count == 11:
#         ax.set_title('Neg control')
#     else:
#         dc = drug_conc[count]

#         ax.set_title('Conc = ' + str(dc))

#     count+=1

# # ax_list[0][3].legend()
# fig7.tight_layout()

#%%
# fig1.savefig('figures/all_hill_fits.pdf',bbox_inches='tight')
# fig3.savefig('figures/new_ecoli_seascape.pdf',bbox_inches='tight')
# fig4.savefig('figures/gr_v_ic50.pdf',bbox_inches='tight')
# fig5.savefig('figures/weinreich_MIC_comparison.pdf',bbox_inches='tight')
# fig8.savefig('figures/mutation_tradeoffs.pdf',bbox_inches='tight')

# df = pd.DataFrame(seascape_lib)
# df.to_excel('results/seascape_library.xlsx')

#%%
#%%

cmap = cm.get_cmap('tab10',12)
cmap = cmap.colors

genotype = 7

plate_t = all_timeseries[str(genotype)]
log_params_t = all_log_params[str(genotype)]

replicates = ['A','B','C','D','E','F','G','H']

# for col in replicates:
col = 'A'
fig6,ax = plt.subplots()
for row in range(12):
    od = np.zeros(len(plate_t['Time']))
    gr = 0
    OD_0 = 0
    OD_max = 0
    l = 0
# for col in replicates:
    # col = 'B'
    key = col + str(row+1)
    od += np.array(plate_t[key])

    d = log_params_t[key]
    gr += d['gr']
    OD_0 += d['OD_0']
    OD_max += d['OD_max']
    l += d['lambda']

    if gr == 0:
        # OD_max = 0.1
        OD_0 = 0.1

    # # od = od/8  
    # if max(od) > 0.4: 
    #     od = od/max(od) 
    ax.plot(t_vect,od,'*',color = cmap[row],linewidth=1.5)

    # gr = gr/8
    # OD_0 = OD_0/8
    # OD_max = OD_max/8
    # l = l/8

    # OD_max = 1
    t = np.arange(min(t_vect),max(t_vect))
    od_est = [logistic_growth_with_lag(tt,gr,OD_0,OD_max,l) for tt in t]
    if type(drug_conc[row]) != str and drug_conc[row] != 0:
        label = str(round(np.log10(drug_conc[row]),2))
    else:
        label = str(drug_conc[row])
    ax.plot(t,od_est,color=cmap[row],label=label)
    # ax.set_title(col)

ax.set_ylabel('$OD_{600}$',fontsize=12)
ax.set_xlabel('Time (s)',fontsize=12)
ax.spines["right"].set_visible(False)
ax.spines["top"].set_visible(False)
ax.tick_params(axis='both', labelsize=12)
ax.set_xlim(0,60000)

ax.legend(loc=(1.05,0.1),title='Log Drug conc. \n(ug/mL)',frameon=False)

# fig6.tight_layout()
# fig6.subplots_adjust(right=0.7)
fig6.savefig('figures/example_growth_curves.png',bbox_inches='tight')

#%%
fig7,ax = plt.subplots()
    
gr = gr_lib[str(genotype)]['avg'][0:-1]
gr_err = gr_lib[str(genotype)]['err'][0:-1]

gr_err= [3600*err for err in gr_err]
gr = [3600*g for g in gr]

dc_log = drug_conc
dc_log = [dc for dc in dc_log if type(dc) != str]
dc_log = [np.log10(dc) for dc in dc_log if dc > 0]
dc_log = dc_log + [-3]

ax.errorbar(dc_log,gr,yerr=gr_err,fmt='o')
# ax.scatter(dc_log,gr)
# ax.set_xscale('log')

dc_t = np.logspace(-4,4,num=1000)

key = str(genotype)
ic50 = seascape_lib[key]['ic50']
g_drugless = seascape_lib[key]['g_drugless']
hc = seascape_lib[key]['hc']

g_est = logistic_pharm_curve(dc_t,ic50,g_drugless,hc)

dc_log = [np.log10(dc) for dc in dc_t]

ax.plot(dc_log,g_est)
ax.set_xlim(-3.5,4.5)

ax.set_ylabel('Growth rate ($hr^{-1}$',fontsize=12)
ax.set_xlabel('Drug conc. (log(ug/mL))',fontsize=12)
ax.spines["right"].set_visible(False)
ax.spines["top"].set_visible(False)
ax.tick_params(axis='both', labelsize=12)

fig7.savefig('figures/example_dose_response.png',bbox_inches='tight')

    # return all_timeseries, all_log_params
    # fig7.savefig('spot_check_plate_9.pdf')