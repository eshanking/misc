import numpy as np
import matplotlib.pyplot as plt

K = 4
OD_n = [3.5,3,2.5,2]
OD_0 = 0.95*K

r_0 = 1

def calc_L0(K,p0):
    return (K-p0)/p0

def ratio(K,OD_n,OD_0,p0):

    L0 = calc_L0(K,p0)
    num = np.log10((K-OD_n)/(OD_n*L0))
    denom = np.log10((K-OD_0)/(OD_0*L0))

    return num/denom

p0 = np.linspace(0.001,0.1,num=100)


fig,ax = plt.subplots()
for OD in OD_n:

    rate_est = [ratio(K,OD,OD_0,p)*r_0 for p in p0]

    ax.plot(p0,rate_est,label=str(OD))

    ax.set_xlabel('$p_{0}$',fontsize=15)
    ax.set_ylabel('Estimated rate',fontsize=15)

ax.legend(frameon=False)