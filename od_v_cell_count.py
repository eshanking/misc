import numpy as np
import matplotlib.pyplot as plt

cells = np.array([432,60,(31+17)/2,0])*10**5
cells = cells/50
od = np.array([np.mean([1.28,1.34,1.35]),
               np.mean([0.24,0.29,0.27]),
               np.mean([0.11,0.12,0.12]),
               0.09])

fig,ax = plt.subplots()

ax.scatter(od,cells)

res = np.polyfit(od,cells, deg=1)

x = np.linspace(0,np.max(od),100)
y = res[1] + res[0]*x
ax.plot(x,y,color='r')

# get sci notation
power = np.floor(np.log10(res[0]))
leading = np.log10(res[0]) - power
leading = np.round(10**leading,2)

coeff1 = str(leading) + 'e' + str(int(power))

power = np.floor(np.log10(-res[1]))
leading = np.log10(-res[1]) - power
leading = np.round(10**leading,2)

coeff0 = str(leading) + 'e' + str(int(power))

t = '$\# cells$ = ' + coeff1 + '$*OD_{600}$' + ' - ' + coeff0
ax.annotate(t,(0,3.5*10**7),fontsize=14)


ax.tick_params('both',labelsize=12)

# fig.savefig('OD600_vs_cells.png',bbox_inches='tight',dpi=300)

od_old = od

#%% More OD measurements

proportion = np.array([1,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1,0])
cell_num = (proportion*432*10**5)/50
od = np.array([1.32,1.22,1.13,1.03,0.98,0.866,0.738,0.589,0.44,0.31,0.09])

fig2,ax2 = plt.subplots()

ax2.scatter(od,cell_num)

ydata = np.concatenate((cell_num,cells))
xdata = np.concatenate((od,od_old))

res = np.polyfit(xdata,ydata,deg=2)
x = np.arange(0,max(od)+0.1,step=0.1)
y = res[2] + res[1]*x + res[0]*x**2
ax2.plot(x,y)

ax2.scatter(od_old,cells)

ax2.tick_params('both',labelsize=12)
ax2.set_xlabel('$OD_{600}$',fontsize=14)
ax2.set_ylabel('$\# cells$',fontsize=14)

t = '$Y = 1.48e7X^{2} + 1.35e7X - 4.6e5$'
ax2.annotate(t,(0,3.5*10**7),fontsize=14)

plt.grid(True)
fig2.savefig('OD600_vs_cells.png',bbox_inches='tight',dpi=300)
# %%
