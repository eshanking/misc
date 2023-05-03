#%%
from fears.population import Population
from fears.utils import pharm
import matplotlib.pyplot as plt
import numpy as np
#%%

p = Population(k_elim=0.003,k_abs=0.1,death_model=None)

fig,ax = plt.subplots()

k_abs_list = [0.004,0.008,0.016,0.032]

for k_abs in k_abs_list:
    p.k_abs = k_abs
    p.initialize_drug_curve()
    ax.plot(p.drug_curve,label=k_abs,linewidth=2)

ax.legend(frameon=False)
ax.tick_params(axis='both', labelsize=12)
ax.set_ylabel('Drug concenrtation ($\mu$g/mL)',fontsize=14)
ax.set_xlabel('Time (s)',fontsize=14)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

fig.savefig('pharm_ex.png',bbox_inches='tight')

#%%
p = Population(k_elim=0.01,k_abs=1,death_model=None,curve_type='pulsed')

fig,ax = plt.subplots()
ax.plot(p.drug_curve,color='black',linewidth=2)
ax.set_xlim(-10,600)
ax.set_ylabel('Drug conc (a.u.)',fontsize=14)
ax.set_xlabel('Time (hr)',fontsize=14)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.tick_params(axis='both', labelsize=12)

fig.savefig('figures/simulated_serum_conc.png',bbox_inches='tight')

fig,ax = plt.subplots(figsize=(6.4,1))
ax.plot(p.impulses,color='black',linewidth=2)
ax.set_xlim(-10,600)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['left'].set_visible(False)
ax.set_xlabel('Time (hr)',fontsize=14)
ax.tick_params(axis='both', labelsize=12)

fig.savefig('figures/impulses.png',bbox_inches='tight')

impulses = np.zeros(len(p.impulses))
impulses[0] = 1
p.impulses = impulses
dc = p.convolve_pharm(impulses)
fig,ax = plt.subplots(figsize=(6.4,4))

ax.plot(dc,color='black',linewidth=2)
ax.set_ylabel('Drug conc (a.u.)',fontsize=14)
ax.set_xlabel('Time (hr)',fontsize=14)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.tick_params(axis='both', labelsize=12)
fig.savefig('figures/single_dose.png',bbox_inches='tight')
# %%
