from fears.population import Population
from fears.utils import pharm
import matplotlib.pyplot as plt

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