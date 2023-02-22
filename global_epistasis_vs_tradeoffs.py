#%%
from fears.population import Population
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
import seaborn as sns
import pandas as pd
#%%

np.random.seed(2023)

def plot_global_epistasis(p,mut,conc):
    bf,f = get_epistasis_data(p,mut,conc)
    data = pd.DataFrame({'Background fitness':bf,'Fitness effect':f})
    ax = sns.regplot(x='Background fitness',y='Fitness effect',data=data)


def unif_to_norm(rho):
    return (4118/3163)*np.sin(rho/3) + (3183/3149)*np.sin(2*rho/3) - (145/2391)*np.sin(rho)

def get_background(N: int, p: int):
    """Get the genotypes that exclude a mutation in p.

    Adapted from a solution provided by an AI assistant 
    (https://openai.com/blog/better-language-models/)

    Args:
        N (int): total number of genotypes
        p (int): mutation to exclude

    Returns:
        List[str]: background genotypes that exclude p
    """
    # Initialize an empty list to store the binary numbers
    binary_numbers = []

    # Iterate over the range from 0 to 2^N - 1
    for i in range(N):
        # Convert the current number to binary
        binary = bin(i)[2:].zfill(N)

        # Check if the binary number uses the specific power of 2
        uses_power_of_2 = binary[N-p-1] == "1"
        
        # If the binary number does not use the specific power of 2, append it to the list
        if not uses_power_of_2:
            binary_numbers.append(binary)

    # Convert to integers base 10
    numbers = [int(x,2) for x in binary_numbers]

    return numbers

def get_epistasis_data(pop,focal_mut,conc):
    """Calculate the global epistasis for a focal mutation at a specific drug 
    concentration.

    Args:
        pop (_type_): _description_
        focal_genotype (_type_): _description_
        conc (_type_): _description_
    """

    # get the background genotypes that exclude focal_mut
    background_gens = get_background(pop.n_genotype,focal_mut) 

    # now get the genotypes that include focal mut
    focal_mut_gens = set(np.arange(pop.n_genotype))

    focal_mut_gens = focal_mut_gens - set(background_gens)
    focal_mut_gens = list(focal_mut_gens)

    conc = 10**conc

    fl = pop.gen_fit_land(conc)

    background_fitness = fl[background_gens]
    focal_fitness = fl[focal_mut_gens]

    background_gens = np.array(background_gens)
    focal_mut_gens = np.array(focal_mut_gens)

    fitness_effect = focal_fitness - background_fitness

    return background_fitness,fitness_effect

def get_min_max_dge(p):
    # get the max ic50:
    max_ic50 = np.max(p.ic50)

    max_dge = 0
    min_dge = 0

    for a in range(p.n_allele):
        bg,f = get_epistasis_data(p,a,0)
        res0 = np.polyfit(bg,f,deg=1)[0]

        bg,f = get_epistasis_data(p,a,max_ic50/5)
        res1 = np.polyfit(bg,f,deg=1)[0]

        if res0-res1 > max_dge:
            max_dge = res0-res1
        if res0-res1 < min_dge:
            min_dge = res0-res1
        
    return min_dge, max_dge

def get_max_delta_f(p):
    max_ic50 = np.max(p.ic50)

    df_list = []

    for a in range(p.n_allele):
        bg,f0 = get_epistasis_data(p,a,-4)

        bg,f1 = get_epistasis_data(p,a,max_ic50)

        df = np.mean(f0)-np.mean(f1)

        df_list.append(df)

    return np.max(np.abs(df_list))

def get_max_ge_slope(p):
    max_ic50 = np.max(p.ic50)

    slope_list = []

    for a in range(p.n_allele):
        bg0,f0 = get_epistasis_data(p,a,-4)
        res0 = np.corrcoef(bg0,f0)[0][1]

        bg1,f1 = get_epistasis_data(p,a,max_ic50)
        res1 = np.corrcoef(bg1,f1)[0][1]

        slope_list.append(res0-res1)
    
    indx = np.argwhere(np.abs(slope_list) == np.max(np.abs(slope_list)))[0][0]

    return slope_list[indx]

#%% 1 allele

# p = Population(fitness_data='random',n_allele=1)

# delta_global_epistasis = []
# tradeoff_strength = []

# for i in range(1000):
#     ic50 = np.random.uniform(low=-3,high=3,size=2)
#     drugless_rates = np.random.uniform(low=0.5,high=1.5,size=2)

#     ic50 = np.sort(ic50)

#     p.ic50 = ic50
#     p.drugless_rates = drugless_rates

#     # get the global epistasis data at conc = 0:
#     bg,f = get_epistasis_data(p,0,0)
#     res0 = f[0]

#     # get the global epistasis data at conc = 10^3:
#     bg,f = get_epistasis_data(p,0,10**3)
#     res1 = f[0]

#     # compute the difference:
#     dge = res0 - res1

#     ts_t = (drugless_rates[0]-drugless_rates[1])/(ic50[0]-ic50[1])

#     delta_global_epistasis.append(dge)
#     tradeoff_strength.append(ts_t)

#     # fig,ax = p.plot_fitness_curves()
#     # ax.set_title(str(dge) + ' ' + str(ts_t))

# fig,ax = plt.subplots()

# ax.scatter(tradeoff_strength,delta_global_epistasis)
# ax.set_xlim(-1,1)
# ax.set_xlabel('IC$_{50}$/growth rate correlation',fontsize=14)
# ax.set_ylabel('$\Delta$f',fontsize=14)

# ax.spines['top'].set_visible(False)
# ax.spines['right'].set_visible(False)
# %% 2 alleles (4 genotpyes)

p = Population(fitness_data='random',n_allele=2)

df = []
tradeoff_strength = []

for i in range(10000):
    ic50 = np.random.uniform(low=-3,high=3,size=p.n_genotype)
    drugless_rates = np.random.uniform(low=0.5,high=1.5,size=p.n_genotype)

    # ic50 = np.sort(ic50)

    p.ic50 = ic50
    p.drugless_rates = drugless_rates

    df.append(get_max_delta_f(p))

    ts = np.corrcoef(ic50,drugless_rates)
    # ts = stats.spearmanr(ic50,drugless_rates)
    tradeoff_strength.append(ts[0,1])
    # tradeoff_strength.append(ts.correlation)

    # fig,ax = p.plot_fitness_curves()
    # ax.set_title(str(df[-1]))

fig,ax_list = plt.subplots(ncols=2,figsize=(10,4))

ax = ax_list[0]

data_2_allele = pd.DataFrame({'tradeoff_strength':tradeoff_strength,'df':df})
ax = sns.regplot(x='tradeoff_strength',y='df',data=data_2_allele,scatter_kws={'alpha':0.1},
                 line_kws={'color':'darkorange','linewidth':3},ax=ax)

# ax.scatter(tradeoff_strength, df,alpha=0.1)
ax.set_ylabel('Max mean $\Delta GE$',fontsize=14)
ax.set_xlabel('IC$_{50}$ and drugless \ngrowth rate correlation',fontsize=14)
ax.set_title('N$_{allele}$ = 2',fontsize=14)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
# ax.set_ylim(-50,50)

corr_matrix = np.corrcoef(tradeoff_strength, df)
corr = corr_matrix[0,1]
R_sq = corr**2

# print('$R^{2}$ = ' + str(round(R_sq,3)))
# %%

p = Population(fitness_data='random',n_allele=3)

df = []
tradeoff_strength = []
max_tradeoff = []

for i in range(10000):
    ic50 = np.random.uniform(low=-3,high=3,size=p.n_genotype)
    drugless_rates = np.random.uniform(low=0.5,high=1.5,size=p.n_genotype)

    # ic50 = np.sort(ic50)

    p.ic50 = ic50
    p.drugless_rates = drugless_rates

    df.append(get_max_delta_f(p))

    ts = np.corrcoef(ic50,drugless_rates)
    # ts = stats.spearmanr(ic50,drugless_rates)
    tradeoff_strength.append(ts[0,1])
    # tradeoff_strength.append(ts.correlation)

# fig,ax = plt.subplots()

ax = ax_list[1]

data_3_allele = pd.DataFrame({'tradeoff_strength':tradeoff_strength,'df':df})
ax = sns.regplot(x='tradeoff_strength',y='df',data=data_3_allele,scatter_kws={'alpha':0.1},
                 line_kws={'color':'darkorange','linewidth':3},ax=ax)

# ax.scatter(tradeoff_strength, df,alpha=0.1)
ax.set_ylabel('Max mean $\Delta GE$',fontsize=14)
ax.set_xlabel('IC$_{50}$ and drugless \ngrowth rate correlation',fontsize=14)
ax.set_title('N$_{allele}$ = 3',fontsize=14)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

fig.savefig('correlation_plot.pdf',bbox_inches='tight')
# %% 2-allele cartoon

fig2,ax = plt.subplots(figsize=(5,4))

ic50 = [0,-2,2,3]

dr = [1,1.2,0.7,0.6]

p = Population(n_allele=2,fitness_data='random')
p.ic50 = ic50
p.drugless_rates = dr
fig_t,ax = p.plot_fitness_curves(linewidth=4,ax=ax)

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

ax.legend(loc=(1.05,0.1),fontsize=12,frameon=False)

fig2.savefig('ex_dr_curve.pdf',bbox_inches='tight')
# %%

# data = data_2_allele

# q25 = data_2_allele.quantile(q=0.90) # top 25%

# data25 = data.loc[(data['df'] > q25.df),['tradeoff_strength','df']]

# # df = np.array(data25['df'])
# # ts = np.array(data25['tradeoff_strength'])

# # # bin df over tradeoff strength

# # df_bin = []

# # for i in np.arange(-1,1,0.25):
# #     df_bin.append(df[(ts>i)*(ts<i+0.25)])

# fig,ax = plt.subplots()

# ax.scatter(data25['tradeoff_strength'],data25['df'])
# %%

fig,ax = plt.subplots()

ax = sns.jointplot(x='tradeoff_strength',y='df',data=data_3_allele,ax=ax)
# %%
# p = Population(fitness_data='random',n_allele=3)

# dge = []
# tradeoff_corr = []

# for i in range(10000):
#     ic50 = np.random.uniform(low=-3,high=3,size=p.n_genotype)
#     drugless_rates = np.random.uniform(low=0.5,high=1.5,size=p.n_genotype)

#     # ic50 = np.sort(ic50)

#     p.ic50 = ic50
#     p.drugless_rates = drugless_rates

#     dge.append(get_max_ge_slope(p))

#     # ts = stats.pearsonr(ic50,drugless_rates)
#     ts = stats.spearmanr(ic50,drugless_rates)
#     # tradeoff_corr.append(ts.statistic)
#     tradeoff_corr.append(ts.correlation)

#     # fig,ax = p.plot_fitness_curves()
#     # ax.set_title(str(df[-1]))

# fig,ax = plt.subplots(ncols=1,figsize=(5,4))

# data_2_allele = pd.DataFrame({'tradeoff_corr':tradeoff_corr,'dge':dge})
# # ax = sns.regplot(x='tradeoff_corr',y='dge',data=data_2_allele,scatter_kws={'alpha':0.1},
#                 #  line_kws={'color':'darkorange','linewidth':3},ax=ax)
# ax.scatter(tradeoff_corr,dge)
# # ax.scatter(tradeoff_strength, df,alpha=0.1)
# ax.set_ylabel('Max mean $\Delta GE$',fontsize=14)
# ax.set_xlabel('IC$_{50}$ and drugless \ngrowth rate correlation',fontsize=14)
# ax.set_title('N$_{allele}$ = 2',fontsize=14)
# ax.spines['top'].set_visible(False)
# ax.spines['right'].set_visible(False)
# # ax.set_ylim(-50,50)

# corr_matrix = np.corrcoef(tradeoff_strength, df)
# corr = corr_matrix[0,1]
# R_sq = corr**2

# ax.set_ylim(-0.1,0.1)

