from fears.population import Population
import numpy as np
from typing import List
import matplotlib.pyplot as plt

def get_background(N: int, p: int) -> List[int]:
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

    # not get the genotypes that include focal mut
    focal_mut_gens = set(np.arange(pop.n_genotype))

    focal_mut_gens = focal_mut_gens - set(background_gens)
    focal_mut_gens = list(focal_mut_gens)

    fl = pop.gen_fit_land(conc)

    background_fitness = fl[background_gens]
    focal_fitness = fl[focal_mut_gens]

    # genotype 3 always has fitness zero-- want to remove it

    background_gens = np.array(background_gens)
    focal_mut_gens = np.array(focal_mut_gens)

    if 3 in background_gens:
        # figure out which index correponds to gen 3
        indx = np.argwhere(background_gens==3)[0][0]
        background_fitness = np.delete(background_fitness,indx)
        focal_fitness = np.delete(focal_fitness,indx)
    elif 3 in focal_mut_gens:
        # figure out which index correponds to gen 3
        indx = np.argwhere(focal_mut_gens==3)[0][0]
        background_fitness = np.delete(background_fitness,indx)
        focal_fitness = np.delete(focal_fitness,indx)

    fitness_effect = focal_fitness - background_fitness

    return background_fitness,fitness_effect

p = Population(fitness_data='two-point')

p.ic50 = p.ic50+6

# compute_global_epistasis(p,3,10)

dc = np.logspace(-2,3,num=6)
dc = np.concatenate(([0],dc))
corr_list = []

fig,ax_list = plt.subplots(nrows=4,ncols=7,layout='constrained',figsize=(8,4))

focal_gen_list = [2,0,3,1]

for r in range(4):
    indx = 0
    ymax = 0
    ymin = 1
    for conc in dc:
        focal_gen = focal_gen_list[r]
        bg,f = get_epistasis_data(p,focal_gen,conc)            

        # bg = bg - np.min(bg)
        # bg = bg/np.max(bg)

        # f = f - np.min(f)
        # f = f/np.max(f)
        
        ax = ax_list[r,indx]
        ax.scatter(bg,f,color='black',s=10)
        # ax.ticklabel_format(axis='y',style='sci',scilimits=(0,0))

        # ax.set_xlim(-0.1,1.1)
        # ax.set_ylim(-0.1,1.1)
        res = np.polyfit(bg,f,deg=1)
        xl = ax.get_xlim()
        xl = np.array(xl)
        yfit = xl*res[0] + res[1]

        ax.plot(xl,yfit)

        # ax.set_ylim(-0.2,0.7)

        ax.set_yticks([])
        
        yl = ax.get_ylim()

        ymax_t = yl[1]
        ymin_t = yl[0]

        if ymax_t > ymax:
            ymax = ymax_t
        if ymin_t < ymin:
            ymin = ymin_t

        indx+=1
    
    for ax in ax_list[r,:]:
        ax.set_ylim(ymin-0.1,ymax+0.1)



# for ax in ax_list[:,0]:
#     ax.set_yticks([0,0.5])