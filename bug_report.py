import numpy as np
import matplotlib.pyplot as plt

#%% Chunk 1
x = np.arange(15)

y = x**2

fig,ax = plt.subplots()

ax.plot(x,y)
#%% Chunk 2
plt.show()

xtl = ax.get_xticklabels()
xtl[1].set_text('New Label')

ax.set_xticklabels(xtl)

plt.show()

fig.savefig('bug2.png')