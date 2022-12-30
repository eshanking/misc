import numpy as np
import matplotlib.pyplot as plt

def steady_state(I,L,m,x):

    u_x = np.zeros(len(x))

    for x_t in x:
        for n in range(m):

            u_x[x_t] += I*np.cos(((1+2*n)*np.pi*x_t)/(2*L))

    return u_x

x = np.arange(100)
m = 1000
L = 100
I = 1

u_x = steady_state(I,L,m,x)

fig,ax = plt.subplots()

ax.plot(u_x)
ax.set_xlim([0,100])