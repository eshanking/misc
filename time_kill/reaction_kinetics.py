import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

def kinetic_eqn(y,t,ka,N,kd):
    r,R = y
    dydt = [-ka*r*N,ka*r*N - kd*R]
    return dydt

ka = 10**-8
kd = 0.1*10**-4

y0 = [1,0]

t = np.linspace(0,80000)

fig,ax = plt.subplots()

N_vect = np.logspace(3,8,7)

for N in N_vect:
    sol = odeint(kinetic_eqn,y0,t,args=(N,ka,kd))
    ax.plot(t,sol[:,1],label='%.1E' % N)

ax.legend(frameon=False,title='Population size',loc='upper right')
ax.set_ylabel('Normalized intensity')
ax.set_xlabel('Time (s)')