import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

def kinetic_eqn(y,t,ka,N,kd):
    r,R = y
    dydt = [-ka*r*N,ka*r*N - kd*R]
    return dydt

N = 10**5
ka = 10**-9
kd = 10**-4

y0 = [10**10,0]

t = np.linspace(0,50000)

fig,ax = plt.subplots()

N_vect = np.logspace(4,8,10)

for N in N_vect:
    sol = odeint(kinetic_eqn,y0,t,args=(N,ka,kd))
    ax.plot(t,sol[:,1],label=N)