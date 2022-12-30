import matplotlib.pyplot as plt
import numpy as np
from scipy import fft

# estimating a noisy signal with scipy and matplotlib

N = 1 # seconds, length of measurement
T = 1/10**4 # sampling freq

np.random.seed(2022)

x = np.linspace(0.0, N, num = int(N/T), endpoint=False)

y = np.sin(100*2*np.pi*x) + np.sin(60*2*np.pi*x) + np.random.normal(loc=0,scale=1,size=len(x))
y_predicted = np.sin(100*2*np.pi*x) + np.sin(60*2*np.pi*x)

fig, ax = plt.subplots(nrows=2,figsize=(6,4))

ax[0].plot(x,y,label='Observed')
ax[0].plot(x,y_predicted,color='red',label='Predicted')
ax[0].set_xlim(0,.1)

yf = fft.fft(y)
xf = np.arange(yf.shape[0])*(1/T)/(len(yf))

ax[1].plot(xf,(T/N)*np.abs(yf))
ax[1].set_xlim(0,200)

ax[0].set_ylabel('Magnitude')
ax[0].set_xlabel('Time (s)')

ax[1].set_ylabel('Magnitude')
ax[1].set_xlabel('Frequency ($s^{-1}$)')

ax[0].legend(frameon=False,ncol=2,loc='upper right')

fig.tight_layout()

fig.savefig('signal_processing.pdf')