import numpy as np
from scipy.signal import convolve2d
import matplotlib.pyplot as plt

f1 = np.zeros((100,100))

f1[20,50] = 1
f1[80,50] = 1

fig,ax = plt.subplots(nrows=1,ncols=3)

ax[0].imshow(f1)

def exp_decay_2d(x,y):

    r = np.sqrt(x**2+y**2)

    return np.exp(-r)

def make_origin_center(x,y,size):

    x_t = int(x-size/2)
    y_t = int(y-size/2)

    return x_t,y_t

