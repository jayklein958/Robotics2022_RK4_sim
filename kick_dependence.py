#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 25 16:33:35 2022

@author: jonathanklein
"""

#!/usr/bin/env python
# coding: utf-8

# In[1]:
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import numpy as np
import scipy.signal
from scipy import fftpack
import math
import numpy
import pygame
import time
import sys
from mpl_toolkits import mplot3d
import numpy as np
import matplotlib.pyplot as plt
from DP2_RK4_single_kick import DoublePendulumLagrangian_kick

def plot(x_data,y_data1):
    plt.plot(x_data,y_data1,linewidth=0.5)
    plt.title("kick from stationary of pi/9")
    plt.ylabel("\u03B81")
    plt.xlabel("time")
    plt.show
    
    
def main(t2a,time):
    dt = 1/1000
    g = -9.81
    m1 = 0.762
    m2 = 8.891
    w1 = 0
    L1 = 1.516
    L2 = 0.191
    t1 = -np.pi/12
    kick_time = 10.405 + time  
    S = DoublePendulumLagrangian_kick(g, m1, m2, t1,w1,L1, L2,t2a,-t2a,kick_time)
    
    step = 0
    x_time = np.linspace(0,20,20000)
    y_t1 = []
    for i in range(20000):
        line = "%f,%f,%f,%f\n"
        y_t1.append(S.t1)
        sys.stdout.write(line % (step * dt, S.t1, S.t2, S.mechanical_energy()))
        S.time_step(dt)
        step += 1
    #plot(x_time,y_t1)
    #print(x_time[np.where(y_t1 == (min(y_t1[8000:10000])))])
    #print(ratio(x_time,y_t1))
    return ratio(x_time,y_t1)
    


def ratio(time,y):
    popt1, pcov1 = curve_fit(cos, time[:8000], y[:8000], p0=[10,2.5,1,-1])
    popt2, pcov2 = curve_fit(cos, time[12000:], y[12000:], p0=[10,2.5,1,-1])
    #plt.plot(x_time,cos(x_time,*popt), color = 'r',linewidth=1)
    ratio = np.abs(popt2[0])/np.abs(popt1[0])
    return(ratio)

def cos(t,A,T,x,y):
    return A*np.cos(2*np.pi*t/T + x) + y
    
def potential_energy2(t1,t2):
    g = 9.81
    m1 = 0.762
    m2 = 8.891
    L1 = 1.516
    L2 = 0.191
   
    y1 = -L1 * 0.5 * np.cos(t1)
    y2 = y1*2 - L2 * np.cos(t1+t2)
    return m1 * g * y1 + m2 * g * y2

def t1_start(t2a):
    x = np.linspace(-np.pi, np.pi, 100001)
    Z = potential_energy2(x,t2a)
    return(x[np.where(Z ==(min(Z)))])
# In[6]:



if __name__ == "__main__":
    start_time = time.time()
    t = np.linspace(-1.297,1.297,100)
    func = np.vectorize(main)
    x = func(np.pi/9,t)
    plot(t,x)
    print(max(x))
    print(t[np.where(x == max(x))])
    print(min(x))
    print(t[np.where(x == min(x))])
    print("--- %s seconds ---" % (time.time() - start_time))
