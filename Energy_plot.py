#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 25 12:35:57 2022

@author: jonathanklein
"""
import numpy as np
import matplotlib.cm as cm
import matplotlib.pyplot as plt


def fdp_energy(t1,t2):
    g = 9.81
    m1 = 0.762
    m2 = 8.891
    w1 = 0
    w2 = 0
    L1 = 1.516
    L2 = 0.191
    y1 = -L1 * 0.5 * np.cos(t1)
    y2 = y1*2 - L2 * np.cos(t2)
    K1 = 0.5 * m1 * (0.5*L1 * w1)**2
    K2 = 0.5 * m2 * ((L1 * w1)**2 + (L2 * w2)**2 + 2 * L1 * L2 * w1 * w2 * np.cos(t1 - t2))
    
    return (K1 + K2 + m1 * g * y1 + m2 * g * y2)/(m1*g*L1)


print(fdp_energy(0,0))
print(fdp_energy(np.pi,np.pi))
print(fdp_energy(np.pi,0))
print(fdp_energy(0,np.pi))


def plot(fx):
    x = np.linspace(-np.pi, np.pi, 100)
    y = np.linspace(-np.pi, np.pi, 100)

    X, Y = np.meshgrid(x, y)
    Z = fx(X, Y)
    levels  = [-13,-12,-10.69793679319109,-9,-5,-1,3,7,10.69793679319109,13]
    contour = plt.contour(X, Y, Z,levels=levels,cmap=cm.plasma)
    plt.clabel(contour, inline=True, fontsize=7)
    plt.xlabel('t1')
    plt.ylabel('t2')
    plt.title("Energy projection onto t1-t2 plane for Double pendulum ")
    plt.colorbar()
    

plot(fdp_energy)