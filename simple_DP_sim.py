#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 25 10:09:01 2022

@author: jonathanklein
"""
import numpy as np
import pygame
import time
import sys
import matplotlib.pyplot as plt

from Simple_DP_RK4 import DoublePendulumLagrangian

def plot(x_data,y_data,x_title,y_title,title):
    plt.plot(x_data,y_data,linewidth=0.5)
    plt.xlabel(x_title)
    plt.ylabel(y_title)
    plt.title(title)
    plt.show

def draw(S,window,Nx,Ny,dt):
    """    
    Draws the double pendulum system on a window.
    S - The double pendulum object.
    window - The window where the double pendulum will be shown.
    Nx - The window width (in pixels).
    Ny - The window height (in pixels).
    dt - The simulation time step.
    """
    
    m1 = S.m1  # m1 - The mass of bob #1.
    m2 = S.m2  # m2 - The mass of bob #2
    t1 = S.t1  # t1 - The initial angle of bob #1.
    t2 = S.t2  # t2 - The initial angle of bob #2.
    L1 = S.L1  # L1 - The length of the rod for bob #1.
    L2 = S.L2  # L2 - The length of the rod for bob #2

    # radius (in pixels) of each bob (min/max: 3/12 pixels)
    #R1 = max(3, int(12 * (m1 / (m1 + m2))))
    R1 = 0
    R2 = max(3, int(12 * (m2 / (m1 + m2))))

    # length (in pixels) of each rod
    P1 = (0.85 * min(Nx / 2, Ny / 2) * (L1 / (L1 + L2)))
    P2 = 0.85 * min(Nx / 2, Ny / 2) * (L2 / (L1 + L2))

    # positions (in (pixels,pixels)) of each bob
    X0 = np.array([int(Nx / 2), int(Ny / 2)])
    X1 = X0 + np.array([int(P1 * np.sin(t1)), int(P1 * np.cos(t1))])
    X2 = X1 + np.array([int(P2 * np.sin(t2)), int(P2 * np.cos(t2))])

    # color: rods and bobs
    color_L1 = (0,0,0)
    color_L2 = (128, 128, 128)
    color_m1 = (255, 0, 0)
    color_m2 = (0, 0, 255)

    # clear the window
    window.fill((255, 255, 255))

    # draw the rods and the bobs
    pygame.draw.line(window, color_L1, X0, X1, 3)
    pygame.draw.line(window, color_L2, X1, X2, 3)
    pygame.draw.circle(window, color_m1, X1, int(R1))
    pygame.draw.circle(window, color_m2, X2, int(R2))

    # write the time step value on the window
    myfont = pygame.font.SysFont("Arial", 15)
    label = myfont.render("dt = %.3g" % dt, 1, (128, 128, 128))
    window.blit(label, (10, 10))

    # update the screen
    pygame.display.flip()
    
def main1():
    """
    Live simulation using Pygame of a double pendulum
    """
    # default simulation parameter values
    dt = 1/200
    g = 9.81
    m1 = 0.762
    m2 = 8.891
    L1 = 1.516
    L2 = 0.191
    t1 = 1
    t2 = 2
    w1 = 0
    w2 = 0

    # default window dimensions
    Nx = Ny = 500
    
    S = DoublePendulumLagrangian(g, m1, m2, t1, t2, w1, w2, L1, L2)
    
    # E0 = initial mechanical energy of the system
    step = 0
    # maximum energy change (compared to E0): too large => unstable simulation
    pygame.init()
    clock = pygame.time.Clock()
    window = pygame.display.set_mode((Nx, Ny), pygame.RESIZABLE)
    pygame.display.set_caption("double pendulum")
    
    # keep running the simulation until the user closes the window
    for i in range(10000):

        # redraw the double pendulum at a maximum rate of 25 fps
        draw(S, window, Nx, Ny, dt)
        clock.tick(200)
        line = "%f,%f,%f,%f\n"
        sys.stdout.write(line % (step * dt, S.t1, S.t2, S.mechanical_energy()))
        S.time_step(dt)
        step += 1    

def main2(t1,t2):
    """Plotting simulation of the simple double pendulum with damping long run"""
    # default simulation parameter values
    g = 9.81
    dt = 1/1000
    m1 = 0.762
    m2 = 8.891
    w1 = 0
    w2 = 0
    L1 = 1.516
    L2 = 0.191
    

   
    S = DoublePendulumLagrangian(g, m1, m2, t1, t2, w1, w2, L1, L2)
   
    step = 0
    x_time = np.linspace(0,100,100000)
    y_t1a = []
    

    for i in range(100000):
        line = "%f,%f,%f,%f\n"
        #x_time.append(step * dt)
        y_t1a.append(S.t1)
        sys.stdout.write(line % (step * dt, S.t1, S.t2, S.mechanical_energy()))
        S.time_step(dt)
        step += 1
    x_title = "time"
    y_title = "t1"
    title = "damped pendulum over time"
    plot(x_time, y_t1a, x_title, y_title, title)
    

      
    
    
if __name__ == "__main__":
    start_time = time.time()
    #main1()
    main2(1,0)
    print("--- %s seconds ---" % (time.time() - start_time))