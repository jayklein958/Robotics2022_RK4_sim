#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 25 09:58:44 2022

@author: jonathanklein
"""

import math
import numpy as np
import pygame
import time
import sys

class  DoublePendulumLagrangian:

    def __init__(self, g, m1, m2, t1, t2, w1, w2, L1, L2):
        """
        Constructs a double pendulum simulator based on its
        Euler-Lagrange equations. Bob #1 is the one attached to the
        fixed pivot.
        g - The gravitational acceleration.
        m1 - The mass of bob #1.
        m2 - The mass of bob #2.
        t1 - The initial angle of bob #1.
        t2 - The initial angle of bob #2.
        w1 - The initial angular velocity of bob #1.
        w2 - The initial angular velocity of bob #2.
        L1 - The length of the rod for bob #1.
        L2 - The length of the rod for bob #2.
        """
        self.g = g
        self.m1 = m1
        self.m2 = m2
        self.t1 = t1
        self.t2 = t2
        self.w1 = w1
        self.w2 = w2
        self.L1 = L1
        self.L2 = L2

    def potential_energy(self):
        """Computes the potential energy of the system."""
        y1 = -self.L1 * 0.5 * np.cos(self.t1)
        y2 = y1*2 - self.L2 * np.cos(self.t2)

        return self.m1 * self.g * y1 + self.m2 * self.g * y2

    def kinetic_energy(self):
        """Computes the kinetic energy of the system."""
        K1 = 0.5 * self.m1 * (0.5*self.L1 * self.w1)**2
        K2 = 0.5 * self.m2 * ((self.L1 * self.w1)**2 + (self.L2 * self.w2)**2 +
                         2 * self.L1 * self.L2 * self.w1 * self.w2 * np.cos(self.t1 - self.t2))

        return K1 + K2

    def mechanical_energy(self):
        """
        Computes the mechanical energy (total energy) of the
        system.
        """

        return self.kinetic_energy() + self.potential_energy()


    def lagrange_rhs(self, t1, t2, w1, w2):
        """
        Computes the right-hand side of the Euler-Lagrange equations
        for the double pendulum and returns it as an array.
        t1 - The angle of bob #1.
        t2 - The angle of bob #2.
        w1 - The angular velocity of bob #1.
        w2 - The angular velocity of bob #2.
        """
        a1 = (self.L2 / self.L1) * (4*self.m2 / (self.m1 + 4*self.m2)) * math.cos(self.t1 - self.t2)
        a2 = (self.L1 / self.L2) * math.cos(self.t1 - self.t2)
        f1 = -(self.L2 / self.L1) * (4*self.m2 / (self.m1 + 4*self.m2)) * (self.w2**2) * math.sin(self.t1 - self.t2) - (self.g / self.L1)* ((2*self.m1+4*self.m2)/(self.m1+4*self.m2))* math.sin(self.t1)
        f2 = (self.L1 / self.L2) * (self.w1**2) * math.sin(self.t1 - self.t2) - (self.g / self.L2) * math.sin(self.t2)

        g1 = (f1 - a1 * f2) / (1 - a1 * a2)
        g2 = (f2 - a2 * f1) / (1 - a1 * a2)
        

        return np.array([self.w1, self.w2, g1 - 2*0.006*self.w1,g2])

    def time_step(self, dt):
        """
        Advances one time step using RK4 (classical Runge-Kutta
        method).
        """

        # y is an array with the generalized coordinates (angles +
        # angular velocities)
        y = np.array([self.t1, self.t2, self.w1, self.w2])

        # compute the RK4 constants
        k1 = self.lagrange_rhs(*y)
        k2 = self.lagrange_rhs(*(y + dt * k1 / 2))
        k3 = self.lagrange_rhs(*(y + dt * k2 / 2))
        k4 = self.lagrange_rhs(*(y + dt * k3))

        # compute the RK4 right-hand side
        R = 1.0 / 6.0 * dt * (k1 + 2.0 * k2 + 2.0 * k3 + k4)

        # update the angles and angular velocities
        self.t1 += R[0]
        self.t2 += R[1]
        self.w1 += R[2]
        self.w2 += R[3]