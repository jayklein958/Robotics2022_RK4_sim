#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 25 11:26:50 2022

@author: jonathanklein
"""


import numpy 
import numpy as np

class DoublePendulumLagrangian_kick:

    def __init__(self, g, m1, m2, t1,w1,L1, L2,t2a,t2b,kick_time):
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
        self.kt = kick_time*1000
        self.g = g
        self.time = 0
        self.m1 = m1
        self.m2 = m2
        self.t1 = t1
        self.t2 = t2a
        self.w1 = w1
        self.L1 = L1
        self.L2 = L2
        self.t2a = t2a
        self.t2b = t2b

    def potential_energy(self):
        """Computes the potential energy of the system."""
        # compute the height of each bob
        y1 = -self.L1 * 0.5 * np.cos(self.t1)
        y2 = y1*2 - self.L2 * np.cos(self.t2+self.t2)

        return self.m1 * self.g * y1 + self.m2 * self.g * y2

    def kinetic_energy(self):
        """Computes the kinetic energy of the system."""
        # compute the kinetic energy of each bob
        K1 = 0.5*self.m1*(self.w1)**2*(self.L1)**2
        K2 = 0.5*self.m2*(self.w1)**2*((self.L1)**2+(self.L2)**2+2*self.L1*self.L2*np.cos(self.t2))

        return K1 + K2

    def mechanical_energy(self):
        """
        Computes the mechanical energy (total energy) of the
        system.
        """

        return self.kinetic_energy() + self.potential_energy()
    


    def lagrange_rhs(self, t1, w1):
        """
        Computes the right-hand side of the Euler-Lagrange equations
        for the double pendulum and returns it as an array.
        t1 - The angle of bob #1.
        t2 - The angle of bob #2.
        w1 - The angular velocity of bob #1.
        w2 - The angular velocity of bob #2.
        """
        g1 = (self.m1*0.5 + self.m2)*self.g*self.L1*np.sin(t1) + self.m2*self.g*self.L2*np.sin(t1+self.t2)
        g2 = 0.25*self.L1*self.m1+self.m2*((self.L1)**2+(self.L2)**2+2*self.L1*self.L2*np.cos(self.t2))

        #return numpy.array([self.w1,(g1/g2 - 2*0.006*self.w1)])
        return numpy.array([self.w1,(g1/g2)])

    def time_step(self, dt):
        """
        Advances one time step using RK4 (classical Runge-Kutta
        method).
        """
        # y is an array with the generalized coordinates (angles +
        # angular velocities)
        kick_time = np.abs(self.t2b - self.t2a)/6.6
        kick_angle = (self.t2b - self.t2a)/(kick_time*1/dt)
        if self.kt > self.time:
            self.t2 = self.t2a      
        if self.kt < self.time < (self.kt + (kick_time*1/dt)):
            self.t2 += kick_angle
            
        #self.t1 = (self.t1)*np.exp(-0.0000001*self.time/1000)
        #self.w1 = (self.w1)*np.exp(-0.0000001*self.time/1000)
        #*np.exp(-(0.0000009)*(self.time)*1/1000)
        
        y = numpy.array([self.t1,self.w1])
        # compute the RK4 constants
        k1 = self.lagrange_rhs(*y)
        k2 = self.lagrange_rhs(*(y + dt * k1 / 2))
        k3 = self.lagrange_rhs(*(y + dt * k2 / 2))
        k4 = self.lagrange_rhs(*(y + dt * k3))

        # compute the RK4 right-hand side
        R = 1.0 / 6.0 * dt * (k1 + 2.0 * k2 + 2.0 * k3 + k4)

        # update the angles and angular velocities
        
        self.t1 += R[0]
        self.w1 += R[1]
        #self.t1=self.t1*np.exp(-(0.007)*self.time/1000)
        #self.w1=self.w1*(0.06/1000)*np.exp(-(0.06)*self.time/1000)
        self.time +=1