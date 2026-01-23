# -*- coding: utf-8 -*-
"""
Created on Fri Aug  3 20:03:54 2018

@author: Ahmad

This script calculates P,v, and T as given by the Van der Waals equation of
state for a single component.
All input units are SI.
"""
import numpy as np
import pandas as pd
import warnings

'''
Universal constants
'''
R = 8.314/100 # m^3*bar/(kmol*K)

class component(object):
    def __init__(self,name,MW,Tc,Pc):
        self.name = name
        self.MW = MW # g/mol
        self.Tc = Tc # K
        self.Pc = Pc # bar
        
        self.a = 27/64*(R*Tc)**2/Pc # m^6*bar/kmol^2
        self.b = R*Tc/(8*Pc) # m^3/kmol
    # VdW EOS implementations
    def P(self,v,T):
        # returns pressure at a given specific volume and temperature
        a = self.a; b = self.b
        return R*T/(v - b) - a/(v**2)
    
    def v(self,P,T):
        # returns specific volume at a given pressure and temperature
        a = self.a; b = self.b
        coeffs = [P, -(R*T + P*b), a, -a*b]
        sols = np.roots(coeffs)
        sols = np.real(sols[np.isreal(sols)]) # take only the real roots
        if len(sols) == 1:
            return [sols[0]]
        elif len(sols) != 0:
            print("Liquid-Vapour Region for component " + comp.name)
            return [min(sols), max(sols)]
        else:
            warnings.warn("Warning! No solution found for specific volume.")
            return
    
    def T(self,P,v):
        # returns the temperature at a given pressure and specific volume
        a = self.a; b = self.b
        return (P + a/(v**2))*(v - b)/R

# ideal gas EOS implementations
def idealP(v,T):
    return R*T/v
def idealv(P,T):
    return R*T/P
def idealT(P,v):
    return P*v/R
    
# sample test on ethane
ethane = component("ethane",30.069,305.4,48.74) # parameters from Koretsky
v = 10 # m^3/kmol
T = 273 # K
P = ethane.P(v,T)

print("Ideal gas pressure is: %.2f" % idealP(v,T))
print("VdW gas pressure is: %.2f" % ethane.T(P,v))

print("Ideal gas specific volume is: %.2f" % idealv(P,T))
print("VdW gas specific volume is: %.2f" % ethane.v(P,T)[-1])

print("Ideal gas temperature is: %.2f" % idealT(P,v))
print("VdW gas temperature is: %.2f" % ethane.T(P,v))