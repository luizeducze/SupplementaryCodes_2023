#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Luiz Eduardo Czelusniak

Compute the theoretical droplet thermodynamically consistent equilibrium 
densities

In this code you enter with the Carnahan-Starling EOS parameters:
    EOS.a, EOS.b, EOS.R

Then, just control the reduced temperature in Dens.T and the number of 
points for integration INT.N. Then, define the droplet radius and surface 
tension with variables "radius" and "gamma".
After that run the code to obtain the equilibrium densities:
    rho_v, rho_l

This code must be in the same folder as Miscellaneous.py
"""

import matplotlib.pyplot as plt
import numpy as np
from Miscellaneous import CS_EOS, Integration

# We give this name because is the code that compute the densities
class DensityFunc:
    
    # Variables
    T = 1
    
    # Maximum and minimum
    rho_max = 1
    rho_min = 0
    
    # Pointer to the desire EOS
    pEOS = 0
    
    # Function we want to integrate
    def DensFunc( self, rho ):
        # Defining necessary parameters
        T = self.T
        P_eos = self.pEOS.P_from_rho( rho, T )
        # Computing the function
        with np.errstate(all='raise'): # Avoid division by zero
            try:
                f = ( P_eos ) / rho / rho
            except:
                f = ( P_eos ) 
        return f
    
#****************************************************************************

# Find the Condition of the Maxwell Rule
def MAXWELL( EOS, Dens, INT, gamma, radius ):
    # Defining range of search
    rho_max = Dens.rho_max
    P_MAX = EOS.P_from_rho( rho_max, T )
    P_MIN = 0
    # Bissection Method
    tol = 1e-8
    Int = 1
    while ( P_MAX - P_MIN > tol ) | ( abs(Int) > tol ):
        P_out = ( P_MAX + P_MIN ) / 2
        P_in = P_out + gamma / radius
        # Updating the integration limits
        INT.x_i = EOS.rhov_from_P( P_out, Dens.T )
        INT.x_f = EOS.rhol_from_P( P_in, Dens.T )
        # Integrating the base function
        Int = INT.Simpson( Dens.DensFunc ) + P_in / INT.x_f - P_out / INT.x_i
        # Evaluating the results
        if Int > 0:
            P_MIN = P_out
        else:
            P_MAX = P_out
    rho_v = EOS.rhov_from_P( P_out, Dens.T )
    rho_l = EOS.rhol_from_P( P_in, Dens.T )
    return P_out, P_in, rho_v, rho_l
    
#****************************************************************************
# Code Start

# Initialization

# Storing EOS information
EOS = CS_EOS()
# Equation of state parameters
EOS.a = 0.25
EOS.b = 4
EOS.R = 1
# Computing reduced properties
T_c = EOS.a / EOS.b / EOS.R / 2.6502
P_c = 0.18727 ** 2 / 0.4963 * EOS.a / EOS.b ** 2
rho_c = 0.5218 / EOS.b
EOS.T_c = T_c
EOS.P_c = P_c
EOS.rho_c = rho_c

# Droplet information
gamma = 0.004633485900242 # Surface tension
radius = 40 # Droplet radius

# Storing function information
Dens = DensityFunc()
Dens.T = 0.60 * T_c # Control the reduced temperature
Dens.pEOS = EOS
#
T = Dens.T
rho_i = 0.1*rho_c
rho_f = 3.2*rho_c

# Defining max and min for the EOS
Dens.rho_max = EOS.MAX_EOS( T ) # point of maximum for EOS
print("rho_max: ", Dens.rho_max)
Dens.rho_min = EOS.MIN_EOS( T ) # point of minimum for EOS
print("rho_min: ", Dens.rho_min)

# Defining some dynamics properties
rho = np.linspace( rho_i, rho_f, 100 )
P_eos = EOS.P_from_rho( rho, T )
Dens.P_0 = EOS.P_from_rho( rho_i, T )
P_0 = Dens.P_0

# Defining integration properties
INT = Integration()
INT.N = 20000 # Number of points for integration


P_out, P_in, rho_v, rho_l = MAXWELL( EOS, Dens, INT, gamma, radius )
print("P_out: ", P_out)
print("P_in: ", P_in)
print("rho_v: ", rho_v)
print("rho_l: ", rho_l)

Int_REC = INT.Rectangular( Dens.DensFunc ) + P_in / INT.x_f - P_out / INT.x_i
print("Integral O1: ", Int_REC)

Int_SMP = INT.Simpson( Dens.DensFunc ) + P_in / INT.x_f - P_out / INT.x_i
print("Integral O3: ", Int_SMP)

# Teste
print( "P_in - P_out: ", P_in - P_out )
print( "gamma/radius: ", gamma / radius )

plt.plot( rho, P_eos )
plt.ylabel('Pressure')
plt.xlabel('density')
# specifying horizontal line type
plt.axhline(y = P_in, color = 'r', linestyle = '-')
plt.axhline(y = P_out, color = 'r', linestyle = '-')
# mark in liquid density
plt.plot( rho_v, P_out, marker='o', color = 'k' )
plt.plot( rho_l, P_in, marker='o', color = 'k' )
plt.show()
