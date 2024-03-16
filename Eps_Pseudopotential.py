#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Luiz Eduardo Czelusniak

Compute the epsilon parameter of pseudopotential method to match the same
equilibrium densities given by the Maxwell-rule in a planar interface

In this code you enter with the Carnahan-Starling EOS parameters:
    Var.a, Var.b, Var.R

Then, enter the reduced temperature in T 

Enter the planar interface equilibrium densities and pressure obtained
in the code Maxwell_Densities.py for the same EOS parameters and temperaute:
    rho_v, rho_l, P_0

This code must be in the same folder as Miscellaneous.py
"""

import matplotlib.pyplot as plt
import numpy as np
from Miscellaneous import CS_EOS

plt.close('all')

# Function to compute the pseudopotential
def PSI( rho, T ):
    # Importing EOS class
    #Var = CS_EOS()
    # Auxiliary for pseudopotential
    G = -1
    c = 1
    cs2 = 1 / 3
    P_eos = Var.P_from_rho( rho, T )
    # Computing pseudopotential
    psi = np.sqrt( 2 * ( P_eos - rho * cs2 ) / G / c / c  )
    return psi

def dPSI( rho, T ):
    # Importing EOS class
    #Var = CS_EOS()
    # Auxiliary for pseudopotential
    G = -1
    c = 1
    cs2 = 1 / 3
    dP_eos = Var.dP_from_rho( rho, T )
    psi = PSI( rho, T )
    # Computing pseudopotential
    dpsi = 1 / G / c ** 2 / psi * ( dP_eos - cs2 )  
    return dpsi

# Base function to be integrated
def func( rho, P_0, T, eps ):
    # Importing EOS class
    #Var = CS_EOS()
    # Auxiliary computations
    psi = PSI( rho, T )
    P_eos = Var.P_from_rho( rho, T )
    with np.errstate(all='raise'): # Avoid division by zero
        try:
            f = 2* ( P_0 - P_eos ) / psi ** ( 1 + eps )
        except:
            f = 2 * ( P_0 -  P_eos ) 
    return f
  
# Base Function for error calculation
def base_Func( rho, P_0, T, B ):
    # Importing EOS class
    #Var = CS_EOS()
    # Auxiliary computations
    psi = PSI( rho, T )
    P_eos = Var.P_from_rho( rho, T )
    with np.errstate(all='raise'): # Avoid division by zero
        try:
            f = 2 * ( P_0 - P_eos ) / B / psi
        except:
            f = 2 * ( P_0 -  P_eos ) 
    return f

"""
# Integral O2
def Int_O2( rho_f, rho_i, P_0, eps, T, B, N ):
    # Preparing the integral
    rho = rho_i
    drho = ( rho_f - rho_i ) / N
    Int = 0
    Nd = round( N / 2 )
    for i in range( 0, Nd ):
        f_0 = func( rho, P_0, T, eps ) / B
        f_1 = func( rho + drho, P_0, T, eps ) / B
        f_2 = func( rho + 2 * drho, P_0, T, eps ) / B
        dpsi_0 = dPSI( rho, T )
        dpsi_1 = dPSI( rho + drho, T )
        dpsi_2 = dPSI( rho + 2 * drho, T )
        Int = Int + ( f_0 * dpsi_0 + 4 * f_1 * dpsi_1 + f_2 * dpsi_2 ) / 3 * drho
        rho = rho + 2 * drho
    return Int
"""
# Compute the integral of the base function
def Rectangular( P_0, rho_i, rho_f, eps, T, N ):
    # Starting integration
    rho = rho_i
    drho = ( rho_f - rho_i ) / N
    Int = 0
    for i in range(0,N):
        f = func( rho, P_0, T, eps )
        dpsi = dPSI( rho, T )
        Int = Int + f * dpsi * drho
        rho = rho + drho
    return Int

# Compute the integral of the base function
def Simpson( P_0, rho_i, rho_f, eps, T, N ):
    N = N + 1
    # Starting integration
    rho = rho_i
    drho = ( rho_f - rho_i ) / ( N - 1 )
    Int = 0
    Nd = round( (N+1) / 2 )
    for i in range( 0, Nd+1 ):
        f_0 = func( rho, P_0, T, eps )
        f_1 = func( rho + drho, P_0, T, eps )
        f_2 = func( rho + 2 * drho, P_0, T, eps )
        dpsi_0 = dPSI( rho, T )
        dpsi_1 = dPSI( rho + drho, T )
        dpsi_2 = dPSI( rho + 2 * drho, T )
        Int = Int + ( f_0 * dpsi_0 + 4 * f_1 * dpsi_1 + f_2 * dpsi_2 ) / 3 * drho
        rho = rho + 2 * drho
    return Int

#*****************************************************************************
# Code Start

# Initialization

# Creating an object Var to store all variables
Var = CS_EOS()
# Equation of State Parameters
Var.a = 0.25
Var.b = 4
Var.R = 1
# Computing Reduced Properties
T_c = Var.a / Var.b / Var.R / 2.6502
P_c = 0.18727 ** 2 / 0.4963 * Var.a / Var.b ** 2
rho_c = 0.5218 / Var.b
Var.T_c = T_c
Var.P_c = P_c
Var.rho_c = rho_c
# Defining Function Range
rho_v = 0.02172865606779925
rho_l = 0.30717839449894835
P_0 = 0.00032955400929036147
T = 0.80 * T_c

# Integral less accurate
eps = 2
F = Rectangular( P_0, rho_v, rho_l, eps, T, 1000 )
print(F)

tol = 1e-8
deps = 1

while( ( abs(deps) > tol ) | ( abs(F) > tol ) ):
    delta = 1e-8
    Fn = Simpson( P_0, rho_v, rho_l, eps + delta, T, 50000 )
    F = Simpson( P_0, rho_v, rho_l, eps, T, 50000 )
    dF = ( Fn - F ) / delta
    deps = - F / dF
    eps = eps + deps
    print("deps: ", deps)

print( "Int: ", F )
print("eps: ", eps)






