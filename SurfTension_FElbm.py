#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Luiz Eduardo Czelusniak

Compute the theoretical planar interface surface tension for the 
free-energy method

In this code you enter with the Carnahan-Starling EOS parameters:
    a, b, R

Then, enter the reduced temperature in T and the free-energy kappa
parameter

Enter the planar interface equilibrium densities and pressure obtained
in the code Maxwell_Densities.py for the same EOS parameters and temperaute:
    rho_v, rho_l, P_0

Run the code to obtain the surface tension (in lattice units)
"""

import matplotlib.pyplot as plt
import numpy as np

plt.close('all')

# Equation of State Parameters
a = 0.25
b = 4
R = 1
# Computing Reduced Properties
T_c = a / b / R / 2.6502
P_c = 0.18727 ** 2 / 0.4963 * a / b / b
rho_c = 0.5218 / b
# Defining Function Range
rho_v = 0.0006268162514430721
rho_l = 0.45406865373315114
P_0 = 7.311478997389506e-06
T = 0.50 * T_c
# Interface width parameter
kappa = 0.5

# Function to compute the equation of state
def EOS( rho, a, b, R, T ):
    n = b * rho / 4
    num = ( 1 + n + n ** 2 - n ** 3 )
    den = ( 1 - n ) ** 3
    aux = num / den
    Peos = rho * R * T * aux - a * rho ** 2
    return Peos

# Function to compute the derivative of EOS
def dEOS( rho, a, b, R, T ):
    n = b * rho / 4
    T1 = rho * R * T
    T2 = 1 + n + n ** 2 - n ** 3
    T3 = ( 1 - n ) ** 3
    dT1 = R * T
    dT2 = b / 4 + 2 * n * b / 4 - 3 * n * n * b / 4
    r_dT3_T3 = 3 * ( - b / 4 ) / ( 1 - n ) ** 4
    dP = dT1 * T2 / T3 + T1 * dT2 / T3 - T1 * T2 * r_dT3_T3 - 2 * a * rho
    return dP

# Base function
def func( rho, P_0, a, b, R, T ):
    P_eos = EOS( rho, a, b, R, T )
    f = 2 * ( P_eos - P_0 ) / rho / rho
    return f

# Compute the integral of the base function
def INT_O2( P_0, rho_i, rho, a, b, R, T, N ):
    rho_f = rho
    # Starting integration
    rho = rho_i 
    drho = ( rho_f - rho_i ) / N
    Int = 0
    Nc = round( N / 2 )
    # Computing integration
    for i in range(0,Nc):
        f_0 = func( rho, P_0, a, b, R, T )
        f_1 = func( rho + drho, P_0, a, b, R, T )
        f_2 = func( rho + 2 * drho, P_0, a, b, R, T )
        Int = Int + ( f_0 + 4 * f_1 + f_2 ) / 3 * drho
        rho = rho + 2 * drho 
    F = rho / kappa * Int
    return F

# Function that compute the density profile
def dRHO_dx( rho_v, rho, P_0, a, b, R, T, Ni ):
    FuncInt = INT_O2( P_0, rho_v, rho, a, b, R, T, Ni )
    drho_dx = np.sqrt( np.abs( FuncInt ) )
    return drho_dx

# Function that compute the density profile
def SURF_TENSION( rho_v, rho_l, P_0, a, b, R, T, Np, Ni ):
    # Starting integration
    rho = rho_v 
    drho = ( rho_l - rho_v ) / Np
    Int = 0
    Nc = round( Np / 2 )
    # Computing integration
    for i in range(0,Nc):
        f_0 = dRHO_dx( rho_v, rho, P_0, a, b, R, T, Ni )
        f_1 = dRHO_dx( rho_v, rho + drho, P_0, a, b, R, T, Ni )
        f_2 = dRHO_dx( rho_v, rho + 2 * drho, P_0, a, b, R, T, Ni )
        Int = Int + ( f_0 + 4 * f_1 + f_2 ) / 3 * drho
        rho = rho + 2 * drho 
    gamma = kappa * Int
    return gamma
    

gamma = SURF_TENSION( rho_v, rho_l, P_0, a, b, R, T, 200, 4000 )
print( "gamma: ", gamma )