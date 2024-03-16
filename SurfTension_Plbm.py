#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Luiz Eduardo Czelusniak

Compute the theoretical planar interface surface tension for the 
pseudopotential method

In this code you enter with the Carnahan-Starling EOS parameters:
    a, b, R

Then, enter the reduced temperature in T and the free-energy kappa
parameter

Enter the planar interface equilibrium densities and pressure obtained
in the code Maxwell_Densities.py for the same EOS parameters and temperature:
    rho_v, rho_l, P_0
    
Enter de epsilon parameter "eps" obtained in the code Eps_Pseudopotential.py

Enter the kappa which is equivalent to the kappa_p parameter in the 
Modified Pseudopotential Method section of the paper

Enter the kappa_Li parameter to control surface tension

Enter the free energy surface tension gamma_FE

The code returns the surface tension gamma and the kappa_Li necessary to match
the free energy surface tension

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
T = 0.60 * T_c
# Interface width parameter
kappa = 0.8772583900683957
B = - kappa / 4
kappa_Li = -0.9483889821543865
factor = ( kappa - kappa_Li ) / 6
# Epsilon parameter
eps = 1.947427650061337
# Surface Tension of FE-LBM
gamma_FE = 0.01976871747796271

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

# Pseudopotential function
def PSI( rho, a, b, R, T ):
    G = -1
    c = 1
    cs2 = 1 / 3
    P_eos = EOS( rho, a, b, R, T )
    psi = np.sqrt( 2 * ( P_eos - rho * cs2 ) / G / c / c )
    return psi

# Derivative of pseudopotential function
def dPSI( rho, a, b, R, T ):
    G = -1
    c = 1
    cs2 = 1 / 3
    psi = PSI( rho, a, b, R, T )
    dP_eos = dEOS( rho, a, b, R, T )
    dpsi = 1 / psi / G / c ** 2 * ( dP_eos - cs2 )
    return dpsi 

# Base function
def func( rho, P_0, eps, a, b, R, T, B ):
    P_eos = EOS( rho, a, b, R, T )
    psi = PSI( rho, a, b, R, T )
    f = 2 * ( P_0 - P_eos ) / B / psi ** ( 1 + eps )
    return f

# Compute the integral of the base function
def INT_O2( P_0, rho_i, rho, eps, a, b, R, T, B, N ):
    psi = PSI( rho, a, b, R, T )
    rho_f = rho
    # Starting integration
    rho = rho_i 
    drho = ( rho_f - rho_i ) / N
    Int = 0
    Nc = round( N / 2 )
    # Computing integration
    for i in range(0,Nc):
        f_0 = func( rho, P_0, eps, a, b, R, T, B )
        f_1 = func( rho + drho, P_0, eps, a, b, R, T, B )
        f_2 = func( rho + 2 * drho, P_0, eps, a, b, R, T, B )
        dpsi_0 = dPSI( rho, a, b, R, T )
        dpsi_1 = dPSI( rho + drho, a, b, R, T )
        dpsi_2 = dPSI( rho + 2 * drho, a, b, R, T )
        Int = Int + ( f_0 * dpsi_0 + 4 * f_1 * dpsi_1 + f_2 * dpsi_2 ) / 3 * drho
        rho = rho + 2 * drho 
    F = psi ** eps * Int
    return F

# Function that compute the density profile
def dPSI_dx( rho_v, rho, P_0, eps, a, b, R, T, B, Ni ):
    FuncInt = INT_O2( P_0, rho_v, rho, eps, a, b, R, T, B, Ni )
    dpsi_dx = np.sqrt( np.abs( FuncInt ) )
    return dpsi_dx

# Function that compute the density profile
def SURF_TENSION( rho_v, rho_l, P_0, eps, a, b, R, T, B, Np, Ni ):
    # Starting integration
    rho = rho_v 
    drho = ( rho_l - rho_v ) / Np
    Int = 0
    Nc = round( Np / 2 )
    # Computing integration
    for i in range(0,Nc):
        f_0 = dPSI_dx( rho_v, rho, P_0, eps, a, b, R, T, B, Ni )
        f_1 = dPSI_dx( rho_v, rho + drho, P_0, eps, a, b, R, T, B, Ni )
        f_2 = dPSI_dx( rho_v, rho + 2 * drho, P_0, eps, a, b, R, T, B, Ni )
        dpsi_0 = dPSI( rho, a, b, R, T )
        dpsi_1 = dPSI( rho + drho, a, b, R, T )
        dpsi_2 = dPSI( rho + 2 * drho, a, b, R, T )
        Int = Int + ( f_0 * dpsi_0 + 4 * f_1 * dpsi_1 + f_2 * dpsi_2 ) / 3 * drho
        rho = rho + 2 * drho 
    gamma = factor * Int
    return gamma
    

gamma = SURF_TENSION( rho_v, rho_l, P_0, eps, a, b, R, T, B, 200, 4000 )
print( "gamma: ", gamma )

kappa_Li = kappa - 6 * gamma_FE * factor / gamma
print( "kappa_L: ", kappa_Li )

