#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Luiz Eduardo Czelusniak

Compute the theoretical density profile of pseudopotential method for
a planar interface

In this code you enter with the Carnahan-Starling EOS parameters:
    a, b, R

Then, enter the reduced temperature in T 

Enter the planar interface equilibrium densities and pressure obtained
in the code Maxwell_Densities.py for the same EOS parameters and temperature:
    rho_v, rho_l, P_0
    
Enter de epsilon parameter "eps" obtained in the code Eps_Pseudopotential.py

Enter the kappa which is equivalent to the kappa_p parameter in the 
Modified Pseudopotential Method section of the paper

This code must be in the same folder as Miscellaneous.py
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
kappa = 1
B = - kappa / 4
# Epsilon parameter
eps = 1.947427650061337

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
def PROFILE( rho_v, rho_l, P_0, eps, a, b, R, T, B, L, Np, Ni ):
    rho_m = ( rho_v + rho_l ) / 2
    x = np.linspace( -L, L, Np + 1 )
    rho = np.linspace( rho_v, rho_l, Np + 1 )
    dx = 2 * L / Np
    Nc = round( Np / 2 )
    rho[ Nc ] = rho_m
    for i in range(0,Nc):
        F_l = INT_O2( P_0, rho_v, rho[ Nc - i ], eps, a, b, R, T, B, Ni )
        F_r = INT_O2( P_0, rho_v, rho[ Nc + i ], eps, a, b, R, T, B, Ni )
        rho[ Nc - i - 1 ] = rho[ Nc - i ] - dx * np.sqrt( np.abs( F_l ) )
        rho[ Nc + i + 1 ] = rho[ Nc + i ] + dx * np.sqrt( np.abs( F_r ) )
        if( rho[ Nc + i + 1 ] > rho_l ):
            rho[ Nc + i + 1 ] = rho_l
        if( rho[ Nc - i - 1 ] < rho_v ):
            rho[ Nc - i - 1 ] = rho_v
    return x, rho

# Function that compute the interface width
def WIDTH( x, rho, rho_v, rho_l, Np ):
    rho_min = rho_v
    rho_max = rho_l
    rho_nd = ( rho - rho_min ) / ( rho_max - rho_min )
    for i in range( 0, Np ):
        if( rho_nd[ i ] < 0.12 ) & ( rho_nd[ i + 1 ] > 0.12 ):
            dx = x[ i + 1 ] - x[ i ]
            delta = ( 0.12 - rho_nd[ i ] ) / ( rho_nd[ i + 1 ] - rho_nd [ i ] )
            x_init = x[ i ] + dx * delta
        if( rho_nd[ i ] < 0.88  ) & ( rho_nd[ i + 1 ] > 0.88 ):
            dx = x[ i + 1 ] - x[ i ]
            delta = ( 0.88 - rho_nd[ i ] ) / ( rho_nd[ i + 1 ] - rho_nd [ i ] )
            x_final = x[ i ] + dx * delta
    width = x_final - x_init
    return width, x_init, x_final

# Function that finds the required kappa value
def FIND_B( W, rho_v, rho_l, P_0, eps, a, b, R, T, B, L, Np, Ni ):
    B_try = B
    delta_B = 1
    tol = 1e-6
    dB = 1e-8
    while( delta_B > tol ):
        x, rho = PROFILE( rho_v, rho_l, P_0, eps, a, b, R, T, B_try, L, Np, Ni )
        w, x_i, x_f = WIDTH( x, rho, rho_v, rho_l, Np )
        F = W - w
        x, rho = PROFILE( rho_v, rho_l, P_0, eps, a, b, R, T, B_try + dB, L, Np, Ni )
        w, x_i, x_f = WIDTH( x, rho, rho_v, rho_l, Np )
        F_n = W - w
        dF = ( F_n - F ) / dB
        delta_B = - F / dF
        B_try = B_try + delta_B
        print( "B_try: ", B_try )
    return B_try
        
W = 3.4903942124364757
B_new = FIND_B( W, rho_v, rho_l, P_0, eps, a, b, R, T, B, L = 8, Np = 200, Ni = 4000 )
x, rho = PROFILE( rho_v, rho_l, P_0, eps, a, b, R, T, B_new, L = 8, Np = 200, Ni = 4000 )
w, x_i, x_f = WIDTH( x, rho, rho_v, rho_l, Np = 200 )
print( "w: ", w )
print( "kappa: ", 4 * B_new )

plt.plot( x, rho )
plt.plot( x_i, rho_v + 0.12*( rho_l - rho_v ), marker='o', color = 'k' )
plt.plot( x_f, rho_v + 0.88*( rho_l - rho_v ), marker='o', color = 'k' )
plt.axhline(y = rho_v, color = 'r', linestyle = '--')
plt.axhline(y = rho_l, color = 'r', linestyle = '--')
plt.ylabel('density')
plt.xlabel('position')
plt.show()

psi = PSI( rho, a, b, R, T )
plt.plot( rho, psi )
plt.ylabel('Pressure')
plt.xlabel('density')
plt.show()









