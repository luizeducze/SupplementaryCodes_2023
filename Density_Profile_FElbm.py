#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Luiz Eduardo Czelusniak

Compute the theoretical planar interface density profile for the 
free-energy method

In this code you enter with the Carnahan-Starling EOS parameters:
    EOS.a, EOS.b, EOS.R

Then, enter the reduced temperature PROF.T and the free-energy kappa
parameter PROF.kappa

Enter the planar interface equilibrium densities and pressure obtained
in the code Maxwell_Densities.py for the same EOS parameters and temperaute:
    rho_v, rho_l, P_0

Run the code to obtain the density profile and interface thickness 
(in lattice units):
    w_O2 (interface thickness computed with second order of accuracy) 
    w_O4 (interface thickness computed with fourth order of accuracy)

This code must be in the same folder as Miscellaneous.py
"""

import matplotlib.pyplot as plt
import numpy as np
from Miscellaneous import CS_EOS, Integration

plt.close('all')

#*****************************************************************************
# Classes

# This class will be used to compute the density profile
class ProfileClass:
    # Variables
    T = 1 # Temperature
    mu_0 = 1 # bulk chemical potential
    kappa = 1 # surface tension parameter
    
    # Pointer to the desired EOS object
    pEOS = 0
    
    def ProfileFunc( self, rho ):
        # chemical potential
        mu = self.pEOS.Chemical( rho, self.T )
        # profile function
        f = 2 * ( mu - self.mu_0 ) / self.kappa
        return f

#*****************************************************************************

# Density profile of order O(dx)^2
def ProfileO2( EOS, PROF, INT, Np, L, rho_v, rho_l ):
    # Computations
    rho_m = ( rho_v + rho_l ) / 2
    x = np.linspace( -L, L, Np + 1 )
    rho = np.linspace( rho_v, rho_l, Np + 1 )
    dx = 2 * L / Np
    Nc = round( Np / 2 ) # index of the center
    rho[ Nc ] = rho_m
    for i in range( 0, Nc ):
        # Defining limits of integration left
        INT.x_i = rho_v
        INT.x_f = rho[ Nc - i ]
        # Integration left
        Int_l = INT.Simpson( PROF.ProfileFunc )
        # Defining limits of integration right
        INT.x_i = rho_v
        INT.x_f = rho[ Nc + i ]
        # Integral right
        Int_r = INT.Simpson( PROF.ProfileFunc )
        # Density computation
        rho[ Nc - i - 1 ] = rho[ Nc - i ] - dx * np.sqrt( np.abs( Int_l ) )
        rho[ Nc + i + 1 ] = rho[ Nc + i ] + dx * np.sqrt( np.abs( Int_r ) )
        if( rho[ Nc + i + 1 ] > rho_l ):
            rho[ Nc + i + 1 ] = rho_l
        if( rho[ Nc - i - 1 ] < rho_v ):
            rho[ Nc - i - 1 ] = rho_v
    return x, rho

# Density profile of order O(dx)^4
def ProfileO4( EOS, PROF, INT, Np, L, rho_v, rho_l ):
    # Computations
    rho_m = ( rho_v + rho_l ) / 2
    f = 1 / 12 # factor of the error
    x = np.linspace( -L, L, Np + 1 )
    rho = np.linspace( rho_v, rho_l, Np + 1 )
    dx = 2 * L / Np
    Nc = round( Np / 2 ) # index of the center
    rho[ Nc ] = rho_m
    for i in range( 0, Nc ):
        # Left part 
        # Limits of integration
        INT.x_i = rho_v
        INT.x_f = rho[ Nc - i ]
        # Integration
        Int_l = INT.Simpson( PROF.ProfileFunc )
        # Error terms
        mu_l = EOS.Chemical( rho[ Nc - i ], PROF.T )
        dmu_l = EOS.dChemical( rho[ Nc - i ], PROF.T )
        aux_1 = ( f / PROF.kappa ) * ( mu_l ** 2 - PROF.mu_0 ** 2 )
        aux_2 = ( f / PROF.kappa ) * 2 * PROF.mu_0 * ( mu_l - PROF.mu_0 )
        # Total term for density variation
        term_l = 1 / ( 1 + f * dmu_l ) ** 2 * ( Int_l + aux_1 - aux_2 )
        # Density computation
        rho[ Nc - i - 1 ] = rho[ Nc - i ] - dx * np.sqrt( np.abs( term_l ) )
        # Right part 
        # Limits of integration
        INT.x_i = rho_v
        INT.x_f = rho[ Nc + i ]
        # Integration
        Int_r = INT.Simpson( PROF.ProfileFunc )
        # Error terms
        mu_r = EOS.Chemical( rho[ Nc + i ], PROF.T )
        dmu_r = EOS.dChemical( rho[ Nc + i ], PROF.T )
        aux_1 = ( f / PROF.kappa ) * ( mu_r ** 2 - PROF.mu_0 ** 2 )
        aux_2 = ( f / PROF.kappa ) * 2 * PROF.mu_0 * ( mu_r - PROF.mu_0 )
        # Total term for density variation
        term_r = 1 / ( 1 + f * dmu_r ) ** 2 * ( Int_r + aux_1 - aux_2 )
        rho[ Nc + i + 1 ] = rho[ Nc + i ] + dx * np.sqrt( np.abs( term_r ) )
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

#*****************************************************************************
# Code Start

# Initialization

# Creating an object EOS to store EOS parameters
EOS = CS_EOS()
# Equation of State Parameters
EOS.a = 0.25
EOS.b = 4
EOS.R = 1
# Computing Reduced Properties
T_c = EOS.a / EOS.b / EOS.R / 2.6502
P_c = 0.18727 ** 2 / 0.4963 * EOS.a / EOS.b ** 2
rho_c = 0.5218 / EOS.b
EOS.T_c = T_c
EOS.P_c = P_c
EOS.rho_c = rho_c

# Defining Function Range
rho_v = 0.0006268162514430721
rho_l = 0.45406865373315114
P_0 = 7.311478997389506e-06

# Function properties
PROF = ProfileClass()
PROF.T = 0.50 * T_c
PROF.mu_0 = EOS.Chemical( rho_v, PROF.T )
PROF.kappa = 0.5;
PROF.pEOS = EOS

rho = np.linspace( rho_v, rho_l, 100 )
mu = EOS.Chemical( rho, PROF.T )

plt.plot( rho, mu )
plt.ylabel('chemical potential')
plt.xlabel('density')
# specifying horizontal line type
plt.axhline(y = PROF.mu_0, color = 'r', linestyle = '-')
# mark in liquid density
plt.plot( rho_v, PROF.mu_0, marker='o', color = 'k' )
plt.plot( rho_l, PROF.mu_0, marker='o', color = 'k' )
plt.show()

INT = Integration()
INT.x_i = rho_v
INT.x_f = rho_l
INT.N = 10000
Int_SMP = INT.Simpson( PROF.ProfileFunc )
print("Integral: ", Int_SMP)

# Interface caracteristics
L_domain = 16 # size of domain in lattice units
N_domain = 400 # number of points of the domain

x, rho_O2 = ProfileO2( EOS, PROF, INT, N_domain, L_domain, rho_v, rho_l )
w_O2, xO2_i, xO2_f = WIDTH( x, rho_O2, rho_v, rho_l, Np = N_domain )
print( "w O2: ", w_O2 )

x, rho_O4 = ProfileO4( EOS, PROF, INT, N_domain, L_domain, rho_v, rho_l )
w_O4, xO4_i, xO4_f = WIDTH( x, rho_O4, rho_v, rho_l, Np = N_domain )
print( "w O4: ", w_O4 )

plt.plot( x, rho_O2, color = 'r', linestyle = '-' )
plt.plot( x, rho_O4, color = 'k', linestyle = '-' )
plt.plot( xO2_i, rho_v + 0.12*( rho_l - rho_v ), marker='+', color = 'k', markersize=20  )
plt.plot( xO2_f, rho_v + 0.88*( rho_l - rho_v ), marker='+', color = 'k', markersize=20  )
plt.plot( xO4_i, rho_v + 0.12*( rho_l - rho_v ), marker='x', color = 'r', markersize=20  )
plt.plot( xO4_f, rho_v + 0.88*( rho_l - rho_v ), marker='x', color = 'r', markersize=20  )
plt.ylabel('density')
plt.xlabel('position')
plt.show()

# Save the results
strx = 'x'
stry = 'rho_O2'
strrho = 'rho_O4'
"""
#open file
file = open("Profile_FElbm.txt", "w")

file.write(strx + ";" + " " + stry + ";" + " " + strrho + ";" + "\n")

for i in range( 0, len(x) ):
     
    #convert variable to string
    str1 = repr( x[i] )
    str2 = repr( rho_O2[i] )
    str3 = repr( rho_O4[i] )
    file.write(str1 + ";" + " " + str2 + ";" + " " + str3 + ";" + "\n")
     
#close file
file.close()
""" 