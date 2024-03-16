#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Luiz Eduardo Czelusniak
"""

import matplotlib.pyplot as plt
import numpy as np
from Miscellaneous import CS_EOS, Integration, Interpolation

plt.close('all')

#*****************************************************************************
# Classes

# This class will be used to compute the density profile
class ProfileClass:
    # Variables
    T = 1 # Temperature
    P_0 = 1 # equilibrium pressure
    eps = 1 # term to control densities
    B = 1 # term to control interface width
    
    # Pointer to the desired EOS object
    pEOS = 0
    
    def PSI( self, rho ):
        G = -1
        c = 1
        cs2 = 1 / 3
        P_eos = self.pEOS.P_from_rho( rho, self.T )
        psi = np.sqrt( 2 * ( P_eos - rho * cs2 ) / G / c / c )
        return psi
    
    # Derivative of pseudopotential function
    def dPSI( self, rho ):
        G = -1
        c = 1
        cs2 = 1 / 3
        psi = self.PSI( rho )
        dP_eos = self.pEOS.dP_from_rho( rho, self.T )
        dpsi = 1 / psi / G / c ** 2 * ( dP_eos - cs2 )
        return dpsi 
    
    def ProfileFunc( self, rho ):
        # P from equation of state
        P_eos = self.pEOS.P_from_rho( rho, self.T )
        # pseudopotential function
        psi = self.PSI( rho )
        # profile function
        f = 2 * ( self.P_0 - P_eos ) / self.B / psi ** ( 1 + self.eps )
        f = f * self.dPSI( rho )
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
    
    # Computing integral from rho_v to rho_m 
    INT.x_i = rho_v
    INT.x_f = rho_m 
    INT.N = 20000
    Int_m = INT.Simpson( PROF.ProfileFunc )
    Int_l = Int_m
    Int_r = Int_m
    
    INT.N = 1000
    
    for i in range( 0, Nc ):
        if( rho[ Nc - i ] > rho_v ):
            # Defining limits of integration left
            if( i == 0 ):
                INT.x_f = rho_m
            else:
                INT.x_f = rho[ Nc - i + 1 ]
            #INT.x_i = rho_v
            INT.x_i = rho[ Nc - i ]
            # Integration left
            Int_l = Int_l - INT.Simpson( PROF.ProfileFunc )
            F_l = PROF.PSI( rho[ Nc - i ] ) ** PROF.eps * Int_l
            # Density computation
            rho[ Nc - i - 1 ] = rho[ Nc - i ] - dx * np.sqrt( np.abs( F_l ) )
            if( rho[ Nc - i - 1 ] < rho_v ):
                rho[ Nc - i - 1 ] = rho_v
        else:
            rho[ Nc - i - 1 ] = rho_v
            
        if( rho[ Nc + i ] < rho_l ):
            # Defining limits of integration right
            if( i == 0 ):
                INT.x_i = rho_m
            else:
                INT.x_i = rho[ Nc + i - 1 ]
            #INT.x_i = rho_v
            INT.x_f = rho[ Nc + i ]
            # Integral right
            Int_r = Int_r + INT.Simpson( PROF.ProfileFunc )
            F_r = PROF.PSI( rho[ Nc + i ] ) ** PROF.eps * Int_r
            # Density computation
            rho[ Nc + i + 1 ] = rho[ Nc + i ] + dx * np.sqrt( np.abs( F_r ) )
            if( rho[ Nc + i + 1 ] > rho_l ):
                rho[ Nc + i + 1 ] = rho_l
        else:
            rho[ Nc + i + 1 ] = rho_l
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
EOS.a = 0.25 # <- Please, add here 'a' for CS-EOS
EOS.b = 4 # <- Please, add here 'b' for CS-EOS
EOS.R = 1 # <- Please, add here 'R' for CS-EOS
# Computing Reduced Properties
T_c = EOS.a / EOS.b / EOS.R / 2.6502
P_c = 0.18727 ** 2 / 0.4963 * EOS.a / EOS.b ** 2
rho_c = 0.5218 / EOS.b
EOS.T_c = T_c
EOS.P_c = P_c
EOS.rho_c = rho_c

# Defining Function Range
"""
To compute the density profile is computed from an integration process,
so we need to specify the range for integration from the equilibrium vapor
density to the equilibrium liquid density given by the Maxwell rule
for the desired temperature;
we also need to specify the equilibrium pressure;
these information is obtained in the code Maxwell_Densities.py
"""
rho_v = 0.0030824221257919567 # <- Please, add here equilibrium vapor density
rho_l = 0.40619262842616505 # <- Please, add here equilibrium liquid density
P_0 = 4.178249509272372e-05 # <- Please, add here equilibrium pressure

# Function properties
"""
We create a Profile class which is store the function used to generate the
density profile; the profile is generated by the integration of this function
"""
PROF = ProfileClass()
PROF.T = 0.60 * T_c # <- Please, add here the desired temperature
PROF.P_0 = P_0 
PROF.eps = 1.9373335029111907 # <- Please, add here the desired kappa value
kappa = 0.9190714865493067
PROF.B = - kappa / 4
PROF.pEOS = EOS

# Plot Chemical Potential Function
# Defining vectors for plotting
rho = np.linspace( rho_v, rho_l, 100 )
Peos = EOS.P_from_rho( rho, PROF.T )
# Plotting
plt.plot( rho, Peos )
plt.ylabel('Pressure')
plt.xlabel('density')
# specifying horizontal line for equilibrium pressure
plt.axhline(y = PROF.P_0, color = 'r', linestyle = '-')
# mark vapor and liquid equilibrium densities
plt.plot( rho_v, PROF.P_0, marker='o', color = 'k' )
plt.plot( rho_l, PROF.P_0, marker='o', color = 'k' )
plt.show()

# Defining Class Integration
"""
We defined previously the function used to compute the density profile;
we need to integrate this functions, integrations can be done by the 
integration class
"""
INT = Integration()
INT.x_i = rho_v
INT.x_f = rho_l
INT.N = 10000 # <- Please add here the number of nodes desired for integration
Int_SMP = INT.Simpson( PROF.ProfileFunc )
print("Integral: ", Int_SMP)

# Computing the density profile
Nint = 1600 # <- Please, add here the number of nodes for density profile
"""
H is equal to the half of the domain physical size; 
if domain has total high of Ly, then H = Ly/2;
remember, if we want to compute a solution to compare with LBM 
than H must represent a size in lattice units;
"""
H = 50 # <- Please, add here the half size of the domain
"""
The next function compute density profile; the next one compute
the interface thickness
"""
x, rho_O2 = ProfileO2( EOS, PROF, INT, Nint, H, rho_v, rho_l )
w_O2, xO2_i, xO2_f = WIDTH( x, rho_O2, rho_v, rho_l, Np = Nint )
print( "w O2: ", w_O2 )

# Plot the density profile
plt.plot( x, rho_O2, color = 'r', linestyle = '-' )
plt.plot( xO2_i, rho_v + 0.12*( rho_l - rho_v ), marker='+', color = 'k', markersize=20  )
plt.plot( xO2_f, rho_v + 0.88*( rho_l - rho_v ), marker='+', color = 'k', markersize=20  )
plt.ylabel('density')
plt.xlabel('position')
plt.show()

"""
Now, we compute the theoretical velocity profile. 
The equation for the continuous profile is:
    Ux - Ux(0) = - Fx int_{0}^{y} y/mu(y) dy 
                 + mu(0) d/dy(Ux(0)) int_{0}^{y} 1/mu(y) dy
Observe that our solution depends on the two integrals:
    1) int_{0}^{y} y/mu(y) dy 
    2) int_{0}^{y} 1/mu(y) dy
"""

# Computing the velocity profile
"""
We start by defining Fx, and tau (used for nu);
note that we want to compute the integral related with mu,
so we need nu to compute mu, and then integrate
"""
Fx = 2e-7
tau = 1
nu = ( tau - 0.5 ) / 3
"""
Now we create a velocity vector to store our solution,
the y vector is created from x, remember x places the interface
at x = 0, so we add H, to make y vary from zero to 2H
"""
Ux = np.zeros( Nint + 1 )
y = x + H
dy = x[1] - x[0] # Here we compute de space step dy

# Computing dynamic viscosity
mu = nu * rho_O2

"""
We want to compute the integrals related with y and mu;
However, these variables are vectors and our Integral class compute
integrals from functions and not arrays. Then, we first create an interpolation
class;
This is a function that compute the variable at any point by interpolating the
vector;
"""
# Creating class interpolation
f1 = y / mu
f2 = 1 / mu
INTP1 = Interpolation( y, f1 )
INTP2 = Interpolation( y, f2 )

# Plot functions y/mu and 1/mu for physical interpretation
plt.plot( y, y/mu, "r" )
plt.plot( y, 1/mu, "k" )
plt.show()

# Creating class integral
"""
We had the vector y and mu; then we created the interpolation function;
now we pass the interpolation to the class integration;
first we create vectors to store the integral for all nodes of our density
profile;
"""
Int1 = np.zeros( Nint + 1 )
Int2 = np.zeros( Nint + 1 )

# 1: This class integrate y/mu
INT1 = Integration()
INT1.x_i = y[0]
INT1.x_f = y[Nint]
INT1.N = 500

# 2: This class integrate 1/mu
INT2 = Integration()
INT2.x_i = y[0]
INT2.x_f = y[Nint]
INT2.N = 500

# Next we compute the entire value of the integral
"""
One of the variables of the analytical solution is not know; 
which is the term d/dy(Ux(0)); then we integrate until y=H which result in:
    Fx int_{0}^{H} y/mu(y) dy + mu(0) d/dy(Ux(0)) int_{0}^{H} 1/mu(y) dy
"""
Int_f1 = INT1.Simpson( INTP1.Interpol )
Int_f2 = INT2.Simpson( INTP2.Interpol )
# Then we are able to compute the unknown term
mu0_x_dUx_dy_0 = Fx * Int_f1 / Int_f2 

"""
Lets separate the integral in three parts:
    vapor region;
    interface;
    liquid region;
The integral just need more nodes in the interface
"""

INT1.N = 20
INT2.N = 20

# Loop
# Now we compute the integral for each term  
for i in range( 1, Nint + 1 ):
    # 1
    INT1.x_i = y[i-1]
    INT1.x_f = y[i]
    Int1[i] = Int1[i-1] + INT1.Simpson( INTP1.Interpol )
    # 2
    INT2.x_i = y[i-1]
    INT2.x_f = y[i]
    Int2[i] = Int2[i-1] + INT2.Simpson( INTP2.Interpol )
    #
    print(i)
  
Ux = - Fx * Int1 +  mu0_x_dUx_dy_0 * Int2


plt.plot( y, Ux )
# specifying horizontal line type
plt.axhline(y = 0, color = 'r', linestyle = '-')
# mark in liquid density
plt.plot( y[0], 0, marker='o', color = 'k' )
plt.plot( y[Nint-1], 0, marker='o', color = 'k' )
plt.show()

""" 
Next we also plot the discrete solution
"""

# Discrete solution of multiphase poiseulle
mu1 = mu[0]
mu2 = mu[Nint-1]

R = H

du1_dr_0 = ( mu2 - mu1 )/mu1/( mu1 + mu2 ) * Fx * R / 2 
u1 = Fx*R*R/2/mu1 - R*du1_dr_0 

r = np.linspace(0,R,40)

u1r = u1 - Fx/mu1*r*r/2 + du1_dr_0*r
u2r = u1 - Fx/mu2*r*r/2 - mu1/mu2*du1_dr_0*r

plt.plot(-r+R,u1r, color = 'r', linestyle = '--')
plt.plot(r+R,u2r, color = 'r', linestyle = '--')
plt.plot( y, Ux )
# specifying horizontal line type
plt.axhline(y = 0, color = 'r', linestyle = '-')
# mark in liquid density
plt.plot( y[0], 0, marker='o', color = 'k' )
plt.plot( y[Nint-1], 0, marker='o', color = 'k' )
plt.show()

# Save the results
strx = 'y'
stry = 'Ux'
strrho = 'rho'

#open file
file = open("Continuous_Ux_Plbm.txt", "w")

file.write(strx + ";" + " " + stry + ";" + " " + strrho + ";" + "\n")

for i in range( 0, Nint + 1 ):
     
    #convert variable to string
    str1 = repr( y[i] )
    str2 = repr( Ux[i] )
    str3 = repr( rho_O2[i] )
    file.write(str1 + ";" + " " + str2 + ";" + " " + str3 + ";" + "\n")
     
#close file
file.close()

