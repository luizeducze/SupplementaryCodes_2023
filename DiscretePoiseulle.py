#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Luiz Eduardo Czelusniak

Compute the two-phase flow between parallel plates solution considering a 
discrete interface

Enter with the force Fx in lattice units

Enter the tau parameter to define viscosity

Enter the phase densities rho_1 and rho_2

Enter the half of the domain size R = H/2
"""

import numpy as np
import matplotlib.pyplot as plt

# Discrete solution of multiphase poiseulle

Fx = 2e-7

tau = 1
nu = ( tau - 0.5 ) / 3

rho_1 = 0.0030824221257919567
rho_2 = 0.40619262842616505

mu1 = nu * rho_1
mu2 = nu * rho_2

R = 50

du1_dr_0 = ( mu2 - mu1 )/mu1/( mu1 + mu2 ) * Fx * R / 2 
u1 = Fx*R*R/2/mu1 - R*du1_dr_0 

N = 100

r = np.linspace(0,R,N)

u1r = u1 - Fx/mu1*r*r/2 + du1_dr_0*r
u2r = u1 - Fx/mu2*r*r/2 - mu1/mu2*du1_dr_0*r

plt.plot(-r,u1r,"r")
plt.plot(r,u2r,"k")
plt.show

# Save the results
strx = 'y'
stry = 'U1'
strrho = 'U2'

# Uncomment this last part to save the results

"""
#open file
file = open("Discrete_Ux_Tr060.txt", "w")

file.write(strx + ";" + " " + stry + ";" + " " + strrho + ";" + "\n")

for i in range( 0, N ):
     
    #convert variable to string
    str1 = repr( r[i] )
    str2 = repr( u1r[i] )
    str3 = repr( u2r[i] )
    file.write( str1 + ";" + " " + str2 + ";" + " " + str3 + ";" + "\n")
     
#close file
file.close()
"""