#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 18 13:32:39 2023

@author: Luiz Eduardo Czelusniak
"""

import numpy as np

# This class provide all the necessary EOS information
class CS_EOS:
    # EOS parameters a, b and R
    a = 1
    b = 1
    R = 1
    # EOS critical properties
    T_c = 1
    P_c = 1
    rho_c = 1
        
    # Compute EOS pressure from density and temperature
    def P_from_rho( self, rho, T ):
        # Defining EOS varaibles
        a = self.a
        b = self.b
        R = self.R  
        # Auxiliary parameters for EOS
        n = b * rho / 4
        num = ( 1 + n + n ** 2 - n ** 3 )
        den = ( 1 - n ) ** 3
        aux = num / den
        # Compute EOS
        Peos = rho * R * T * aux - a * rho ** 2
        return Peos
    
    # Compute the EOS pressure derivative in respect to rho
    def dP_from_rho( self, rho, T ):
        # Defining EOS variables
        a = self.a
        b = self.b
        R = self.R
        # Auxiliary parameters for dP
        n = b * rho / 4
        T1 = rho * R * T
        T2 = 1 + n + n ** 2 - n ** 3
        T3 = ( 1 - n ) ** 3
        dT1 = R * T
        dT2 = b / 4 + 2 * n * b / 4 - 3 * n * n * b / 4
        r_dT3_T3 = 3 * ( - b / 4 ) / ( 1 - n ) ** 4
        # Compute EOS derivative in respect to rho
        dP = dT1 * T2 / T3 + T1 * dT2 / T3 - T1 * T2 * r_dT3_T3 - 2 * a * rho
        return dP
    
    # Compute EOS vapor density from pressure and temperature
    def rhov_from_P( self, P_0, T ):
        # Defining the range of the search
        rho_right = self.MAX_EOS( T )
        rho_left = 0
        # Bissection method
        tol = 1e-10
        while ( rho_right - rho_left ) > tol:
            rho_med = ( rho_right + rho_left ) / 2
            P_med = self.P_from_rho( rho_med, T )
            if P_med > P_0:
                rho_right = rho_med
            else:
                rho_left = rho_med
        return rho_med
    
    # Compute EOS liquid density from pressure and temperature
    def rhol_from_P( self, P_0, T ):
        # Defining the range of the search
        rho_right = 3 / self.b
        if self.P_from_rho( rho_right, T ) < P_0:
            print("Problems in liquid density, adjust right density")
        rho_left = self.MIN_EOS( T )
        # Bissection method
        tol = 1e-8
        while ( rho_right - rho_left ) > tol:
            rho_med = ( rho_right + rho_left ) / 2
            P_med = self.P_from_rho( rho_med, T )
            if P_med > P_0:
                rho_right = rho_med
            else:
                rho_left = rho_med
        return rho_med
    
    # Find the point of Maximum pressure of the EOS for a certain T
    def MAX_EOS( self, T ):
        # Defining necessary parameters
        b = self.b
        rho_c = self.rho_c
        # Interval of search
        rho_i = 0
        rho_f = 3 / b
        # Incremental Search
        drho = rho_c / 100
        rho = rho_i
        verifier = 0
        while verifier == 0:
            dP = self.dP_from_rho( rho, T  )
            rho_next = rho + drho
            dP_next = self.dP_from_rho( rho_next, T  )
            if dP > 0 and dP_next < 0:
                rho_min = rho
                rho_max = rho_next
                verifier = verifier + 1
            rho = rho_next
            if rho > rho_f:
                print("Problems to find EOS maximum")
                break
        # Bissection Method
        tol = 1e-6
        while ( rho_max - rho_min ) > tol:
            rho_med = ( rho_max + rho_min ) / 2
            dP_med = self.dP_from_rho( rho_med, T  )
            if dP_med > 0:
                rho_min = rho_med
            else:
                rho_max = rho_med
        return rho_med
    
    # Find the point of Minimum pressure of the EOS for a certain T
    def MIN_EOS( self, T ):
        # Defining necessary parameters
        b = self.b
        rho_c = self.rho_c
        # Interval of search
        rho_i = 0
        rho_f = 3 / b
        # Incremental Search
        drho = rho_c / 100
        rho = rho_i
        verifier = 0
        while verifier == 0:
            dP = self.dP_from_rho( rho, T  )
            rho_next = rho + drho
            dP_next = self.dP_from_rho( rho_next, T  )
            if dP < 0 and dP_next > 0:
                rho_min = rho
                rho_max = rho_next
                verifier = verifier + 1
            rho = rho_next
            if rho > rho_f:
                print("Problems to find EOS minimum")
                break
        # Bissection Method
        tol = 1e-6
        while ( rho_max - rho_min ) > tol:
            rho_med = ( rho_max + rho_min ) / 2
            dP_med = self.dP_from_rho( rho_med, T  )
            if dP_med < 0:
                rho_min = rho_med
            else:
                rho_max = rho_med
        return rho_med
    
    # Function to compute the chemical potential
    def Chemical( self, rho, T ):
        # Defining EOS varaibles
        a = self.a
        b = self.b
        R = self.R
        # Auxiliary parameters for chemical potential
        n = b * rho / 4
        num = ( 1 + n + n ** 2 - n ** 3 )
        den = ( 1 - n )
        aux_1 = R * T * num / den ** 3
        aux_2 = 2 * R * T / den
        aux_3 = R * T / den ** 2
        aux_4 = - 2 * a * rho
        aux_5 = R * T * np.log( rho )
        # chemical potential
        mu = aux_1 + aux_2 + aux_3 + aux_4 + aux_5
        return mu
    
    # Function to compute the derivative of chemical potential
    def dChemical( self, rho, T ):
        # Defining EOS varaibles
        a = self.a
        b = self.b
        R = self.R
        # Auxiliary parameters for chemical potential derivative
        n = b * rho / 4
        num = ( 1 + n + n ** 2 - n ** 3 )
        den = ( 1 - n )
        aux_1 = R * T * b / 4 * ( 1 + 2 * n - 3 * n ** 2 ) / den ** 3
        aux_2 = R * T * 3 * b / 4 * num / den ** 4
        aux_3 = R * T * b / 2 / den ** 2
        aux_4 = R * T * b / 2 / den ** 3
        aux_5 = - 2 * a
        aux_6 = R * T / rho
        # Derivative of chemical potential
        dmu = aux_1 + aux_2 + aux_3 + aux_4 + aux_5 + aux_6
        return dmu
    
#*****************************************************************************

class Integration:
    
    # Range of integration
    x_i = 0 
    x_f = 1 
    # Number of nodes for integration
    N = 10
    
    # Compute the integral with rectangular method
    def Rectangular( self, func ):
        # Importing variables
        x_i = self.x_i
        x_f = self.x_f
        N = self.N
        # Preparing integration
        x = x_i
        dx = ( x_f - x_i ) / N
        Int = 0
        # Computing integral
        for i in range( 0, N ):
            f = func( x )
            Int = Int + f * dx
            x = x + dx
        return Int

    # Compute the integral with Simpson's rule
    def Simpson( self, func ):
        # Importing variables
        x_i = self.x_i
        x_f = self.x_f
        N = self.N
        # Preparing integration
        x = x_i
        dx = ( x_f - x_i ) / N
        Int = 0
        Nd = round( N / 2 ) # Number of nodes after division
        # Computing integral
        for i in range( 0, Nd ):
            f_0 = func( x )
            f_1 = func( x + dx )
            f_2 = func( x + 2 * dx )
            Int = Int + ( f_0 + 4 * f_1 + f_2  ) / 3 * dx
            x = x + 2 * dx
        return Int
    
#*****************************************************************************    
    
class Interpolation:
    
    # Class initialization
    def __init__( self, x_input, y_input ):
        self.xv = x_input
        self.yv = y_input
        self.Nv = len( x_input )
        
    # Compute interpolation
    def Interpol( self, x ):
        # Verify if x is a valid point
        if( x < self.xv[0] - 1e-12 or x > self.xv[self.Nv-1] + 1e12 ):
            print( "Incorrect x value for interpolation" )
        
        # Loop for all xv values
        dx = ( self.xv[self.Nv-1] - self.xv[0] ) / ( self.Nv - 1 )
        j = np.floor( ( x - self.xv[0] ) / dx )
        if( j > 1 ):
            j = j - 1
        for i in range( round(j), self.Nv - 1 ):
            if( x >= self.xv[i] - 1e-12 and x <= self.xv[i+1] + 1e12  ):
                dy = self.yv[i+1] - self.yv[i]
                dx = self.xv[i+1] - self.xv[i]
                y = self.yv[i] + ( x - self.xv[i] ) * dy / dx
        
        return y
    

            
    
    
    
    
    
    
    
    
    
    
    

    
    