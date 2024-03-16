# SuppMaterial_Alfa01
This project contains the codes used in our publication comparing P-LBM and FE-LBM

Miscellaneous.py: 
    Contain classes and functions used by all other codes, such as
    EOS related functions, integration and interpolation functions
    
Maxwell_Densities.py:
    Given EOS parameters and reduced temperature, this code return
    the equilibrium densities given by the Maxwell rule
    
Droplet_Consistent_Densities.py:
    Given EOS parameters, reduced temperaute, droplet radius and surface tension
    this code return the thermodynamic consistent equilibrium densities
    
Density_Profile_FElbm.py:
    Given EOS parameters, reduced temperature and kappa parameter, this code
    return the density profile for the free-energy method
    
Surface_Tension_FElbm.py:
    Given EOS parameters, reduced temperature and kappa parameters, this code
    reutrn the surface tension for the free-energy method
    
Eps_Pseudopotential.py:
    Compute the theoretical epsilon parameter of pseudopotential method to 
    match the same planar interface equilibrium densities of Maxwell rule
    
Density_Profile_Plbm.py:
    Compute the theoretical density profile for a planar interface for the 
    pseudopotential method
    
SurfaceTension_Plbm.py:
    Compute the theoretical surface tension for a planar interface for the
    pseudopotential method
    
DiscretePoiseulle.py:
    Compute the theoretical velocity profile for a two-phase flow between 
    parallel plates considering a discrete interface
    
ContinuousPoiseulle_FElbm.py: 
    Compute the theoretical velocity profile for a two-phase flow between 
    parallel plates considering a diffuse interface with the free-energy method
    
ContinuousPoiseulle_Plbm.py:
    Compute the theoretical velocity profile for a two-phase flow between 
    parallel plates considering a diffuse interface with the pseudopotential 
    method
    

    
    
    
    
    
    