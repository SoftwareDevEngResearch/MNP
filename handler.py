# -*- coding: utf-8 -*-
"""
Created on Tue Jan 25 13:42:43 2022

@author: ethan
"""

import numpy as np
import DelayedFunctions as df #will clean up later

# The purpose of this file is to handle calling the functions as necessary to step
# individual cells forward in time.

# For now, each cell is considered as an infinite medium, which has a loss term
# to its neighbors and a source term from its neighbors

#read in input data from XS_G12.npz as an example
data        = np.load('XS_G12.npz')
E           = data['E']
SigmaT      = data["SigmaT"]
SigmaF      = data["SigmaF"]
SigmaS      = data["SigmaS"]
v           = data["v"] * 1000
nu_prompt   = data["nu_prompt"]
nu_delayed  = data["nu_delayed"]
chi_prompt  = data["chi_prompt"]
chi_delayed = data["chi_delayed"]
decay       = data["decay"] #lambda
beta_frac   = data["beta_frac"]

#read program input file to decide what controls should be on the simulation.
with open('input.IN',"r") as f:
    
    T   = f.readline() #final time, initial is assumed to be zero
    dt  = f.readline() #timestep size
    X   = f.readline() #radius of object
    Nx  = f.readline() #number of space grid points
    BCL = f.readline() #left reflectance boundary condition
    BCR = f.readline() #right reflectance boundary condition
    geo = f.readline() #geometry 1-slab, 2-cyl, 3-sphere
    assert geo == 0 or geo == 1 or geo == 2
    
    pass


todo = 'boggle'

class Grid:
    
    def __init__(self,idx):
        local_r_pos = df.create_grid(X,Nx)[1][idx]
        left_leakage_term = todo
        right_leakage_term = todo
        
    pass