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
# probably this should be done with a command line input pointed at a file with 
# the same structure

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
    
    T   = float(f.readline()) #final time, initial is assumed to be zero
    dt  = float(f.readline()) #timestep size
    X   = float(f.readline()) #radius of object
    Nx  = int(f.readline()) #number of space grid points
    BCL = float(f.readline()) #left reflectance boundary condition
    BCR = float(f.readline()) #right reflectance boundary condition
    geo = int(f.readline()) #geometry 1-slab, 2-cyl, 3-sphere
    assert geo == 0 or geo == 1 or geo == 2
    
    pass

################# placeholders #########################

todo = 'boggle'

matLib = np.array(data,dtype = object)

########################################################
class Grid:
    
    def __init__(self,idx,matNo):

        local_r_pos = df.create_grid(X,Nx)[1][idx]
        local_time  = 0
        #TODO: Make material library
        SigmaT      = matLib[matNo]["SigmaT"] #total cross section
        SigmaF      = matLib[matNo]["SigmaF"] #fission cross section
        SigmaS      = matLib[matNo]["SigmaS"] #down scattering cross section
        nu_prompt   = matLib[matNo]["nu_prompt"] #probability a neutron is born fast
        nu_delayed  = matLib[matNo]["nu_prompt"]
        
    def queryNeighbors():

        left_outgoing_term   = todo
        right_outgoing_term  = todo
        left_incoming_term  = todo
        right_incoming_term = todo
        
    pass

test = Grid(0,0)