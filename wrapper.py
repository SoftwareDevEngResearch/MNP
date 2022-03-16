# -*- coding: utf-8 -*-
"""
Created on Thu Mar 10 17:15:16 2022

@author: ethan
"""

import numpy as np
import handler as h
import os
import sys

if len(sys.argv) == 1:
    print('run with $ wrapper.py material.npz input.IN')
    print('see input.in for more details')
    print('open material.npz for the proper ordering')
    print('')
    sys.exit()

if sys.argv[1] == '-h':
    print('run with $ wrapper.py material.npz input.IN')
    print('see input.in for more details')
    print('open material.npz for the proper ordering')
    sys.exit()
    
if sys.argv[1] == 'go left':
    print('eaten by a grue!')
    sys.exit()
    
assert len(sys.argv) == 3

material_name_string    = sys.argv[1]  
input_name_string       = sys.argv[2] 

data        = np.load(material_name_string)
E           = data['E']
SigmaT      = data["SigmaT"]
SigmaF      = data["SigmaF"]
SigmaS      = data["SigmaS"]
v           = data["v"] * 1000
nu_prompt   = data["nu_prompt"]
nu_delayed  = data["nu_delayed"]
chi_prompt  = data["chi_prompt"]
chi_delayed = data["chi_delayed"]
decay       = data["decay"]
beta_frac   = data["beta_frac"]

def stringSplitter(string):
    idx = string.find('#')
    return string[:idx]

with open(input_name_string,"r") as f:
    T0      = float(stringSplitter(f.readline()))   #initial time
    T       = float(stringSplitter(f.readline()))   #final time
    NT      = int(stringSplitter(f.readline()))     #Number of time steps
    method  = int(stringSplitter(f.readline()))     #timestepping method to use
    X       = float(stringSplitter(f.readline()))   #radius of object
    Nx      = int(stringSplitter(f.readline()))     #number of space grid points
    BCL     = float(stringSplitter(f.readline()))   #left reflectance boundary condition
    BCR     = float(stringSplitter(f.readline()))   #right reflectance boundary condition
    geo     = int(stringSplitter(f.readline()))     #geometry 1-slab, 2-cyl, 3-sphere
    verbose = bool(stringSplitter(f.readline()))    #whether or not to store output data for every step
    timeit  = bool(stringSplitter(f.readline()))    #whether or not to benchmark the solution
    saveit  = bool(stringSplitter(f.readline()))    #whether or not to save the output to csv
    alpha   = [BCL,BCR]
assert geo == 0 or geo == 1 or geo == 2

if timeit == True:
    import time
    t_start = time.time()

matLib = np.array((E,SigmaT,SigmaF,SigmaS,v,nu_prompt,nu_delayed,chi_prompt,chi_delayed,decay,beta_frac),dtype = object)

grid = h.Grid(X,Nx,matLib,geo,T0,T,NT,method,alpha)
grid.neutronicsEnergyGridProperties(matLib)
grid.delayedEnergyGridProperties(matLib)
grid.buildCells()
for i in range(Nx):
        grid.cellsList[i].neutronicsMaterialProperties(matLib)
        grid.cellsList[i].delayedMaterialProperties(matLib)
grid.cellsList[0].initializeFlux(v)
if verbose == True:
    grid.Run_verbose()
    phiHist = grid.PHI_History
    np.savez('./Output/phi_Hist.npz',t = grid.t,phi_hist = phiHist)
elif verbose == False:
    grid.Run()
    grid.buildPHI_Total()
    phiFinal = grid.PHI()
    np.savez('./Output/phi_final.npz',phi_final = phiFinal)

if timeit == True:
    print('Runtime: {:.5} sec.'.format(time.time() - t_start))
