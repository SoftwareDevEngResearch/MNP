# -*- coding: utf-8 -*-
"""
Created on Tue Mar 15 00:11:24 2022

@author: ethan
"""

import handler as h
import numpy as np

material_name_string = 'sample/XS_G12.npz'

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
decay       = data["decay"] #lambda
beta_frac   = data["beta_frac"]

input_name_string = 'sample/input.IN'
    
with open(input_name_string,"r") as f:
    
    T0      = float(h.stringSplitter(f.readline()))   #initial time
    T       = float(h.stringSplitter(f.readline()))   #final time
    NT      = int(h.stringSplitter(f.readline()))     #Number of time steps
    method  = int(h.stringSplitter(f.readline()))     #timestepping method to use
    X       = float(h.stringSplitter(f.readline()))   #radius of object
    Nx      = int(h.stringSplitter(f.readline()))     #number of space grid points
    BCL     = float(h.stringSplitter(f.readline()))   #left reflectance boundary condition
    BCR     = float(h.stringSplitter(f.readline()))   #right reflectance boundary condition
    geo     = int(h.stringSplitter(f.readline()))     #geometry 1-slab, 2-cyl, 3-sphere
    verbose = bool(h.stringSplitter(f.readline()))    #whether or not to store output data for every step
    timeit  = bool(h.stringSplitter(f.readline()))    #whether or not to benchmark the solution
    saveit  = bool(h.stringSplitter(f.readline()))    #whether or not to save the output to csv
    alpha   = [BCL,BCR]
    
matLib = np.array((E,SigmaT,SigmaF,SigmaS,v,nu_prompt,nu_delayed,chi_prompt,chi_delayed,decay,beta_frac),dtype = object)

def test_grid_creation():
    
    grid = h.Grid(X,Nx,matLib,geo,T0,T,NT,0,alpha)
    
    assert hasattr(grid,'geo')
    
def test_grid_material_property_assignment():
    
    grid = h.Grid(X,Nx,matLib,geo,T0,T,NT,0,alpha)
    grid.neutronicsEnergyGridProperties(matLib)
    
    assert hasattr(grid,'E') and hasattr(grid,'G') and hasattr(grid,'v')
    
def test_grid_delayed_property_assignment():
    
    grid = h.Grid(X,Nx,matLib,geo,T0,T,NT,0,alpha)
    grid.delayedEnergyGridProperties(matLib)
    
    assert hasattr(grid,'J') and hasattr(grid,'decay')
    
def test_grid_build_cells():
    
    grid = h.Grid(X,Nx,matLib,geo,T0,T,NT,0,alpha)
    grid.buildCells()
    
    assert grid.cells_built == 1
    
def test_grid_build_cell_matrices():
    
    grid = h.Grid(X,Nx,matLib,geo,T0,T,NT,0,alpha)
    grid.neutronicsEnergyGridProperties(matLib)
    grid.delayedEnergyGridProperties(matLib)
    grid.buildCells()
    
    for i in range(grid.NumX):
        grid.cellsList[i].neutronicsMaterialProperties(matLib)
        grid.cellsList[i].delayedMaterialProperties(matLib)
        grid.cellsList[i].buildRow()
    
    for i in range(grid.NumX):
        assert hasattr(grid.cellsList[i],'A_contribution')
        
def test_grid_step():
    
    grid = h.Grid(X,Nx,matLib,geo,T0,T,NT,0,alpha)
    grid.neutronicsEnergyGridProperties(matLib)
    grid.delayedEnergyGridProperties(matLib)
    grid.step()
    
    assert grid.current_time_step != 0
    
def test_output_not_verbose():
    grid = h.Grid(X,Nx,matLib,geo,T0,T,NT,0,alpha)
    grid.neutronicsEnergyGridProperties(matLib)
    grid.delayedEnergyGridProperties(matLib)
    grid.step()
    grid.Run()
    grid.buildPHI_Total()
    phiFinal = grid.PHI
    np.savez('./Output/phi_final.npz',phi_final = phiFinal)
    
    data = np.load('Output/phi_final.npz',allow_pickle = True)
    assert len(data['phi_final']) == grid.G*grid.NumX
    
def test_output_verbose():
    grid = h.Grid(X,Nx,matLib,geo,T0,T,NT,0,alpha)
    grid.neutronicsEnergyGridProperties(matLib)
    grid.delayedEnergyGridProperties(matLib)
    grid.step()
    grid.Run_verbose()
    grid.buildPHI_Total()
    phiHist = grid.PHI_History
    np.savez('./Output/phi_Hist.npz',t = grid.t,phi_hist = phiHist)
    
    data = np.load('Output/phi_hist.npz',allow_pickle = True)
    assert len(data['phi_hist']) == grid.G*grid.NumX and len(data['t']) == grid.numTimePoints