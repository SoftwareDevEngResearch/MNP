# -*- coding: utf-8 -*-
"""
Created on Tue Jan 25 13:42:43 2022

@author: ethan
"""

import numpy as np
import os
import csv
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

matLib = np.array((E,SigmaT,SigmaF,SigmaS,v,nu_prompt,nu_delayed,chi_prompt,chi_delayed,decay,beta_frac))
#                 [0,1     ,2     ,3     ,4,5,       ,6,        ,7         ,8          ,9    ,10]
########################################################
class Grid:
    #creates the computational grid, 1d for now.
    #   initializes a time to be zero and calls a function to define points on the time grid
    #   
    def __init__(self,maxX,NumX,aNeuMatLib,geo,startTime,finalTime,numTimePoints,method): 
        #create the grid
        
        #These are defined for the computataional spatial grid
        self.geo = geo
        self.dr, self.edges,self.centers = df.create_grid(X,Nx)
        self.current_time_step = 0
        #build surface area and volumes for the grid
        assert self.geo == 0 or self.geo == 1 or self.geo == 2
        #0: slab, 1: cylinder, 2: sphere
        self.S = np.zeros(Nx+1) #surface area on the edges
        self.V = np.zeros(Nx)   #volume of the cells
        if self.geo == 0:
            self.S[:] = 1
            self.V = self.dr
        elif self.geo == 1:
            self.S = 4*np.pi*self.edges**2
            self.V = 4/3*np.pi*(self.edges[1:]**3 - self.edges[:-1]**3)
        elif self.geo == 2:
            self.S = 2*np.pi*self.edges
            self.V = np.pi*(self.edges[1:]**2 - self.edges[:-1]**2)
        assert method == 1 or method == 0
        if method == 0:
            self.t = np.linspace(startTime,finalTime,numTimePoints)
        elif method == 1:
            self.t = np.logspace(startTime,finalTime,numTimePoints)
        
        self.global_time    = self.t[self.current_time_step]
        self.dts            = self.t[1:] - self.t[:-1]
        self.current_dt     = self.dts[self.current_time_step]
        pass
            
    def neutronicsEnergyGridProperties(self,aNeuMatLib):
        self.E      = aNeuMatLib[0] #property of the neutronics grid, discrete energy bounds
        self.v      = aNeuMatLib[4] #property of the neutronics grid, discrete speed corresponding to energy
        self.G      = len(self.E) - 1  #number of points in the energy grid
        
    def delayedEnergyGridProperties(self,aNeuMatLib):
        self.decay  = aNeuMatLib[7] #property of the delayed neutron grid, discrete decay constant vector
        self.J      = len(self.decay)  #number of points in the precursor flavor grid
        
class Cell:
    #creates an individual cell of the grid, needs to have a parent grid.
        
    def __init__(self,idx,Grid):
        
        self.idx                = idx
        self.left_edge          = Grid.edges[idx]
        self.right_edge         = Grid.edges[idx + 1]
        self.center             = Grid.centers[idx]
        self.current_time_step  = Grid.current_time_step
        self.dt                 = Grid.current_dt
        if hasattr(Grid, 'G'):
            self.neutronics     = 1
            self.G              = Grid.G
            self.v              = Grid.v
        else:
            self.neutronics     = 0
        if hasattr(Grid, 'J'):
            self.delayed        = 1
            self.J              = Grid.J
        else:
            self.delayed        = 0
            
    def queryNeighbors(self):

        self.left_outgoing_term   = todo
        self.right_outgoing_term  = todo
        self.left_incoming_term   = todo
        self.right_incoming_term  = todo
        
    def buildLocalMatrix(self):
        if hasattr(self, 'J'): #if delayed neutronics is enabled
            self.A = df.GetAMultiGroupDelayedINF(self.G,self.J,self.SigmaT,
                                                 self.SigmaF,self.SigmaS,self.v,
                                                 self.nu_prompt,self.nu_delayed,
                                                 self.chi_prompt,self.chi_delayed,
                                                 self.decay,self.beta_frac,self.dt,0)
        if not hasattr(self, 'J'):
            self.A = df.GetAMultiGroupDelayedINF(self.G,0,self.SigmaT,
                                                 self.SigmaF,self.SigmaS,self.v,
                                                 self.nu_prompt,0,
                                                 self.chi_prompt,0,
                                                 0,0,self.dt,0)            
    def neutronicsMaterialProperties(self,aNeuMatLib):
    #initialize neutronics
        self.SigmaT     = aNeuMatLib[1]
        self.SigmaF     = aNeuMatLib[2]
        self.SigmaS     = aNeuMatLib[3]
        self.nu_prompt  = aNeuMatLib[5]
        self.chi_prompt = aNeuMatLib[7]
    
    def delayedMaterialProperties(self,aNeuMatLib):
    #initialize delayed neutronics, requires neutronics to be initialized
        
        self.nu_delayed = aNeuMatLib[6]
        self.beta_frac  = aNeuMatLib[10]
        self.chi_delayed = aNeuMatLib[8]
        self.decay      = aNeuMatLib[9] 
        
    def readStateFromFile(self,file_name_str):
        print('readStateFromFile:todo')
        
        f = open(file_name_str,'r')
        self.soln_vector = f.read()
        pass
    
    def writeStateToFile(self):
        print('writeStateToFile:todo')
        path = './Output'
        pathExists = os.path.exists(path)
        if pathExists:
            pass
        else:
            os.makedirs(path)
        pass
    
        file_name_str = r'out_{}_{}.csv'.format(self.current_time_step,self.idx)    
        with open(file_name_str,'w') as csvfile:
            csvwriter = csv.writer(csvfile)
            
            csvwriter.writerow(self.DataLabels)
            csvwriter.writerows(self.output_data)
       
    def prepareOutputData(self):
        #check what physics is initialized, prepare .csv file labels and data accordingly
        self.DataLabels = []
        num_cols = 0
        if self.neutronics == 1:
            self.DataLabels.append('phi')
            num_cols += 1
            
        if self.delayed == 1:
            self.DataLabels.append('C')
            num_cols += 1
        
        #Likely, one of these will be longer than the other, this needs to be handled
        
        if len(self.phi) > len(self.C):
            larger_vector_len = len(self.phi)
        elif len(self.C) > len(self.phi):
            larger_vector_len = len(self.C)
        else:
            larger_vector_len = len(self.phi)
        
        self.output_data = np.zeros([num_cols,larger_vector_len])
        if num_cols == 1:
            self.output_data[0,:] = self.phi
        elif num_cols == 2:
            try:
                self.output_data[0,:] = self.phi
            except IndexError:
                for g in range(self.G):
                    self.output_data[0,g] = self.phi[g]
            try:
                self.output_data[1,:] = self.C
            except IndexError:
                for j in range(self.J):
                    self.output_data[1,j] = self.C[j]
                    
print("init grid")
test_grid = Grid(10,50,matLib,geo,0,10,20,0)
print("assign eneryg grid properties")
test_grid.neutronicsEnergyGridProperties(matLib)
print("assign delayed neutron properties")
test_grid.delayedEnergyGridProperties(matLib)
print("init a cell")
test_cell = Cell(7,test_grid)
print("assign neutron material properties to the cell")
test_cell.neutronicsMaterialProperties(matLib)
print('assign delayed properties to the cell')
test_cell.delayedMaterialProperties(matLib)
print('build local A matrix for the cell')
test_cell.buildLocalMatrix()