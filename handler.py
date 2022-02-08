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
        self.dr, self.edges,self.centers = df.create_grid(X,NumX)
        self.current_time_step = 0
        self.NumX = NumX
        #build surface area and volumes for the grid
        assert self.geo == 0 or self.geo == 1 or self.geo == 2
        #0: slab, 1: cylinder, 2: sphere
        self.S = np.zeros(NumX+1) #surface area on the edges
        self.V = np.zeros(NumX)   #volume of the cells
        if self.geo == 0:
            self.S[:] = 1
            self.V = np.ones(NumX)*self.dr
        elif self.geo == 1:
            self.S = 4*np.pi*self.edges**2
            self.V = 4/3*np.pi*(self.edges[1:]**3 - self.edges[:-1]**3)
        elif self.geo == 2:
            self.S = 2*np.pi*self.edges
            self.V = np.pi*(self.edges[1:]**2 - self.edges[:-1]**2)
        
        #Define surface area div by Volume for each cell
        self.S_over_V = np.zeros(NumX + 1,dtype = object)
        
        for i in range(NumX + 1):
            if i == 0: #On the left boundary
                if geo == 1 or geo == 2:
                    self.S_over_V[i] = np.zeros(self.G)
                else:
                    self.S_over_V[i] = 2*(self.S[i])/self.V[i]
        
            elif i == self.NumX: #On the right boundary
                self.S_over_V[i] = 2*(self.S[i])/self.V[i - 1]
        
            else: #In between
                if i == 1 and (geo == 1 or geo == 2):
                    self.S_over_V[i] = 2/(self.S[i])/self.V[i]
                else:
                    self.S_over_V[i] = 2/((self.V[i - 1]/(self.S[i-1])) + self.V[i]/(self.S[i]))
        
        #Time Grid
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
        self.decay  = aNeuMatLib[9] #property of the delayed neutron grid, discrete decay constant vector
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
        #initialize different physics
        if hasattr(Grid, 'G'):
            self.neutronics     = 1
            self.G              = Grid.G
            self.v              = Grid.v
            self.phi            = np.zeros(self.G)
        else:
            self.neutronics     = 0
        if hasattr(Grid, 'J'):
            self.delayed        = 1
            self.J              = Grid.J
            self.C              = np.zeros(self.J)
        else:
            self.delayed        = 0
        #determine if this cell is a boundary cell
        if idx == 0:
            self.is_left_boundary   = 1
            self.is_right_boundary  = 0
        elif idx == Grid.NumX:
            self.is_right_boundary  = 1
            self.is_left_boundary   = 0
        else:
            self.is_left_boundary   = 0
            self.is_right_boundary  = 0
            
    def queryNeighborsFlux(self):
        
        if self.is_left_boundary == 1:
            self.right_soln_vector  = self.readStateFromFile(self, \
                 r'out_{}_{}.csv'.format(self.current_time_step,self.idx - 1))
            pass
        
        elif self.is_right_boundary == 1:
            self.left_soln_vector   = self.readFluxStateFromFile(self, \
                 r'out_{}_{}.csv'.format(self.current_time_step,self.idx - 1))
            pass
        else:            
            self.left_soln_vector   = self.readFluxStateFromFile(self, \
                 r'out_{}_{}.csv'.format(self.current_time_step,self.idx - 1))
            
            self.right_soln_vector  = self.readStateFromFile(self, \
                 r'out_{}_{}.csv'.format(self.current_time_step,self.idx - 1))
        
    def buildLocalInfMedMatrix(self):
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
      
    def buildLocalFullMatrixBC(self,neighborA,neighborVec):
        if self.is_left_boundary or self.is_right_boundary == 1:
            
            self.bigA = np.zeros([2*self.G,2*self.G])
            if self.is_left_boundary == 1:
                self.bigA[:self.G,:self.G] = self.A
                self.bigA[self.G:,self.G:] = neighborA # for right neighbor
                self.bigA[self.G:,:self.G] = np.diag(neighborVec)
                pass
            elif self.is_right_boundary == 1:
                pass
            else:
                print('Error in Cell.buildLocalMatrix(), cell is both boundary nodes!')
                assert 0
        else:
            self.bigA = np.zeros([3*self.G,3*self.G])
            
    def buildLocalFullMatrix(self,LeftNeighborA,LeftNeighborvec,RightNeighborA,RightNeighborVec):
        self.bigA = np.zeros([3*self.G,3*self.G])
        
        pass
    
    def neutronicsMaterialProperties(self,aNeuMatLib):
    #initialize neutronics
        self.SigmaT     = aNeuMatLib[1]
        self.SigmaF     = aNeuMatLib[2]
        self.SigmaS     = aNeuMatLib[3]
        self.nu_prompt  = aNeuMatLib[5]
        self.chi_prompt = aNeuMatLib[7]
        self.D          = 1/(3*self.SigmaT)
    
    def delayedMaterialProperties(self,aNeuMatLib):
    #initialize delayed neutronics, requires neutronics to be initialized
        
        self.nu_delayed = aNeuMatLib[6]
        self.beta_frac  = aNeuMatLib[10]
        self.chi_delayed = aNeuMatLib[8]
        self.decay      = aNeuMatLib[9] 
        
    def readFluxStateFromFile(self,flux_file_name_str):
        print('readFluxStateFromFile:todo')
        
        f = open(flux_file_name_str,'r')
        self.phi = f.read()
        pass
    
    def writeFluxStateToFile(self):
        print('writeFluxStateToFile:todo')
        path = './Output/'
        pathExists = os.path.exists(path)
        if pathExists:
            pass
        else:
            os.makedirs(path)
            
        self.flux_file_name_str = r'Flux_out_{}_{}.dat'.format(self.current_time_step,self.idx)
        self.flux_file_name_str = path + self.flux_file_name_str
        with open(self.flux_file_name_str,'w') as outfile:
            for line in self.phi:
                outfile.write("{}\n".format(line))
                
    def readPrecursorStateFromFile(self,precursor_file_name_str):
        print('readPrecursorStateFromFile: todo')
        f = open(precursor_file_name_str,'r')
        self.C = f.read()
        pass
    
    def writePrecursorStateToFile(self):
        print('writePrecursorStateToFile: Todo')
        path = './Output/'
        pathExists = os.path.exists(path)
        if pathExists:
            pass
        else:
            os.makedirs(path)
            
        self.precursor_file_name_str = r'Pre_out_{}_{}.dat'.format(self.current_time_step,self.idx)
        self.precursor_file_name_str = path + self.precursor_file_name_str
        with open(self.precursor_file_name_str,'w') as outfile:
            for line in self.C:
                outfile.write("{}\n".format(line))
maxX = X
NumX = Nx
aNeuMatLib = matLib
geo = geo
startTime = 0
finalTime = T
numTimePoints = 20
method = 0

print("init grid")
test_grid = Grid(maxX,NumX,matLib,geo,startTime,finalTime,numTimePoints,method)
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
test_cell.buildLocalInfMedMatrix()