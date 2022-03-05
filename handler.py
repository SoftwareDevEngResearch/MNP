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
    
    T0  = float(f.readline())   #initial time
    T   = float(f.readline())   #final time
    NT  = int(f.readline())     #Number of time steps
    X   = float(f.readline())   #radius of object
    Nx  = int(f.readline())     #number of space grid points
    BCL = float(f.readline())   #left reflectance boundary condition
    BCR = float(f.readline())   #right reflectance boundary condition
    geo = int(f.readline())     #geometry 1-slab, 2-cyl, 3-sphere
    alpha = [BCL,BCR]
    assert geo == 0 or geo == 1 or geo == 2
    
    pass

################# placeholders #########################

todo = 'boggle'

matLib = np.array((E,SigmaT,SigmaF,SigmaS,v,nu_prompt,nu_delayed,chi_prompt,chi_delayed,decay,beta_frac),dtype = object)
#                 [0,1     ,2     ,3     ,4,5,       ,6,        ,7         ,8          ,9    ,10]
########################################################
class Grid:
    #creates the computational grid, 1d for now.
    #   initializes a time to be zero and calls a function to define points on the time grid
    #   
    def __init__(self,maxX,NumX,aNeuMatLib,geo,startTime,finalTime,numTimePoints,method,alpha): 
        #create the grid
        
        #These are defined for the computataional spatial grid
        self.geo = geo
        self.dr, self.centers,self.edges = df.create_grid(X,NumX)
        self.current_time_step = 0
        self.numTimePoints = numTimePoints
        self.NumX = NumX
        #build surface area and volumes for the grid
        assert self.geo == 0 or self.geo == 1 or self.geo == 2
        #0: slab, 1: cylinder, 2: sphere
        self.S = np.zeros(NumX+1) #surface area on the edges
        self.V = np.zeros(NumX)   #volume of the cells
        self.cells_built = 0
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
        self.alpha = alpha
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
        
    def buildCells(self):
        self.cellsList = []
        for i in range(self.NumX):
            self.cellsList.append(Cell(i,self))
        self.cells_built = 1
        pass

    def buildA_Total(self):
        if self.cells_built != 1:
            self.buildCells()
        assert self.cells_built == 1
        self.A_total = np.zeros([self.NumX*self.G,self.NumX*self.G])
        for i in range(self.NumX):
            if i == 0:
                self.cellsList[i].buildRow()
                self.A_total[:self.G,:2*self.G] = self.cellsList[i].A_contribution.copy()
            elif i == self.NumX - 1:
                self.cellsList[i].buildRow()
                self.A_total[(self.NumX*self.G) - self.G:,(self.NumX*self.G) - 2*self.G:] = \
                    self.cellsList[i].A_contribution.copy()
            else:
                self.cellsList[i].buildRow()
                self.A_total[i*self.G:(i + 1)*self.G, (i-1)*self.G:(i + 2)*self.G] = self.cellsList[i].A_contribution.copy()
        pass
    
    def buildRHS_Total(self):
        if self.cells_built != 1:
            self.buildCells()
        assert self.cells_built == 1
        
        self.RHS_total =  np.zeros(self.NumX*self.G)
        for i in range(self.NumX):
            self.cellsList[i].buildRHS()
            self.RHS_total[i*self.G:(i+1)*self.G] = self.cellsList[i].RHS.copy()
        pass
    
    def step(self):
        
        assert self.current_time_step < self.numTimePoints - 1
        
        self.buildA_Total()
        self.buildRHS_Total()
        
        newFlux = np.linalg.solve(self.A_total, self.RHS_total)
        for i in range(self.NumX):
            
            self.cellsList[i].initializeFlux(newFlux[i*self.G:(i+1)*self.G])
            self.cellsList[i].C = df.MultiGroupCStepper(\
                self.J, self.cellsList[i].nu_delayed, self.cellsList[i].SigmaF,\
                    newFlux[i*self.G:(i+1)*self.G], self.current_dt, self.decay\
                    , self.cellsList[i].C, self.cellsList[i].beta_frac)
            
        
        self.current_time_step += 1
        self.global_time += self.current_dt

        
        if self.current_time_step >= self.numTimePoints - 1:
            print('reached end of time steps, t = {}'.format(self.global_time))    
        else:
            self.current_dt = self.dts[self.current_time_step]
            for i in range(self.NumX):
                self.cellsList[i].advance(self.current_dt)
                self.cellsList[i].buildLocalInfMedMatrix()
        
        
    def buildPHI_Total(self):
        self.PHI = np.zeros(self.NumX*self.G)
        for i in range(self.NumX):
            self.PHI[i*self.G:(i+1)*self.G] = self.cellsList[i].phi.copy()
    
    def buildC_Total(self):
        self.C = np.zeros(self.NumX*self.J)
        for i in range(self.NumX):
            self.C[i*self.J:(i+1)*self.J] = self.cellsList[i].C.copy()
        
    def buildN_Total(self):
        self.N_Total = np.zeros(self.NumX)
        for i in range(self.NumX):
            self.N_Total[i] = sum(np.divide(self.cellsList[i].phi,self.v))
    
class Cell:
    #creates an individual cell of the grid, needs to have a parent grid.
        
    def __init__(self,idx,Grid):
        
        self.idx                = idx
        self.left_edge          = Grid.edges[idx]
        self.right_edge         = Grid.edges[idx + 1]
        self.center             = Grid.centers[idx]
        self.current_time_step  = Grid.current_time_step
        self.dt                 = Grid.current_dt
        self.dr                 = Grid.dr
        self.S_over_V_L         = Grid.S_over_V[idx]
        self.S_over_V_R         = Grid.S_over_V[idx + 1]
        #initialize different physics
        if hasattr(Grid, 'G'):
            self.neutronics     = 1
            self.G              = Grid.G
            self.v              = Grid.v
            self.phi            = np.zeros(self.G)
            self.RHS            = np.zeros(self.G)
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
            self.alpha              = Grid.alpha
        elif idx == Grid.NumX - 1:
            self.is_right_boundary  = 1
            self.is_left_boundary   = 0
            self.alpha              = Grid.alpha
        else:
            self.is_left_boundary   = 0
            self.is_right_boundary  = 0
        

        
    def buildLocalInfMedMatrix(self):
        if not hasattr(self, 'SigmaT'):
            print('***************')
            print('Initialize physics with neutronicsMaterialProperties() and if necessary delayedMaterialProperties()')
            print('***************')
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
      
    def buildRow(self):
        #builds the local row, to send to the global matrix
        if not hasattr(self, 'A'):
            self.buildLocalInfMedMatrix()
        else:
            pass
        if self.is_left_boundary == 0 and self.is_right_boundary == 0:
            diag_idx_left         = np.array(np.diag_indices(self.G)[0])
            diag_idx_middle       = np.array(np.diag_indices(self.G)[0]) + self.G
            diag_idx_right        = np.array(np.diag_indices(self.G)[0]) + 2*self.G
            self.A_contribution = np.zeros([self.G,3*self.G])
            self.A_contribution[:,self.G:2*self.G]                    = self.A
            self.A_contribution[diag_idx_left,diag_idx_middle]        += self.SigmaJ_L + self.SigmaJ_R
            self.A_contribution[diag_idx_left,diag_idx_left]          = -1*self.SigmaJ_L
            self.A_contribution[diag_idx_left,diag_idx_right]         = -1*self.SigmaJ_R
        
        if self.is_left_boundary == 1 and self.is_right_boundary == 0:
            diag_idx_left   = np.array(np.diag_indices(self.G)[0])
            diag_idx_middle = np.array(np.diag_indices(self.G)[0]) + self.G
            self.A_contribution = np.zeros([self.G,2*self.G])
            self.A_contribution[:,:self.G] = self.A
            self.A_contribution[diag_idx_left,diag_idx_middle] = -1*self.SigmaJ_R
            self.A_contribution[diag_idx_left,diag_idx_left] += self.SigmaJ_R + self.SigmaJ_L
            self.A_contribution[diag_idx_left,diag_idx_left] -= np.multiply(self.SigmaJ_L,self.B_L)
        
        if self.is_right_boundary == 1 and self.is_left_boundary == 0:
            diag_idx_left = np.array(np.diag_indices(self.G)[0])
            diag_idx_right = np.array(np.diag_indices(self.G)[0]) + self.G
            self.A_contribution = np.zeros([self.G,2*self.G])
            self.A_contribution[diag_idx_left,diag_idx_left]  = -1*self.SigmaJ_L
            self.A_contribution[:,self.G:] = self.A
            self.A_contribution[diag_idx_left,diag_idx_right] += self.SigmaJ_R + self.SigmaJ_L
            self.A_contribution[diag_idx_left,diag_idx_right] -= np.multiply(self.SigmaJ_R,self.B_R)
            
    def buildRHS(self):
        dtInv = 1/self.dt
        for g in range(self.G):
            Delayed_contrib = np.zeros(self.G)
            for j in range(self.J):
                Delayed_contrib[g] += self.chi_delayed[g,j] * self.decay[j] * \
                    (dtInv + self.decay[j])**-1 * dtInv * self.C[j]
            self.RHS[g] = ( dtInv/self.v[g] * self.phi[g] + Delayed_contrib[g] )
        pass
            
    def initializeFlux(self,vector = ''):
        if vector == '':
            vector = np.zeros(self.G)
        self.phi = vector
        
    def initializeC(self,vector = ''):
        if vector == '':
            vector = np.zeros(self.J)
        self.C = vector        
    
    def advance(self, dt):
        self.current_time_step += 1
        self.dt = dt
    
    def neutronicsMaterialProperties(self,aNeuMatLib):
    #initialize neutronics
        self.SigmaT         = aNeuMatLib[1]
        self.SigmaF         = aNeuMatLib[2]
        self.SigmaS         = aNeuMatLib[3]
        self.nu_prompt      = aNeuMatLib[5]
        self.chi_prompt     = aNeuMatLib[7]
        self.D              = 1/(3*self.SigmaT)
        self.DS_over_V_L    = self.S_over_V_L * self.D
        self.DS_over_V_R    = self.S_over_V_R * self.D
        self.SigmaJ_L       = self.DS_over_V_L / self.dr
        self.SigmaJ_R       = self.DS_over_V_R / self.dr 
        if self.is_left_boundary == 1 and self.is_right_boundary == 0:
            self.B_L = (1-self.alpha[0])/(1+self.alpha[0]) * self.dr/self.D/4 + 1
            self.B_L = self.B_L**-1
        if self.is_right_boundary == 1 and self.is_left_boundary == 0:
            self.B_R = (1-self.alpha[1])/(1+self.alpha[1]) * self.dr/self.D/4 + 1
            self.B_R = self.B_R**-1
            
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
if __name__ == "__main__":
    
    maxX = X
    NumX = Nx
    aNeuMatLib = matLib
    geo = geo
    startTime = T0
    finalTime = T
    numTimePoints = NT
    method = 1
    
    print("init grid")
    test_grid = Grid(maxX,NumX,matLib,geo,startTime,finalTime,numTimePoints,method,alpha)
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
    print('build row of global A matrix')
    test_cell.buildRow()
    print('build cells with method of Grid')
    test_grid.buildCells()
    print('initalize the cells with material properties')
    for i in range(NumX):
        test_grid.cellsList[i].neutronicsMaterialProperties(aNeuMatLib)
        test_grid.cellsList[i].delayedMaterialProperties(aNeuMatLib)
    print('build A_total with method')
    test_grid.buildA_Total()
    t = test_grid.A_total
    print('set value of RHS[0] with a method')
    test_grid.cellsList[0].initializeFlux(v)
    print('build RHS with method')
    test_grid.buildRHS_Total()
    print('step forward in time twice with a method')
    
    print(test_grid.current_dt,test_grid.cellsList[0].dt)
    
    phiHist = np.empty([test_grid.NumX*test_grid.G,test_grid.numTimePoints])
    qHist   = np.zeros([test_grid.NumX*test_grid.G,test_grid.numTimePoints])
    matHist = np.zeros(3,dtype = object)
    CHist   = np.zeros([test_grid.NumX*test_grid.J,test_grid.numTimePoints])
    
    test_grid.buildPHI_Total()
    test_grid.buildC_Total()
    phiHist[:,0] = test_grid.PHI
    matHist[0] = test_grid.A_total
    qHist[:,0] = test_grid.RHS_total
    CHist[:,0] = test_grid.C
    test_grid.step()
    print(test_grid.current_dt, test_grid.cellsList[0].dt)
    
    test_grid.buildPHI_Total()
    test_grid.buildC_Total()
    phiHist[:,1] = test_grid.PHI
    matHist[1] = test_grid.A_total
    qHist[:,1] = test_grid.RHS_total
    CHist[:,1] = test_grid.C
    test_grid.step()
    print(test_grid.current_dt,test_grid.cellsList[0].dt)
    
    test_grid.buildPHI_Total()
    test_grid.buildC_Total()
    qHist[:,2] = test_grid.RHS_total
    phiHist[:,2] = test_grid.PHI
    matHist[2] = test_grid.A_total
    CHist[:,2] = test_grid.C
