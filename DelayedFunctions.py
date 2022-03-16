from numpy import zeros,pi,dot,copy,shape, arange,logspace,nonzero, \
    size, append, matrix, diag, flip, array, diff, cumsum
from numpy import sum as sum2
from numpy.linalg import solve
from scipy.linalg import svd, eig
from numba import jit

@jit(nopython = True)
def GetAMultiGroupDelayedINF(G,J,SigmaT,SigmaF,SigmaS,v,nu_prompt,nu_delayed,
                             chi_prompt,chi_delayed,decay,beta_frac,dt,Radius): 
    #Multi Group with Delayed Neutrons for an infinite medium problem with buckling
    #set radius = 0 for no buckling
    A           = zeros((G,G))
    R           = zeros((G,G))
    S           = zeros((G,G))
    F           = zeros((G,G))
    SigmaL      = zeros(G)
    sumFactor3  = zeros(J)
    dtInv = dt**-1 #dt is a float, not an array
    
    if Radius != 0:
        for g in range(G):
            SigmaL[g] = 1/3/SigmaT[g] * (pi/Radius)**2
    else:
        pass
    
    for g in range(G):
        R[g,g] = (v[g]*dt)**-1 + (SigmaT[g] + SigmaL[g])
        for gprime in range(G):
            S[gprime,g] = SigmaS[gprime,g]
            sumFactor3 = 0
            for j in range(J):
                sumFactor3 += chi_delayed[g,j]*decay[j]*(dtInv + decay[j])**-1*beta_frac[j]*nu_delayed[gprime]
            F[g,gprime] = ( chi_prompt[g,gprime] * nu_prompt[gprime] + sumFactor3 ) * SigmaF[gprime]
        
    A = R - S - F
    
    return A #gets called to build local A matrix.

@jit(nopython = True)
def MultiGroupCStepper(J,nu_delayed,SigmaF,phiNew,dt,decay,cOld,beta_frac):
    
    #Generates Precursors after the flux has been updated
    #To be run after all the other physics is done
    #e.g. Timestep -> Diffusion -> Precursor Generation
    
    total = (nu_delayed*SigmaF).dot(phiNew)
    cNew = zeros(J)
    for j in range(J):
        cNew[j] = (1 + dt*decay[j])**-1 * cOld[j] + dt/(1 + decay[j]*dt)*beta_frac[j] * total
            
    return cNew

@jit(nopython = True)
def create_grid(R,I):
    """Create the cell edges and centers for a 
    domain of size R and I cells
    Args:
        R: size of domain
        I: number of cells
        
    Returns:
        Delta_r: the width of each cell
        edges: the cell edges of the grid
        centers: the cell centers of the grid
    """
    Delta_r = float(R)/I
    centers = arange(I)*Delta_r + 0.5*Delta_r
    edges   = arange(I+1)*Delta_r
    
    return Delta_r, centers, edges

