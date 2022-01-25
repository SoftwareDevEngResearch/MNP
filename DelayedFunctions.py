from numpy import zeros,pi,dot,copy,shape, arange,logspace,nonzero, \
    size, append, matrix, diag, flip, array, diff, cumsum
from numpy import sum as sum2
from numpy.linalg import solve
from scipy.linalg import svd, eig

def GetAMultiGroupDelayedINF(G,J,SigmaT,SigmaF,SigmaS,v,nu_prompt,nu_delayed,
                             chi_prompt,chi_delayed,decay,beta_frac,dt,Radius): 
    #Multi Group with Delayed Neutrons for an infinite medium problem with buckling
    #set radius = 0 for no buckling
    A           = zeros([G,G])
    R           = zeros([G,G])
    S           = zeros([G,G])
    F           = zeros([G,G])
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

def MultiGroupDelayedStepperINF(G,J,SigmaT,SigmaF,SigmaS,v,nu_prompt,nu_delayed,
                               chi_prompt,chi_delayed,decay,beta_frac,dt,phi,C,
                               Radius,qFunction):
    
    #inputs:
        # G - Number of Energy Groups, int
        # J - Number of Precursor Groups, int
        # SigmaT - Array of Total Cross Sections, array of double
        # SigmaS - Matrix of Scattering Cross Section, Matrix of double
        # v - Array of Neutron Speed in the Energy Groups
        # nu_prompt - Array of probability of collision creating a prompt neutron, float
        # nu_delayed - Array of probability of collision creating a delayed neutron, float
        # chi_prompt - Array of 
        # chi_delayed - 
        # decay - Array of decay time constants for precursor groups, float
        # beta_frac - Array of fractions of collisions that create a delayed neutron
        # dt - Timestep to use for the backwards euler timestep
        # phi - array of float time solutions of the neutron flux
        # C - array of float, time solutions of the precursor concentration
        # t - array of float, time vector
        # Radius - Radius to use for buckling approximation, float, cm
        # qFunction - Function that returns some prescribed flux source in energy group g
        
    #outpts:
        # phiNew    - New neutron flux after one time step
        # cNew      - New Precursor solution after one time step
    
    phiOld  = phi.copy()
    phiNew  = zeros(len(phiOld))
    cOld    = C.copy()
    qStar   = zeros(G)
    
    dtInv = dt **-1
    A = GetAMultiGroupDelayedINF(G,J,SigmaT,SigmaF,SigmaS,v,nu_prompt,nu_delayed,
                         chi_prompt,chi_delayed,decay,beta_frac,dt,phi,C,Radius)
    
    for g in range(G):
        sumFactor2 = zeros(G)
        for j in range(J):
               sumFactor2[g] += chi_delayed[g,j] * decay[j] * (dtInv + decay[j])**-1 * dtInv * cOld[j]
        qStar[g] =  qFunction(g) + ( dtInv/v[g] * phiOld[g] + sumFactor2[g] )
        
    
    phiNew = solve(A,qStar)
        
    return phiNew

def MultiGroupCStepper(J,nu_delayed,SigmaF,phiNew,dt,decay,cOld,beta_frac):
    
    #Generates Precursors after the flux has been updated
    #To be run after all the other physics is done
    #e.g. Timestep -> Diffusion -> Precursor Generation
    
    total = (nu_delayed*SigmaF).dot(phiNew)
    cNew = zeros(J)
    for j in range(J):
        cNew[j] = (1 + dt*decay[j])**-1 * cOld[j] + dt/(1 + decay[j]*dt)*beta_frac[j] * total
            
    return cNew

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

def StackA(ASoFar,AToBeAdded,Sigma_J_Plus,Sigma_J_Minus,G,space_step):
    #Given a new A matrix, incorpirate it into the existing A matrix along the diagonal.
    #needs the Sigma_J plus for the last energy group and minus for the first energy group
    
    offset = G*space_step
    
    if space_step == 0:
        NewA = AToBeAdded
        NewA[G-1,2*G - 1] = -1*Sigma_J_Plus
        NewA[2*G - 1,G-1] = -1*Sigma_J_Minus
    else:
        NewA = zeros([(space_step + 2)* G,(space_step + 2)* G ])
        
        NZIdx = nonzero(ASoFar)
        NewA[NZIdx] = ASoFar[NZIdx]
        
        for i in range(G):
            for j in range(G):
                NewA[i + space_step*G,j + space_step*G] = AToBeAdded[offset + i, offset + j] # Add values in the GxG removal/scattering/fission matrix
        for i in range(G):
            NewA[offset + i, offset + G + i] = AToBeAdded[offset + i,offset + G + i] #Add values in the diagonal Sigma_J matrixes.
            NewA[offset + G + i, offset + i] = AToBeAdded[offset + G + i, offset + i] #Sigma_J_Minus
                
        NewA[G + offset - 1,G + offset + G - 1] = -1*Sigma_J_Plus
        NewA[G + offset + G - 1,G + offset - 1] = -1*Sigma_J_Minus 
    
    return NewA

def BuildPhiMatrix(ListOfPhi):
    
    #Takes a list of Phi vectors and creates a matrix
    #e.g.
    #                      [  -  -  - ]
    #  [[y0],[y1],[y2]] -> [ y0 y1 y2 ]
    #                      [  -  -  - ] 
    # where y0 = y[y00,y01,y02, ... ,y0G], 0 is the cell number, G is the number of energy bins
    #
    numCells = len(ListOfPhi)
    numGroups = len(ListOfPhi[0])
    PhiMatrix = zeros([numCells,numGroups])
    for Cell in range(numCells):
        PhiMatrix[Cell,:] = ListOfPhi[Cell].copy()
        
    return PhiMatrix

def Generic_1D_Diffusion(phiMatrix,D,geometry,R,I,):
    #
    #Diffuses the flux solution
    #
    #inputs:
        # PhiMatrix - I x G matrix of the flux solution in each cell
        # D         - 1 x G Diffusion Constant for the flux energy groups
        # geometry  - Integer flag for problem geometry
            # 0 = rectangular
            # 1 = cylindrical
            # 2 = spherical
        # R         - Characteristic Length of the problem
        # I         - Number of Cells in the Problem
        #
    # The left boundary condition is always reflecting, it's a line of symmetry
    # The right boundary condition depends on the values in the vector BC
        # BC = [a,b,c]
        # For Reflecting:
            #   
        # For Albedo, with reflectance a:
            #
        # For vacuum, albedo with reflectance 0:
            #
    
    Delta_r, centers, edges = create_grid(R,I)
    
    #define surface areas and volumes
    
    assert( (geometry==0) or (geometry == 1) or (geometry == 2))
    if (geometry == 0):
        #in slab it's 1 everywhere except at the left edge
        S = 0.0*edges+1
        S[0] = 0.0 #to enforce Refl BC
        #in slab its dr
        V = 0.0*centers + Delta_r
    elif (geometry == 1):
        #in cylinder it is 2 pi r
        S = 2.0*pi*edges
        #in cylinder its pi (r^2 - r^2)
        V = pi*( edges[1:(I+1)]**2 
                   - edges[0:I]**2 )
    elif (geometry == 2):
        #in sphere it is 4 pi r^2
        S = 4.0*pi*edges**2
        #in sphere its 4/3 pi (r^3 - r^3)
        V = 4.0/3.0*pi*( edges[1:(I+1)]**3
                   - edges[0:I]**3 )
    
    phiMatrixNew = zeros(shape(phiMatrix))
    
    
    
    return phiMatrixNew

def MultiGroupPhi(PhiArray):
    #Like MultiGroupPhi, but doesn't have zeroes between the groups
    G = len(PhiArray[0])
    length = size(PhiArray)
    # Phi = zeros(G*length)
    # for i in range(length):
    #     for g in range(G):
    #         Phi[i+length*g] = PhiArray[i][g]
    # for some reason np.append is extremely slow, cast these vectors to
    # regular lists and do regular python append
    Phi = []
    for i in range(length):
        # Phi = append(Phi,PhiArray[i])
        Phi.extend(PhiArray[i].tolist())
        # print(Phi,PhiArray[i])  
        
    #cast back to a numpy array and return
    return array(Phi)

def MultiGroupPhiInv(Phi,G):
    #Like MultiGroupPhiInv but for the corresponding MultiGroupPhi2
    length = int((len(Phi))/G)
    # PhiArray = zeros([length,G])
    # for g in range(G):
    #     for i in range(length):
    #         PhiArray[i,g] = Phi[i + length*g]
    
    PhiArray = Phi.flatten(order = 'F').reshape((length,G),order='C')
    
    return PhiArray


def DMD(phiHist, cHist, t):
    #stack the flux and precursor concentration data into a (G+J) x N matrix

    phi = append(phiHist,cHist,0)
    # print(phi.shape)
    #compute dt, we assume constant for this algorithm
    dt = t[1] - t[0]
    #USV decomposition
    u,s,v = svd(flip(phi[:,0:-1]),full_matrices = False) # y_
    u = matrix(u)  
    #check what values in S are small enough to be neglected
    sPos = s[s/cumsum(s)>1e-14]#s[(1-np.cumsum(s)/np.sum(s)) > 1e-15]
    # print(sPos,s,np.cumsum(s)/np.sum(s))
    
    u = u[:,0:len(sPos)]
    v = v[0:len(sPos),:]
    
    #cast v to matrix to use getH
    v = matrix(v)
    #build S matrix by inverting values in sPos and making diagonal from vector
    sInv = sPos**-1
    S = diag(sInv)
    
    Stilde0 = dot(u.getH(),flip(phi[:,1:])) #y+
    Stilde1 = dot(Stilde0,v.getH())
    Stilde  = dot(Stilde1,S)
    
    eigs,vecs = eig(Stilde)
    eigs = (1.0 - 1.0/eigs)/dt
    
    return eigs, vecs

def DMD_Prompt(phiHist, t):

    phi = phiHist.copy()
    # print(phi.shape)
    #compute dt, we assume constant for this algorithm
    dt = t[2] - t[1]
    #USV decomposition
    u,s,v = svd(flip(phi[:,0:-1]),full_matrices = False)
    u = matrix(u)  
    #check what values in S are small enough to be neglected
    sPos = s[s/cumsum(s) > 1e-14]#s[(1-np.cumsum(s)/np.sum(s)) > 1e-15]
    # print(sPos,s,np.cumsum(s)/np.sum(s))
    
    u = u[:,0:len(sPos)]
    v = v[0:len(sPos),:]
    
    #cast v to matrix to use getH
    v = matrix(v)
    #build S matrix by inverting values in sPos and making diagonal from vector
    sInv = sPos**-1
    S = diag(sInv)
    
    Stilde0 = dot(u.getH(),flip(phi[:,1:]))
    Stilde1 = dot(Stilde0,v.getH())
    Stilde  = dot(Stilde1,S)
    
    eigs,vecs = eig(Stilde)
    eigs = (1.0 - 1.0/eigs)/dt
    
    return eigs, vecs

def DMD_VTS(phiHist, cHist, t):
    #Use dynamic mode decomposition to determine the eigenvalues of a system
    #using a data matrix.
    # This version permits a Variable Time Step.
    
    #stack the flux and precursor concentration data into a (G+J)*I x N matrix
    dt = diff(t)
    phi = append(phiHist,cHist,0)
    
    dPhi = diff(phi,axis=1)/dt
    # print(phi.shape)
    #USV decomposition
    u,s,v = svd(phi[:,1:],full_matrices = False)
    u = matrix(u)  
    #check what values in S are small enough to be neglected
    sPos = s[s/cumsum(s)>1e-14]
    # print(sPos,s,np.cumsum(s)/np.sum(s))
    
    u = u[:,0:len(sPos)]
    v = v[0:len(sPos),:]
    
    #cast v to matrix to use getH
    v = matrix(v)
    #build S matrix by inverting values in sPos and making diagonal from vector
    sInv = sPos**-1
    S_function = diag(sInv)
    
    Stilde0 = dot(u.getH(),dPhi) #[:,1:]
    Stilde1 = dot(Stilde0,v.getH())
    Stilde  = dot(Stilde1,S_function)
    
    eigs,vecs = eig(Stilde)
    # eigs = (1.0 - 1.0/eigs)/dt
    
    return eigs, vecs