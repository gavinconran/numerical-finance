from scipy.linalg import solve

def generateMatrix(N, sigma):
    """ Computes the matrix for the diffusion equation with backward Euler
        Dirichlet condition at i=0, Neumann at i=-1
    
    Parameters:
    ----------
    T: array of float
        Temperature at current time step
    sigma: float 
        alpha*dt/dx^2
    
    Returns:
    -------
    A: 2D numpy array of float
        Matrix for diffusion equation
    """
   
    # Setup the diagonal
    d = np.diag(np.ones(N-2)*(2+1./sigma))
    
    # Consider Neumann BC
    #d[-1,-1] = 1+1./sigma
    
    # Setup upper diagonal
    ud = np.diag(np.ones(N-3)*-1, 1)
    
    # Setup lower diagonal
    ld = np.diag(np.ones(N-3)*-1, -1)
    
    A = d + ud + ld
    
    return A

def generateRHS(T, sigma, qdx):
    """ Computes right-hand side of linear system for diffusion equation
        with Forward Euler
    
    Parameters:
    ----------
    T: array of float
        Temperature at current time step
    sigma: float
        alpha*dt/dx^2
    qdx: float
        flux at right boundary * dx
    
    Returns:
    -------
    b: array of float
        Right-hand side of diffusion equation with backward Euler
    """
    
    b = T[1:-1]*1./sigma
    # Consider Dirichlet BC
    b[0] += T[0]
    # Consider Neumann BC
    #b[-1] += qdx
    
    return b

def implicit_btcs(T, A, nt, sigma, qdx):
    """ Advances diffusion equation in time with implicit central scheme
   
    Parameters:
    ----------
    T: array of float
        initial temperature profile
    A: 2D array of float
        Matrix with discretized diffusion equation
    nt: int
        number of time steps
    sigma: float
        alpha*td/dx^2
        
    qdx: float
        flux at right boundary * dx
    Returns:
    -------
    T: array of floats
        temperature profile after nt time steps
    """
    
    for t in range(nt):
        Tn = T.copy()
        b = generateRHS(Tn, sigma, qdx)
        # Use numpy.linalg.solve
        T_interior = solve(A,b)
        T[1:-1] = T_interior 
        # Enforce Neumann BC (Dirichlet is enforced automatically)
        #T[-1] = T[-2] + qdx

    return T

