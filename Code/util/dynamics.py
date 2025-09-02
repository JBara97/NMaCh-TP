import numpy as np

### Functions for dynamics!
def environmental_update(D,A,pc,pr):
    """ 
    environmental_update updates the probability(ies) to be degraded, D, by one time step, 
    given the interaction matrix A, collapse rate pc and recovery rate pr. 
    Quicker version to go into wrappers and other fns

    Parameters
    ---
    D : numpy.array
        The initial array of degraded probabilities
    A : numpy.array
        The interaction matrix (i.e. adjacency matrix), must be of shape (len(D), len(D)). 
        No self loops, i.e. all diagonal elements should be 0. Must be indexed as A[target,source]
    pc : numpy.array or something broadcastable
        The array/list/collection of collapse rates. Must be in [0,1]. 
    pr : numpy.array or something broadcastable
        The array/list/collection of collapse rates. Must be in [0,1]
        
    Returns
    ---
    D_out : numpy.array

    Examples
    ---
    For an initial vector of zeros, regardless of the interaction matrix, D_out will always be pc
    >>> update_equation(np.zeros(2), np.array( [[0, 4], [-2, 0]] ),np.array([0.2, 0.9]))
    np.array([0.2, 0.9]) 

    For numerical examples see update_equation
    """
    return (1+np.dot(np.maximum(A,0),D))*pc*(1-D)/(1+np.dot(np.maximum(-A,0),D)) + (1-pr)*D

def temp_response(T,T_lim,beta):
    """ 
    Returns the response in collapse rate to a temperature T, 
    for a tipping element with threshold T_lim and temperature sensitivity (i.e. logistic parameter) beta
    
    Parameters
    ---
    T : numpy.array or float
        Temperature(s)
    T_lim : numpy.array or float
        Tipping/critical threshold(s)
    beta : numpy.array or something broadcastable
        Temperature sensitivity(ies) (i.e. logistic parameter)

    Returns
    ---
    out : numpy.array
    """
    return 1/(1+np.exp(beta*(T_lim-T)))

### evolve doesn't currently work... so use something else
def evolve(D_init, A, no_time_steps, 
           order_ST=True, 
           pc = None, tau_c = None, 
           T = None, beta = None, T_lim = None,
           pr = None, tau_r = None,
           ):
    """ 
    NEEDS FIXING, DO NOT USE.
    Wrapper to evolve system up to and including t = no_time_steps, 
    with an initial value of D_init at time t=0. 

    Parameters (required)
    ---
    D_init : numpy.array
        The initial array of degraded probabilities at time step 0
    A : numpy.array
        The interaction matrix (i.e. adjacency matrix), must be of shape (len(D_init), len(D_init)). 
        No self loops, i.e. all diagonal elements should be 0. 
    no_time_steps : int
        The number of time steps to simulate over, 
        i.e. the number of times the update equation is successively applied. 
        Must be >= 1
    
    Parameters (optional):
    order_ST : Bool
        The indexing order of the adjacent matrix. 
        If True, then the adjacency matrix will be read with an index-order that gives A[source, target].
        Otherwise the order will be assumed as A[target, source]
        Default True.
    pc : numpy.array or something broadcastable
        The array/list/collection of collapse rates (without interactions). 
        Will supercede any calls to tau_c or T
        Must be in [0,1]. 
    tau_c :  numpy.array or something broadcastable
        The array/list/collection of collapse time scales. 
    T : array-like or float
        The global temperature(s) at which to evaluate. 
        Must be of length >= no_time_steps. Must have beta and T_lim supplied as well if not None
    beta : array-like or float
        Temperature sensitivities (i.e. logistic parameters) for individual tipping points.
        If array-like must be of length len(D_init). 
    T_lim : array-like or float
        Critical temperature thresholds for individual tipping points.
        If array-like must be of length len(D_init)
    pr : numpy.array or something broadcastable
        The array/list/collection of collapse rates. Will supercede tau_r
    tau_r : numpy.array or something broadcastable
        The array/list/collection of collapse time scales
    
        
    Returns
    ---
    D_out : numpy.array
        The output will be a (len(D_init), no_time_steps+1)-shaped array where the D_out[:,0] = D_init
    Examples
    ---
    """
    ## Initialising D_out
    D_out = np.zeros((len(D_init), no_time_steps))
    D_out[:,0] = D_init
    
    ## Changing pc if needed
    if pc is None:
        if tau_c is None:
            raise ValueError("No collapse time scales given")
            pc = np.ones(len(D_init))/2 ## if nothing is supplied pc will be set to 1/2
        elif (T is None) or (beta is None) or (T_lim is None):
            # if any of T,T_lim,beta are None, then pc = 1/(tau_c+1)
            pc = np.ones((len(D_init), no_time_steps))/(tau_c.reshape(len(D_init),1)+1)
        else:
            pc = np.ones((len(D_init), no_time_steps))
            pc /= (tau_c.reshape(len(D_init),1)+1)
            for t in range(no_time_steps):
                pc[:,t] *= temp_response(T[t], T_lim, beta)
    if np.any(pc>1) or np.any(pc<0):
        raise ValueError("Collapse probability pc not in [0,1]")
    ## Changing pr if needed
    if pr is None:
        if tau_r is None:  
            raise ValueError("No recovery time scales given")
        else:
            pr = 1/(tau_r+1)
    if np.any(pr>1) or np.any(pr<0):
        raise ValueError("Recovery probability pr not in [0,1]")
    if order_ST:
        A = A.transpose()

    for t in range(1,no_time_steps):
        D_out[:,t] = environmental_update(D_out[:,t-1],A,pc[:,t-1],pr)

    return D_out