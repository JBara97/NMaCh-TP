import numpy as np
import colorsys
import matplotlib.colors as mc

### Functions for dynamics!
def environmental_update(D,A,pc,pr):
    """ 
    update_equation updates the probability(ies) to be degraded, D, by one time step, 
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

    For numnerical examples see update_equation
    """
    return (1+np.dot(np.maximum(A,0),D))*pc*(1-D)/(1+np.dot(np.maximum(-A,0),D)) + (1-pr)*D

### Utility functions for plotting
def lighten_color(color, amount=0.5):
    """
    Lightens the given color by multiplying (1-luminosity) by the given amount.
    Input can be matplotlib color string, hex string, or RGB tuple. 
    If amount >> 1 then an error might be raised as the corresponding hex string would not be valid

    Parameters
    ---
    color : color string, hex string, or RGB tuple. 
        The color to be modified
    amount : float
        The amount to lighten the colour if amount is in [0,1] otherwise if amount > 1 then the colour is darkened instead

    Returns
    ---
    out : RGB colour
    
    Examples
    ---
    >> lighten_color('g', 0.3)
    >> lighten_color('#F034A3', 0.6)
    >> lighten_color((.3,.55,.1), 0.5)
    """
    try:
        c = mc.cnames[color]
    except:
        c = color
    c = colorsys.rgb_to_hls(*mc.to_rgb(c))
    return colorsys.hls_to_rgb(c[0], 1 - amount * (1 - c[1]), c[2])