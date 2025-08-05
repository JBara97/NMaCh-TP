import numpy as np
import networkx as nx
import colorsys
import matplotlib.colors as mc
### Before ../12/24
def environmental_update_oldest(D,A,p_c=0.4,p_r=0.2,q_c=None,q_r=None,alpha=1):
    """ D - array-like of length N
        A - adjacency matrixb NxN; a_ij is i-> j
        p_c,p_r - probability to collapse and recover, resp., [0,1]
        proportional collapse and recovery """
    x = np.dot(np.transpose(A),D)
    k_in = np.sum(A,axis=0)
    K_bar=np.zeros_like(x)
    K_bar[k_in>0] = x[k_in>0]/k_in[k_in>0]
    if q_c != None:
        p_c=q_c*D
    if q_r != None:
        p_r=q_r*(1-D)
    phi = p_c+D*(1-p_c-p_r)#q_c*D + (1-q_c*D-q_r*(1-D))*D#(1+q_c-q_r)*D + (q_r-q_c)*D**2
    # have only included the provision when alpha is independent of the upstream neighbour (i.e. alpha_{ji} = alpha_i)
    D_updated = phi + K_bar*(alpha-phi)
    return D_updated

### Version 17/12/24
def environmental_update_old(D,A,p_c=0.4,p_r=0.2,q_c=None,q_r=None,alpha=1):
    """ D - array-like of length N
        A - adjacency matrix NxN; a_ij is i-> j. OR a tuple of (A^+,A^-)
        p_c,p_r - probability to collapse and recover, resp., [0,1]
        alpha - number or array-like of length N OR a tuple of (alpha_pos,alpha_neg) """
    if type(A) != tuple:
        A_pos = A
        A_neg = A
    else:
        A_pos,A_neg = A
    if type(alpha) != tuple:
        alpha_pos = alpha
        alpha_neg = alpha
    else:
        alpha_pos,alpha_neg = alpha
    x_pos = np.dot(np.transpose(A_pos),D)
    d_pos = A_pos.sum(axis=0)
    x_pos[d_pos>0] = x_pos[d_pos>0] / d_pos[d_pos>0]
    x_neg = np.dot(np.transpose(A_neg),D)
    d_neg = A_neg.sum(axis=0)
    x_neg[d_neg>0] = x_neg[d_neg>0] / d_neg[d_neg>0]
    return (p_c + alpha_pos*x_pos)*D + (1-p_r-alpha_neg*(1-x_neg))*(1-D)

### Version 6/1/25
def environmental_update(D,A,p_c=0.4,p_r=0.2,q_c=None,q_r=None,alpha=1):
    """ 
    According to the environmental dynamics returns a vector of 
        probabilities to be degraded in the next time step, given
        previous values D, the adjacency matrix of interactions A,
        probabilities to collapse and recover p_c,p_r
    Inputs:
        D - array-like of length N
        A - adjacency matrix NxN; a_ij is i-> j
        p_c,p_r - probability to collapse and recover, resp., [0,1]"""
    x_pos = np.dot(np.transpose(np.where(A>=0,A,0)),D)
    x_neg = np.dot(np.transpose(np.where(A<0,-A,0)),D)
    p_ceff = np.clip(p_c*(1+x_pos)/(1+x_neg),a_min=0,a_max=1) # clipped to [0,1]
    # p_reff = np.clip(p_r*(1+x_neg),a_min=0,a_max=1) # clipped to [0,1]
    return p_ceff + (1-p_ceff-p_r)*D

### Functions for later use in pycrld
def bit_vector(x,width=5):
    """ 
    For integer x, returns a np array of length 'width' 
        of bits in the binary form of x"""
    return np.array([int(xb) for xb in format(x,'0'+str(width)+'b')])

def transition_matrix(p_c,p_r,A):
    """ 
    The `global' state transition matrix T for n nodes 
        where each state 
        is represented essentially as a bit string of
        whether a node is degraded (1) or prosperous (0).
    Output:
        T - a 2^n by 2^n stochastic matrix. 
            The order is given by the binary of the index,
            eg index 7 in a n=3 network corresponds to
            a state of 011, i.e. node 0 prosperous, 
            nodes 1 & 2 degraded
    """
    n = A.shape[0]
    if type(p_c) == float:
        p_c = np.array([p_c]*n)
    if type(p_r) == float:
        p_r = np.array([p_r]*n)
    N = 2**n
    T = np.ones([N]*2)
    for source_state in range(N):
        source_bit = bit_vector(source_state,width=n)
        x_pos = np.dot(np.transpose(np.where(A>=0,A,0)),source_bit)
        x_neg = np.dot(np.transpose(np.where(A<0,-A,0)),source_bit)
        p_ceff = np.clip(p_c*(1+x_pos)/(1+x_neg),a_min=0,a_max=1)
        for target_state in range(N):
            target_bit = bit_vector(target_state,width=n)
            T[source_state,target_state] *= np.prod(1-p_ceff[(1-source_bit&1-target_bit).astype(np.bool)]) # 0 -> 0
            T[source_state,target_state] *= np.prod(p_ceff[(1-source_bit&target_bit).astype(np.bool)]) # 0 -> 1
            T[source_state,target_state] *= np.prod(p_r[(source_bit&1-target_bit).astype(np.bool)]) # 1 -> 0
            T[source_state,target_state] *= np.prod(1-p_r[(source_bit&target_bit).astype(np.bool)]) # 1 -> 1
    return T

def sys_stat_dist_from_Tmatrix(T):
    eigen = np.linalg.eig(np.transpose(T))
    evector = eigen.eigenvectors[:,np.argmin(np.abs(np.abs(eigen.eigenvalues)-1))]
    if len(np.unique(np.sign(evector))) == 1:
        return np.abs(evector)
    else:
        return evector

def stat_dist_T(p_c,p_r,A):
    T = transition_matrix(p_c,p_r,A)
    sys_stat = np.real(sys_stat_dist_from_Tmatrix(T))
    sys_stat = sys_stat/sys_stat.sum()
    stat_dist = np.zeros(A.shape[0])
    for i in range(len(sys_stat)):
        bit_vs = bit_vector(i,width=A.shape[0])
        stat_dist[[bool(b) for b in bit_vs]] += sys_stat[i]
    return stat_dist

def stat_dist_theory_res(D,*data):
    A,q_c,q_r = data
    x_pos = np.dot(np.transpose(np.where(A>=0,A,0)),D)
    x_neg = np.dot(np.transpose(np.where(A<0,-A,0)),D)
    # out = np.zeros(len(D))
    # for i in range(len(D)):
    #     out[i] = D[i] - (1+x_pos[i])*q_c[i]/((1+x_pos)*q_c[i] + (1+x_neg)*q_r[i])
    out = D - (1+x_pos)*q_c/((1+x_pos)*q_c + (1+x_neg)*q_r)
    return out

def stationary_dist_from_graph(G,q_c,q_r,nodelist=None):
    """ 
    Requires:
        - G: a directed acrylic graph   
        - q_c: individual collapse leverages, either float or np.array of size N
        - q_r: individual recovery leverages, either float or np.array of size N
    """
    if not nx.is_directed_acyclic_graph(G):
        raise TypeError("only DAGs allowed")
    if nodelist is None:
        nodelist = np.array(G.nodes)
    N = len(nodelist)
    D = np.zeros(N)
    node_to_ind = {v:ind for ind,v in enumerate(nodelist)}
    
    if type(q_c) == float:
        q_c = np.array([q_c]*N)
    if type(q_r) == float:
        q_r = np.array([q_r]*N)
        
    A = nx.to_numpy_array(G,nodelist=nodelist)
    A_pos = np.where(A>=0,A,0)
    A_neg = np.where(A<0,-A,0)
    k_in = A.sum(axis=0)
    source_mask = k_in==0
    S = list(nodelist[source_mask])
    D[source_mask] = (q_c/(q_c+q_r))[source_mask]
    while len(S) < len(nodelist):
        V_ns = nodelist[~source_mask]
        G_ns = G.subgraph(V_ns)
        for i in V_ns:
            if G_ns.in_degree(i) == 0:
                ind = node_to_ind[i]
                x_pos = np.sum(A_pos[:,ind]*D[:])
                x_neg = np.sum(A_neg[:,ind]*D[:])
                D[ind] = (1+x_pos)*q_c[ind]/((1+x_pos)*q_c[ind] + (1+x_neg)*q_c[ind])
                S += [i]
        source_mask = np.array([v in S for v in nodelist])
    return D

### Utility functions for plotting
def lighten_color(color, amount=0.5):
    """
    Lightens the given color by multiplying (1-luminosity) by the given amount.
    Input can be matplotlib color string, hex string, or RGB tuple. 
    If amount > 1 then the colour is darkened instead

    Examples:
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

### Utility for tipping elements (TE)
def logistic(T,T0,beta):
    return 1/(1+np.exp(beta*(T0-T)))

def piecewise_temp_response(T,T0,T_min,T_max):
    theta = np.zeros_like(T)
    mask_b = (T_min<=T)&(T<T0)
    mask_c = (T0<=T)&(T<T_max)
    theta[mask_b] = ((T-T_min)/(T0-T_min)/2)[mask_b]
    theta[mask_c] = ((T+T_max-2*T0)/(T_max-T0)/2)[mask_c]
    theta[T>=T_max]=1
    return theta

def generate_TE_graph_base(df,atts=['Tipping Element','Threshold (°C) - Est.','Timescale (y) - Est.','Sigmoid Beta']):
    # Inputs:
        # df - Pandas dataframe
        # atts - list of attributes, must be valid columns in df
    # Outputs:
        # G - an empty (edgeless) nx.DiGraph with nodes whose names are from df.index
    G = nx.DiGraph()
    G.add_nodes_from(df.index.values.tolist())
    node_attr = df[atts].to_dict('index')
    nx.set_node_attributes(G, node_attr)
    return G

def generate_TE_graph(df_int, weights, df_ind,
                      edge_atts=['Dist. Min.', 'Dist. Max.'], 
                      atts=['Tipping Element','Threshold (°C) - Est.','Timescale (y) - Est.','Sigmoid Beta']):
    # Inputs:
        # weights - list or array of weights
        # df_ind - Pandas dataframe to generate the base graph
        # df_int - Pandas dataframe of the TE-TE edgelists. 
    G = generate_TE_graph_base(df_ind)
    a = df_int[edge_atts[0]].to_numpy()
    b = df_int[edge_atts[1]].to_numpy()
    edge_weights = a + (b-a)*weights
    G.add_edges_from([(df_int.iloc[i]['Source'],
                       df_int.iloc[i]['Target'], 
                       {'weight':edge_weights[i]}) for i in range(len(df_int))])
    return G