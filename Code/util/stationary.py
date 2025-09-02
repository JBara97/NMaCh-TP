import numpy as np
import scipy as sp
### Stationary distribution
def stat_dist_theory_res(D,*data):
    A,q_c,q_r = data
    x_pos = np.dot(np.transpose(np.where(A>=0,A,0)),D)
    x_neg = np.dot(np.transpose(np.where(A<0,-A,0)),D)
    # out = np.zeros(len(D))
    # for i in range(len(D)):
    #     out[i] = D[i] - (1+x_pos[i])*q_c[i]/((1+x_pos)*q_c[i] + (1+x_neg)*q_r[i])
    out = D - (1+x_pos)*q_c/((1+x_pos)*q_c + (1+x_neg)*q_r)
    return out

def stat_dist_theory(A,q_c,q_r,T,T_lim,beta):
    return sp.optimize.fsolve(stat_dist_theory_res, np.ones(9), args=(A,temp_response(T, T_lim,beta)*q_c,q_r))