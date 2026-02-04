import numpy as np

def TIPMIP(T_peak=2.5, T_conv=0, t_peak=200, t_end = 1000):
    """ 
    Stylised TIPMIP overshoot scenarios (https://egusphere.copernicus.org/preprints/2025/egusphere-2025-3604/)
    are characterised by an initial ramp up equivalent to +0.2C/decade from 0C 
    up to some peak temperature, T_peak, at which temperatures stay constant for some 
    peak duration, t_peak, (for example this might be caused by net-zero emissions) 
    followed by a ramp down equivalent of -0.2C/decade down towards the convergent temperature, 
    T_conv. 
    t_end controls how many time-steps to include in the output

    Parameters
    ---
    T_peak : int or float. default 2.5C
        The temperature at the peak
    t_peak : int. default 200 years
        The length of time spent precisely at temperature T_peak (NOT the return time to T_conv)
    T_conv : int or float. default 0C
        The convergent temperature, i.e. after peaking the temperature at which the scenario will 
        converge to and remain at at long times. T_conv >= T_peak
    t_end : int. default 5000 years
        The total number of time-steps the output scenario will have
        
    Returns
    ---
    temp : numpy.array
        The output will be a numpy.array of length t_end

    Examples
    ---
    The usual example is star from 0C -> 2C -> 0 and stay at 2C for 100 years
    >>> TIPMIP(T_peak = 2, T_conv = 0, t_peak = 100)
    
    With the alternative going up to 4C instead, (but still starting and ending at 0)
    >>> TIPMIP(T_peak = 4, t_peak = 100)

    Alternatively following the Overshoots Conference, more people are interested in 
    overshooting 1.5C (and thus returning to 1.5C)
    >>> TIPMIP(T_peak = 2.3, T_conv = 1.5)
    """
    # hard coded speed of ramp up and ramp down as +- 0.2C/decade
    # also hard coded each time step as 1 year
    if T_peak < T_conv:
        raise ValueError("T_peak is smaller than T_conv, this would be an undershoot not overshoot.")

    temp = np.concatenate((np.arange(0,T_peak+.01,step=0.02), # ramp up from 0 -> T_peak at a rate of +0.2C/dec
                           np.ones(t_peak)*T_peak, # remain at peak for t_peak years
                           np.arange(T_conv,T_peak,step=0.02)[::-1])) # ramp down from T_peak -> T_conv at a rate of -0.2C/dec
    if t_end <= len(temp):
        return temp[:t_end]                          
    else:
        return np.concatenate((temp, np.ones(t_end-len(temp))*T_conv)) # remain at T_conv for t_at_conv years

def moving_avg(temps,N=20):
    """
    Simple implementation of a moving-average (mean) over a window of N time-steps.
    Specifically for the moving_avg of temps at time-step $i >= N-1$ is given by:
        moving_avg(i) = (temps(i-N+1) + temps(i-N+2) + ... + temps(i))/N
    For the edge cases i < N-1:
        moving_avg(i) = (temps(0) + temps(1) + ... + temps(i))/(i+1)
    
    Reasoning for the edge case is since we are primarily concerned about the 20-yr average for 
     global warming temperatures for Paris Agreement compliance, then it doesn't make sense to 
     include future temperatures. e.g. the 20-yr average at 2025 should not include 2026- temps

    Parameters (required)
    ---
    temps : numpy.array
        The temperature trajectory
    
    Parameters (optional)
    N : int
        The size of the moving window in time-steps. Default 20.
    
    Returns
    ---
    moving_avg : numpy.array
        The output will be a numpy.array of length len(temps)
    """
    moving_avg = np.cumsum(temps)
    moving_avg[N:] -= moving_avg[:-N]
    moving_avg[N:] /= N
    moving_avg[:N] /= np.arange(1,N+1)
    return moving_avg

def PA_compliant(temps):
    """ 
    Returns whether a given temperature trajectory is compliant 
     with the Paris Agreement:
    "Holding the increase in the global average temperature to well below 2°C 
        above pre-industrial levels and to pursue efforts to limit the 
        temperature increase to 1.5°C above pre-industrial levels..."
    Here, we have interpreted this to mean two conditions must be met:
    1) the maximum 20-yr avg value is <= 2C
    2) The final (long term) 20-yr avg is <= 1.5C

    Parameters
    ---
    temps : numpy.array
        The temperature trajectory, must have len(temps) >= 20

    Returns
    ---
    out : Bool
    """
    if np.max(temps) > 2:
        return False
    elif np.mean(temps[-20:]) <= 1.5:
        return True
    else:
        return False

def OS_ref(temps,ref=1.5,window=1,tol=0.05):
    """ 
    Returns whether a given temperature trajectory overshoots a 
     reference temperature, temp. In particular we do this by checking
     whether the time between consecutive crossings (if any) of the ref 
     temperature (that is being withing the tolerance tol of ref) is bigger 
     than a reference duration, window

    Parameters (required)
    ---
    temps : numpy.array
        The temperature trajectory, must have len(temps) >= 20

    Parameters (optional)
    ---
    ref : float
        The reference temperature w.r.t which the overshoot occurs. Default 1.5C
    window : float
        How far apart in time two consecutive crossing of ref must occur to count as 
        two separate instances rather than, say, being a contiguous period during which 
        the temperature hovers around ref. Default 1.
    tol : float
        The tolerance level at which to accept a temperature is near enough the ref. 
        Default 0.05C

    Returns
    ---
    out : Bool
    """
    stat_yrs = np.where(np.abs(df[s].to_numpy()-ref)<=tol)[0]
    if len(stat_yrs) < 2:
        return False
    elif np.diff(stat_yrs).max()>window:
        return True
    else:
        return False 
    return 

def OS_model(t, *args):
    """ 
    Stylised overshoot scenario, similar to TIPMIP, but parametrised by 6-parameters.

    Parameters
    ---
    t : np.array
        Vector of time-points
    T_0 : int or float
        Initial temperature
    T_peak : int or float
        The temperature at the peak
    T_os : int or float
        The temperature to return to after the overshoot
    t_up : int or float
        The temporal duration of the ramp-up
    t_plat : int or float
        The temporal duration of the peak
    t_down : int or float
        The temporal duration of the ramp down
        
    Returns
    ---
    temp : numpy.array
        The output will be a numpy.array of the same shape as t
    """
    T_0,T_peak,T_os,t_up,t_plat,t_down=args
    out = T_0 + t[t<=t_up]*(T_peak-T_0)/t_up
    out = np.concatenate([out,[T_peak]*len(t[(t_up<t)&(t<=t_up+t_plat)])])
    out = np.concatenate([out, T_peak - (t[(t_up+t_plat<t)&(t<=t_down+t_up+t_plat)]-t_up-t_plat)*(T_peak-T_os)/(t_down)])
    return np.concatenate([out, [T_os]*(len(t)-len(out))])

def local_safety(D,safety_threshold=0.1,
              dimension_kw=dict(year=-1,sample=1)):
    return np.max(D,axis=dimension_kw['year'],keepdims=True).mean(axis=dimension_kw['sample'],keepdims=True)<=safety_threshold
    
def global_safety(D,safety_threshold=0.1,
                  dimension_kw=dict(year=-1,sample=1,tipping_element=2)):
    return local_safety(D,safety_threshold=0.1,
                  dimension_kw=dimension_kw).mean(axis=tipping_element,keepdims=True)

def critical_parameter(D,pars,axis=0,**kwargs):
    return np.max(np.squeeze(local_safety(D,**kwargs))*pars.reshape(len(pars),1),
                  axis=axis)