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