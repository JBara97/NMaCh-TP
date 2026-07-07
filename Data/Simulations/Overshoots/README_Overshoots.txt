Contains simulation data sufficient to reproduce figures included in the manuscript.

Note that due to data size 
- (500 interaction matrices X 9 tipping elements X 451 time steps X 2237 PROVIDE ensemble members) -
simulation data for the ensemble dynamics must be run at the point of use.

For a given scenario, s, structure for Data/Simulations/Overshoots/s/

OAT/
	# Contains data for the one-at-a-time (OAT) sensitivity analysis
	# Namely, contains the 2300 tipping (ie aggregate) risk EXCLUDING the reefs
	# 	and the 2300 individual tipping risk for the coral reefs
	# Parameter ranges are as follows:
	#	temp_init: [0,0.1,...,5] (°C)
	#	temp_peak: [0,0.1,...,5] (°C)
	#	temp_stab: [0,0.1,...,5] (°C)
	#	time_up:   [0,10,..,200] (yr)
	#	time_peak: [0,10,..,200] (yr)
	#	time_down: [0,10,..,200] (yr)
TAT/
	# Contains data for the two-at-a-time (TAT) sensitivity analysis
	# Namely, contains the 2300 tipping (ie aggregate) risk EXCLUDING the reefs
	# 	and the 2300 individual tipping risk for the coral reefs
	# Parameter ranges are as follows:
	#	temp_peak (°C), temp_stab (°C): [0,0.1,...,5] X [0,0.1,...,3]
	#	temp_stab (°C), time_down (yr): [0,0.1,...,3] X [0,5,...,200]
