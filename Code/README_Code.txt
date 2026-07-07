Folder structure for Code/

Notebooks/
	# Jupyter notebooks used for development, simulations, analysis and visualisations across a number of projects.
	ESD/
		# Code for Bara et al. (ESD, 2026) https://doi.org/10.5194/esd-17-333-2026
	NCC Submission/
		# Code for overshoots project
	OLD/
		# Early dev code, can be ignored.
scripts/
	# useful (though old) scripts used particularly for the 
	# Bara et al. (ESD, 2026) https://doi.org/10.5194/esd-17-333-2026 project
util/
	# Python package for the whole project. Contains multiple modules namely:
	dynamics.py
		# Code to simulate the interacting tipping dynamics
	stationary.py
		# Code for equilibrium/stationary behaviour. Mostly used in ESD paper
	plot.py
		# Helper/wrapper code for easier plotting
	overshoots.py
		# Particularly useful for the overshoots project. 
		# Contains eg the stylised overshoot model, and other useful wrappers
	
