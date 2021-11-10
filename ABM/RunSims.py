# file_name.py
# Python 3.6
"""
Author: Jacob Marglous
Created: Thu Feb 25 13:17:34 2021
Modified: Thu Feb 25 13:17:34 2021

This is a script I wrote for convenience to make running the simulation easier. Run it in a terminal and it will walk you through the process of inputting the selection coefficents. The selection coefficients are a si.txt table generated from the other scripts. The default parameters are d = .05, b0 = 2, K - 5e5, 20,000 timesteps. Inputting a table of calculated equilibrium frequencies is necessary because the ABM is populated with genotype frequencies based on the initial equilibrium distribution, so that the simulation and equilibrium distributions can be compared. 
"""
import numpy as np
import Simulation_Draft2 as st

siLoc = input("Location of selection coefficients to simulate on: ")
si = np.loadtxt(siLoc, delimiter = ",")

XeqLoc = input("location of equilibrium frequencies to compare to: ")
Xeq = np.loadtxt(XeqLoc, delimiter = ",")
savefile = input("Desired location of output file")

test_pop32 = st.population(current_genotype_pops = Xeq[0, :], fitnesses = si, mut_matrix = st.mutdat32, deathrate = .05, b0 = 2, K = 5e5, genotypes = 32)

test_pop32.run(20000, si)
np.savetxt(savefile, test_pop32.genotype_pops, delimiter = ",")    

# parameters = [("Ti", start_T), ("Tf", end_T), ("C_i", start_drug), ("C_f", end_drug), ("T slope", slope_T), ("C slope", slope_drug), ("temp ramp time", ramptime_T), ("drug ramp time", ramptime_drug)]
    
# parameters_array = np.array(parameters, dtype = str)
    
# save_parameters = input("Enter desired filename for parameter values: ")

# np.savetxt(save_parameters, parameters_array, fmt='%s')

