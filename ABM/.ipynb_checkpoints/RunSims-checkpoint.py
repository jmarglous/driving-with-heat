# file_name.py
# Python 3.6
"""
Author: Jacob Marglous
Created: Thu Feb 25 13:17:34 2021
Modified: Thu Feb 25 13:17:34 2021

Description
"""
import numpy as np
import Simulation_Draft2 as st

for i in range (10):
    st.test_pop32.reset(15)
    st.test_pop32.run(20000, st.s32_flat)
    np.savetxt("32GenotypesFlat_start15_" +str(i)+ ".csv", st.test_pop32.genotype_pops, delimiter = ",")    



