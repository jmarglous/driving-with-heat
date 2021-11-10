import numpy as np
from scipy import interpolate as interp
import tempConcAnalysis as tca
import matplotlib.pyplot as plt

def growthrate_tempConc(allele, exp_temp, conc):
    MIC = tca.MIC_fromT(allele, exp_temp)
    g_drugless = tca.growthrate_temp(T = exp_temp)
    growthrate = tca.growthrate_conc(MIC, g_drugless = g_drugless, concentrations = conc)
    
    return growthrate

def interpolation_func(allele):
    temps_exp = np.array([20, 25, 30, 35, 37, 41])
    temps_exp = np.array(temps_exp)
    concentrations = np.arange(-10, 10, .01, dtype = float)
    CC, TT = np.meshgrid(concentrations, temps_exp)
    
    
    growthrates = np.zeros((CC.shape[0], TT.shape[1]))
    
    for i in range(growthrates.shape[0]):
        for j in range(growthrates.shape[1]):
            growthrates[i, j] = growthrate_tempConc(allele, TT[i, j], CC[i, j])
    interpFunc = interp.RectBivariateSpline(temps_exp, concentrations, growthrates, s = 0, kx = 3, ky =3)
    
    return interpFunc

def interpolation_func_alt(allele):
    temps_exp = np.array([20, 25, 30, 35, 37, 41])
    temps_exp = np.array(temps_exp)
    concentrations = np.arange(-10, 10, .01, dtype = float)
    
    temps_concs_growthrates = []
    for conc in concentrations:
        for temp in temps_exp:
            growthrate = growthrate_tempConc(allele, temp, conc)
            temps_concs_growthrates.append((temp, conc, growthrate))
    
    temps_concs_growthrates = np.array(temps_concs_growthrates)
    
    temps = temps_concs_growthrates[:, 0]
    concs = temps_concs_growthrates[:, 1]
    growthrates = temps_concs_growthrates[:, 2]
    
    interpFunc = interp.SmoothBivariateSpline(temps, concs, growthrates, s = 0)
    return interpFunc