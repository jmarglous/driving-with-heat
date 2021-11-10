from __future__ import absolute_import, division, print_function, unicode_literals

import matplotlib.pyplot as plt
import numpy as np
import scipy.interpolate as inp
from tqdm import tqdm

np.seterr(all='raise');

import sys
np.set_printoptions(threshold=sys.maxsize)

import interpolationFunction as InpFunc 

interpolation_functions = []
for i in range(32):
    interpolation_functions.append(InpFunc.interpolation_func(i))
interpolation_functions = np.array(interpolation_functions)

mutdat32_fromscript = np.loadtxt("mutdat32_fromscript.csv", delimiter = ",")

mm32 = np.transpose(mutdat32_fromscript)

for ii in range(mutdat32_fromscript.shape[0]):
    mm32[ii,ii]=-(np.sum(mm32[ii])-mm32[ii,ii])

    
change_T = input("Do you want to change temperature? Y/N: ")
if change_T == "Y" or change_T == "y":
    start_T = float(input("Starting temperature? (enter w/in 20-41; default is 27): "))
    end_T = float(input("Ending temperature? (enter w/in 20-41; default is 37): "))
    slope_T = float(input("Rate of temperature change? (must be positive; default is .002): "))
    ramptime_T = float(input("Time of temperature change? (Default is 10110): "))

if change_T == "N" or change_T == "n":
    start_T = float(input("What would you like the (constant) temperature to be? (Default is 27): "))
    end_T = start_T
    slope_T = 0
    ramptime_T = 0
    
change_drug = input("Do you want to change drug concentration? Y/N: ")
if change_drug == "Y" or change_drug == "y":
   start_drug = float(input("Starting drug concentration? (default is 1e-9): "))
   end_drug = float(input("Ending drug concentration? (default is 1e2): "))
   slope_drug = float(input("Rate of drug concentration change? (must be positive; default is .002): "))
   ramptime_drug = float(input("Time of drug concentration change? (Default is 10110): "))


    
if change_drug == "N" or change_drug == "n":
    start_drug = float(input("What would you like the (constant) drug concentration to be? (Default is 1e2): "))
    end_drug = start_drug
    slope_drug = 0
    ramptime_drug = 0
    
def tt(t):
    if change_T == "Y" or change_T == "y":
        ttt = (end_T-start_T)/(1. + np.exp(-slope_T*(t - ramptime_T)))+start_T
    else:
        ttt = start_T
    return ttt

def cc(t) :
        
    if change_drug == "Y" or change_drug == "y":
        cct=(end_drug-start_drug)/(1. + np.exp(-slope_drug*(t - ramptime_drug)))+start_drug
    else:
        cct = start_drug
    return cct
    
def dosage(times):
    dosage_schedule = []
    for i in times:
        dosage_schedule.append((tt(i), cc(i)))
    dosage_schedule = np.array(dosage_schedule)
    return dosage_schedule

selection_time_start = float(input("At what time would you like to start? Default is 0: "))
selection_time_end = float(input("At what time would you like to end? Defult is 45000: "))
selection_time_interval = float(input("How frequently would you like to pull fitnesses from the landscape? Default is 1 time step: "))
    
selection_times = np.arange(selection_time_start, selection_time_end, selection_time_interval)

plot_TD = input("Do you want to plot the temperature-drug dosage plot? Y/N: ")
if (plot_TD == "Y") or (plot_TD =="y"):
    fig = plt.figure()
        
    dosage_schedule = dosage(selection_times[::10])
    ax1 = fig.add_subplot(221)
    ax1.plot(selection_times[::10], dosage_schedule[:, 0])
    ax1.set_xlabel("time t ")
    ax1.set_ylabel("temperature tt[t]")
        
    ax2 = fig.add_subplot(222)
    ax2.plot(selection_times[::10], dosage_schedule[:, 1])
    ax2.set_xlabel("time t ")
    ax2.set_ylabel("concentration cc[t]")
        
    ax3 = fig.add_subplot(212)
    ax3.scatter(dosage_schedule[::10, 1], dosage_schedule[::10, 0], cmap = "viridis", c = selection_times[::100])
    ax3.set_ylabel("temperature ($\degree$ C)")
    ax3.set_xlabel("drug concentration")
    ax3.set_xscale("log")
    ax3.set_xlim([1e-10, 1e10])
    fig.tight_layout()
        
    savefig_TD = input("Do you want to save the temperature-drug dosage plot? Y/N: ")
    if (savefig_TD == "Y") or (savefig_TD == "y"):
        TDplot_name = input("Enter desired temperature-drug dosage plot filename: ")
        fig.savefig(TDplot_name, dpi = 300)

def sftc(temp, conc):
    ut = np.zeros(32)
        
    ut31 = interpolation_functions[31].__call__(temp, conc)[0]
        
    for i in range(32):
        growthrate_i = interpolation_functions[i].__call__(temp, conc)[0]
        
        if growthrate_i < 0:
            growthrate_i = 0
            
        try:
            ut[i] = growthrate_i/ut31 - 1
        except FloatingPointError:
            ut[i] = 100000
    
    return ut
    
print("Calculating selection coefficients: ")
st=[]
for i in tqdm(selection_times):
    st.append((sftc(tt(i), cc(i))))
    
st = np.asarray(st) 

dst=np.diff(st, axis = 0) 
    
si=[];dsi=[]
for i in tqdm(range(32)):
    ysi=st[:,i] 
    ydsi=dst[:,i]
    tck_si = inp.splrep(selection_times, ysi,k=3,s=0) 
    tck_dsi = inp.splrep(selection_times[1:,], ydsi,k=3,s=0) 
    
    ysin = inp.splev(selection_times, tck_si, der=0)
    ydsin = inp.splev(selection_times[1:,], tck_dsi, der=0)#
    si.append(ysin)
    dsi.append(ydsin)

si=np.transpose(np.asarray(si)) #putting time row wise
dsi=np.transpose(np.asarray(dsi))
    
plot_s = input("Would you like to plot the selection coefficients? Y/N: ")
if (plot_s == "Y") or (plot_s == "y"):
    fig2 = plt.figure(figsize = [8, 8])
    ax21 = fig2.add_subplot(221)
    ax22 = fig2.add_subplot(222)
        #ax23 = fig2.add_subplot(212)
    for i in range(32):
        ax21.plot(selection_times[::10], si[::10, i])
        ax22.plot(selection_times[1::10], dsi[::10, i])
    
    ax21.set_xlabel("Time")
    ax22.set_xlabel("Time")
    ax21.set_ylabel("Selection coefficients")
    ax22.set_ylabel("ds/dt")
        
    fig2.tight_layout()

        #ax21.set_xlim(np.min((ramptime_drug, ramptime_T))-2000, np.max((ramptime_drug, ramptime_T))+2000)
        #ax22.set_xlim(np.min((ramptime_drug, ramptime_T))-2000, np.max((ramptime_drug, ramptime_T))+2000)
        
        
    savefig_s = input("Would you like to save the plot of selection coefficients? Y/N: ")
    if (savefig_s == "Y") or (savefig_s == "y"):
        splot_name = input("enter desired filename of selection coefficient plot: ")
        fig2.savefig(splot_name, dpi = 300)
save_s = input("Would you like to save the table of selection coefficients? Y/N: ")    
if save_s == "Y" or "y":
    s_name = input("Enter desired filename for table of selection coefficients: ")
    np.savetxt(s_name, st, delimiter = ",")
parameters = [("Ti", start_T), ("Tf", end_T), ("C_i", start_drug), ("C_f", end_drug), ("T slope", slope_T), ("C slope", slope_drug), ("temp ramp time", ramptime_T), ("drug ramp time", ramptime_drug)]
    
parameters_array = np.array(parameters, dtype = str)
    
save_parameters = input("Enter desired filename for parameter values: ")

np.savetxt(save_parameters, parameters_array, fmt='%s')

