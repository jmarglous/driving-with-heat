import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from mpl_toolkits import mplot3d
from scipy import interpolate as inp
#import seaborn as sns


def growthrate_temp(T, Tmin = 3.99, Tmax = 43.7, b3 = .0410, c3 = .161):
    '''
    This function takes an input of a single temperature and outputs a growthrate 
    value at that temperature based on the Ratkowsky 3 growthrate vs. temperature 
    model and the parameters described in Zwiettering et al. 
    '''    
    muM = np.power(b3*(T-Tmin), 2)*(1-np.power(np.e, (c3*(T-Tmax))))
    
    if (muM < 0):
        muM = 0
    return muM

def temp_curve (T, Tmin = 3.99, Tmax = 43.7, b3 = .0410, c3 = .161):
    '''
    This function inputs an array of temperature values and outputs a growthrate
    curve by calculating the growthrate at each temperature value via the 
    Ratkowsky3 function. The output is a nx2 array with n the number of elements 
    in the T aray and each row representing a (temperature, growthrate) coordinate. 
    '''
    growthrates = []
    temps_growths = []
    for i in range(len(T)):
        growthrate = growthrate_temp(T[i], Tmin, Tmax, b3, c3)
        temp_growth = [(T[i], growthrate)]
        temps_growths.append(temp_growth)
        
    temps_growths = np.array(temps_growths)
    temps_growths = temps_growths[:, 0, :] #just slices out the 2nd dimension of the temps_growths array -- somehow I have an extra one 
    return temps_growths

def growthrate_conc(MIC, g_drugless, concentrations, c = -.6824968):
    '''
   Takes as input an MIC, g_drugless, and concentration, and fitting constant, and outputs
   a growthrate at that drug concentration as defined by the logistic curve described by 
   Ogbunufor et al. 
   '''
    
    growthrate = g_drugless/(1+19*np.exp((np.log10(MIC)-concentrations)/c))
    return growthrate


tempData = np.loadtxt("./weinreich_temp_modifiedlabels.csv", delimiter = ",", skiprows=1, usecols = (1, 2, 3, 4, 5, 6))

def MICs_fromAllele(allele):
    '''
    Returns array of MICs at 20, 25, 30, 37, 41 degrees Celsisus for the allele inputted in integer notation (i.e. 0-31)
    '''
    MICs = tempData[allele, :]

    return MICs
    

def plot_growthrate_conc_curves(allele, concentrations, c = -.6824968, save = False, fname = "growthrate_conc_curves.png"):
    '''
    Inputs an array from the of MICs at lots growthrate([drug]) at the experimental 
    temperatures (20, 25, 30, 37, 41). Growthrate([drug]) is calculated 
    from the growthrate_conc function, with g_drugless determined from the 
    growthrate_temp function.
    '''
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    temperatures_exp = [20, 25, 30, 35, 37, 41]
    
    MICs = MICs_fromAllele(allele)

    for i in range(len(temperatures_exp)):
        xs = temperatures_exp[i]*np.ones(len(concentrations))
        ys = np.log10(concentrations)
        zs = growthrate_conc(MIC = MICs[i], g_drugless = growthrate_temp(temperatures_exp[i]), concentrations = concentrations, c = -.6824968)
    
        ax.plot(xs, ys, zs)
    ax.set_xlabel("Temperature ($\degree$C)")        
    ax.set_ylabel("$[drug]$")
    ax.set_zlabel("Maximum specific growth rate")
    

    ax.set_title("Growth rate as function of temperature and drug concentration")
    if (save == True):
        fig.savefig(fname)
        
def MIC_fromT(allele, T):
    '''
    Just outputs the MIC of the inputted allele at the inputted temperature. 
    '''
    if (T == 20):
        return tempData[allele, 0]
    if (T == 25):
        return tempData[allele, 1]
    if (T == 30):
        return tempData[allele, 2]
    if (T == 35):
        return tempData[allele, 3]
    if (T == 37):
        return tempData[allele, 4]
    if (T == 41):
        return tempData[allele, 5]

def growthrate_interpolation(allele, concentrations):

    MICs = MICs_fromAllele(allele)

    temperatures_exp = np.array([20, 25, 30, 35, 37, 41])

    temps_concs_growthrates = []
    for i in range(len(MICs)):
        for j in range(len(concentrations)):
            growthrate = growthrate_conc(MIC = MICs[i], g_drugless = growthrate_temp(T = temperatures_exp[i]), concentrations = concentrations[j])
            temps_concs_growthrates.append((temperatures_exp[i], np.log10(concentrations[j]), growthrate))

    temps_concs_growthrates = np.array(temps_concs_growthrates)
    temps_concs = temps_concs_growthrates[:, :-1]
    growthrates = temps_concs_growthrates[:, -1]

    temps_for_grid = np.arange(20,41, .5)
    concs_for_grid =  np.log10(concentrations) #The log and the power of 10 cancel but I just want to keep track of what's linear and what's log scale
    grid_temps, grid_concs = np.meshgrid(temps_for_grid, concs_for_grid)
    
    grid_values = inp.griddata(temps_concs, growthrates, (grid_temps, grid_concs), method='cubic')
    
    for i in range(np.shape(grid_values)[0]):
        for j in range(np.shape(grid_values)[1]):
            if (grid_values[i, j] < 0):
                grid_values[i, j] = 0
    
    temps_concs_grid = []
    for i in range(np.shape(grid_values)[0]):
        for j in range(np.shape(grid_values)[1]):
            temp_conc_grid = [grid_temps[i, j], grid_concs[i, j], grid_values[i, j]]
            temps_concs_grid.append(temp_conc_grid)
    temps_concs_grid = np.array(temps_concs_grid)
    return temps_concs_grid



def plot_growthrate_T_conc(allele, size, plotType = "heatmap", save = False, fname = "growthrate_T_conc.png"):
    '''
    The output from growthrate_interpolation is many (T, conc, growthrate) coordinates. 
    This function takes those values and plots them as a heatmap. 
    '''
    fig = plt.figure(figsize = size)
    concentrations = np.power(10., np.arange(-10, 10, .05))
    temps_concs_grid = growthrate_interpolation(allele, concentrations)
    
    if (plotType == "heatmap"):
        ax = fig.add_subplot(111)
        #cbar_ax = fig.add_axes([.91, .3, .03, .4])
        
        
        df = pd.DataFrame(temps_concs_grid, columns=["temperature", "log_10 drug concentration", "maximum specific growthrate"])
        df = df.pivot("temperature", "log_10 drug concentration", "maximum specific growthrate")
        heatmap = sns.heatmap(df, cmap = "coolwarm", vmin = 0, vmax = 1.2, ax = ax, cbar = True, xticklabels = 99, yticklabels = 10)
        
    if (save == True):
        fig.savefig(fname)

def plot_mult_growthrate_T_conc(save = False, fname = "growthrate_T_conc_alleles.png"):
    fig = plt.figure(figsize = [12, 12])
    concentrations = np.power(10., np.arange(-10, 10, .05))
    cbar_ax = fig.add_axes([.91, .3, .03, .4])
    
    for i in range(32):
        axi = fig.add_subplot(8, 4, i+1)
        
        do_xticks = False
        if (i >= 28):
            do_xticks = 50
        
        do_yticks = False
        if (i%4 == 0):
            do_yticks = 10
            
        temps_concs_grid = growthrate_interpolation(i, concentrations)
        
        df = pd.DataFrame(temps_concs_grid, columns=["temperature", "log_10 drug concentration", "maximum specific growthrate"])
        df = df.pivot("temperature", "log_10 drug concentration", "maximum specific growthrate")
        heatmap = sns.heatmap(df, cmap = "coolwarm", vmin = 0, vmax = 1.2, ax = axi, cbar = i == 0, xticklabels = do_xticks, yticklabels = do_yticks, cbar_ax=None if i else cbar_ax)
        
        
            
        axi.set_ylabel(" ")
        axi.set_xlabel(" ")
    
    fig.text(0.5, -0.1, 'log10 drug concentrations', ha='center')
    fig.text(0.04, 0.5, 'temperatures', va='center', rotation='vertical')

    if (save == True):
        fig.savefig(fname)
    #fig, axs = plt.subplots(8, 4)
    
    #for k in range(32):
        #temps_concs_grid = growthrate_interpolation2(allele = k, concentrations = concentrations)
        #df = pd.DataFrame(temps_concs_grid, columns=["temperature", "log_10 drug concentration", "maximum specific growthrate"])
        #df = df.pivot("temperature", "log_10 drug concentration", "maximum specific growthrate")
        #heatmap = sns.heatmap(df, cmap = "coolwarm", ax = axs[i])
    
