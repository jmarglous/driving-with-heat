import numpy as np
from tqdm import tqdm
import matplotlib.pyplot as plt

si = np.loadtxt("C:/Users/jmarg/OneDrive/Documents/DrivingWithHeat/TempRamp2/si.txt", delimiter = ",")
mm32 = np.loadtxt("C:/Users/jmarg/OneDrive/Documents/DrivingWithHeat/mm32.csv", delimiter = ",")

Xt = .04/31*np.ones(32)
Xt[15] = .96

def step(Xt, sv):
    muts = np.zeros(Xt.shape[0])
    selection = np.zeros(Xt.shape[0])
    deltaX = np.zeros(Xt.shape[0])
    for allele in range(Xt.shape[0]):
        muts[allele] = np.sum(mm32[allele, :]*Xt)
        #print("Xt: "+str(Xt))
        #print("muts: "+str(muts))
        selection[allele] = -Xt[allele]*np.sum(Xt*sv)+Xt[allele]*sv[allele]
        #print("selection: "+str(selection))
        deltaX = muts+selection
    
    #print("Xt: "+str(Xt))
    #print("muts: "+str(muts))    
    #print("selection: "+str(selection))

    return Xt+deltaX

def step2(Xt, sv):
    xg = Xt[:-1]
    muts = np.zeros(Xt.shape[0]-1)
    selection = np.zeros(Xt.shape[0]-1)
    deltaX = np.zeros(Xt.shape[0]-1)
    for allele in range(muts.shape[0]):
        muts[allele] = np.sum(mm32[allele, :]*xg)
        #print("Xt: "+str(Xt))
        #print("muts: "+str(muts))
        selection[allele] = -xg[allele]*np.sum(xg*sv[:-1])+xg[allele]*sv[allele]
        #print("selection: "+str(selection))
        deltaX = muts+selection
    
    print("Xt: "+str(Xt))
    print("muts: "+str(muts))    
    print("selection: "+str(selection))

    return Xt+deltaX

def path(Xt, sv, steps, track = False):
    Xts = []
    
    if track == False:
        for i in range(steps):
            Xt = step(Xt, sv)
            Xts.append(Xt)
            
    if track == True:
        for i in tqdm(range(steps), desc = "path"):
            Xt = step(Xt, sv)
            Xts.append(Xt)
            
    Xts = np.array(Xts)
    return Xts

def approxEquil(Xt, si, stepsPerTime):
    
    approx_Xeqs = []
    for i in tqdm(np.arange(0, si.shape[0], 1)):
        
        path_timestep = path(Xt, si[i], stepsPerTime)
        Xeq_approx = path_timestep[-1, :]
        
        extraCounter = 0
        while np.min(Xeq_approx) < 0:
            path_timestep = path(Xt = Xeq_approx, sv = si[i], steps = 10000, track = True)
            Xeq_approx = path_timestep[-1, :]
            extraCounter += 1
            
            if extraCounter == 100:
                print(i)
                break
        approx_Xeqs.append(Xeq_approx)
        
        Xt = Xeq_approx

    approx_Xeqs = np.array(approx_Xeqs)
    return approx_Xeqs

def plot(Xeqs):
    for i in range(32):
        plt.plot(Xeqs[:, i])
        plt.yscale("log")
        
def checkZeros(Xeqs, stop = True):
    NegTimes = []
    for i in range(Xeqs.shape[0]):
        if np.min(Xeqs[i, :]) < 0:
            print(i)
            NegTimes.append(i)
            if stop == True:
                break
    NegTimes = np.array(NegTimes)
    return NegTimes