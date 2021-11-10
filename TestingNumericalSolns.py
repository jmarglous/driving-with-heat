import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
import scipy.optimize as op
import DeterministicPrediction as det 

mutdat32_fromscript = np.loadtxt("mutdat32_fromscript.csv", delimiter = ",")
mm32 = np.transpose(mutdat32_fromscript)

#siloc = input("Location of set of selection coefficients to solve with: ")
#si = np.loadtxt(siloc, delimiter = ",")

for ii in range(mutdat32_fromscript.shape[0]):
    mm32[ii,ii]=-(np.sum(mm32[ii])-mm32[ii,ii])

def g(xv) :
    x=31;
    mut=np.zeros((x,x))

    for ii in range(x):
        for jj in range(x):
            if ii!=jj:
                mut[ii,jj]=-xv[ii]*xv[jj]
            else:
                mut[ii,jj]=xv[ii]*(1.0 - xv[ii])
    return mut

def losseq(xg, sv):
    #x = np.append(xg, 1.0-np.sum(xg))

    vec_eq = np.zeros(xg.shape[0])
    for i in range(vec_eq.shape[0]):
        vec_eq[i] = np.sum(mm32[i, :]*xg)+xg[i]*sv[i]
        sum_xg_isv_i = np.sum(xg*sv)
        vec_eq[i] = vec_eq[i] - xg[i]*sum_xg_isv_i

    return vec_eq

def losseq2(xg, sv):  
    '''
    Original function from Shamreen
    '''
    sM=np.diag(sv)
    DM=np.identity(32)
    x = np.append(xg, 1.0-np.sum(xg))
    sv31 = sv[0:31] # Leaving out last element for dimensional compatibility
    T1 = np.dot(np.dot(mm32,(DM+sM)),x)[0:31]
    T2 = np.dot(g(xg),sv31)
    return (T1 + T2)

def catch_negs(Xt):
    Xt_min = np.min(Xt)
    if Xt_min < 0:
        return True
    else:
        return False


def solveSequence(xeq0, si):
    
    xg = xeq0
    Xeq = []
    for i in tqdm(range(si.shape[0])): #tqdm here  

        #if i <= 4/9*si.shape[0]:
        Xt = op.fsolve(losseq, xg, args = si[i])
            #Xt = np.append(Xt, 1-np.sum(Xt))
        #diffs = np.abs(Xt - xg)
            
            # if np.argmax(Xt) != np.argmax(xg):
            #     print(i)
            #     X_freqs = det.path(xg, si[i], 1000, track = True)
            #     Xt = op.fsolve(losseq, X_freqs[-1, :], args = si[i])
        while catch_negs(Xt) == True: #or np.abs(np.sum(Xt - xg)) > 1e-15:
                #print("Using approx")
            Xsol = Xt
                #Xt = np.append(xg, 1-np.sum(xg)) ##NOT KOSHER
            X_freqs = det.path(Xsol, si[i], 10000, track = False)
            Xsol = op.fsolve(losseq, X_freqs[-1, :], args = si[i])
            xg = X_freqs[-1, :]
            Xt = Xsol
            
        Xeq.append(Xt)
        xg = Xt
           # print(i)
        #else:
           # Xeq.append(Xt)
        
    Xeq = np.array(Xeq)
    return Xeq

def solveSequenceOldLosseq(xeq0, si):
    
    xg = xeq0[:-1]
    Xeq = []
    for i in tqdm(range(si.shape[0])): #tqdm here  

        #if i <= 4/9*si.shape[0]:
        Xt = op.fsolve(losseq2, xg, args = si[i])
        Xt = np.append(Xt, 1-np.sum(Xt))
        
        #diffs = np.abs(Xt - xg)
            
            # if np.argmax(Xt) != np.argmax(xg):
            #     print(i)
            #     X_freqs = det.path(xg, si[i], 1000, track = True)
            #     Xt = op.fsolve(losseq, X_freqs[-1, :], args = si[i])
        while catch_negs(Xt) == True: #or np.abs(np.sum(Xt - xg)) > 1e-15:
            Xsol = Xt   
            #print("Using approx")
                #Xt = np.append(xg, 1-np.sum(xg)) ##NOT KOSHER
            X_freqs = det.path(Xsol, si[i], 10000, track = False)
            Xsol = op.fsolve(losseq2, X_freqs[-1, :-1], args = si[i]) #xt has length 31
            Xsol = np.append(Xsol, 1-np.sum(Xsol))

            xg = X_freqs[-1, :]
            Xt = Xsol
        
        Xeq.append(Xt)
        xg = Xt[:-1]
           # print(i)
        #else:
           # Xeq.append(Xt)
        
    Xeq = np.array(Xeq)
    return Xeq

def solveSequenceNoApprox(xeq0, si):
    
    xg = xeq0
    Xeq = []
    for i in tqdm(range(si.shape[0])):    
        if i <= 4/9*si.shape[0]:
            Xt = op.fsolve(losseq, xg, args = si[i])
            #Xt = np.append(Xt, 1-np.sum(Xt))
            # if catch_negs(Xt) == True:
            #     print("Using approx")

            #     #Xt = np.append(xg, 1-np.sum(xg)) ##NOT KOSHER
            #     X_freqs = det.path(Xeq[-1], si[i], 100000, track = False)
            #     Xt = op.fsolve(losseq, X_freqs[-1, :], args = si[i])
            Xeq.append(Xt)
            xg = Xt
            
        else:
            Xeq.append(Xt)
        
    Xeq = np.array(Xeq)
    return Xeq