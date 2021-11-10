# -*- coding: utf-8 -*-
"""
Created on Fri Jul 16 10:42:04 2021

@author: jmarg
"""
import numpy as np
import matplotlib.pyplot as plt

#Xeq = np.loadtxt("Xeq.txt", delimiter = ",")

def plot(Xeq):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    
    for i in range(32):
        ax.plot(Xeq[:, i])
    
    ax.set_xlabel("Time")
    ax.set_ylabel("Allele frequency")
    ax.set_yscale("log")
    title = input("Plot title: ")
    ax.set_title(title)

    #fig.savefig("Equilibriumfreqs.png", delimiter = ",")    

def plotDiffs(Xeq):    
    diffs = []
    for i in range(Xeq.shape[0]-1):
        diff = np.abs(np.sum(Xeq[i+1]-Xeq[i]))
        diffs.append(diff)
    diffs = np.array(diffs)

    plt.plot(diffs)
    plt.yscale("log")
    title = input("Plot title")
    plt.title(title)
    plt.xlabel("Time")
    plt.ylim(0, 1e-9)
    