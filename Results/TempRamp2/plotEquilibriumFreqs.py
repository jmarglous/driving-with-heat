# -*- coding: utf-8 -*-
"""
Created on Fri Jul 16 10:42:04 2021

@author: jmarg
"""
import numpy as np
import matplotlib.pyplot as plt

Xeq = np.loadtxt("Xeq.txt", delimiter = ",")

fig = plt.figure()
ax = fig.add_subplot(111)

for i in range(32):
    ax.plot(Xeq[:, i])
    
ax.set_xlabel("Time")
ax.set_ylabel("Allele frequency")
ax.set_yscale("log")
title = input("Plot title: ")
ax.set_title(title)

fig.savefig("Equilibriumfreqs.png", delimiter = ",")