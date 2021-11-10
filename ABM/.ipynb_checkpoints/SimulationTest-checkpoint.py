import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
from numba import jit, int32, float32
import time

sCyc = np.loadtxt("./CDDriving_PythonVersion/sCyc_orig.csv", delimiter = ",")
mutdat = np.loadtxt("./CDDriving_PythonVersion/InputData/mut_matrix_.001.dat", delimiter = ",")

sFlat = np.loadtxt("./Const_st_T27C10e-02.csv", delimiter = ",")
mutdat32 = np.loadtxt("./mutdat32_fromscript.csv", delimiter = ",")

current_pops0 = (.04/16)*100000*np.ones(16)
current_pops0[14] = .96*100000

@jit(nopython=True)
def genotype_birthrate(total_pop, K, b0, sti):
    T = total_pop/K 
    if (T > 1):
        m = 0
    else:
        m = b0*(1-T)
    
    genotype_birthrate = m*(1+sti)
    return genotype_birthrate

@jit(nopython=True)
def mutation(mut_matrix, genotype):
    possibles = []
    
    for i in range(len(mut_matrix[genotype, :])):
        if (mut_matrix[genotype, i] > 0):
            possibles.append((i, mut_matrix[genotype, i]))
    
    possibles = np.array(possibles)
    rnd = np.random.rand()
    
    #print(possibles)
    
    i = 0
    while (rnd >= np.sum(possibles[:i, 1])):
        new_allele = possibles[i, 0]
        i += 1
        
    return int(new_allele) 

         
@jit(nopython=True)
def sim_step(genotypes, current_genotype_pops, deathrate, total_pop, K, b0, mut_matrix, st):
    
    new_genotype_pops = np.zeros(genotypes)
    
    rnd_birth = np.random.rand(genotypes, int(np.max(current_genotype_pops)))
    rnd_death = np.random.rand(genotypes, int(np.max(current_genotype_pops)))    
    
    for i in range(genotypes):
        for j in range(int(current_genotype_pops[i])):
            
            if rnd_death[i, j] - deathrate >= 0:
                new_genotype_pops[i] += 1
            
            if (genotype_birthrate(total_pop, K, b0, st[i]) - rnd_birth[i, j] >= 0):
                child_genotype = mutation(mut_matrix, i)
                new_genotype_pops[child_genotype] += 1
        
    return new_genotype_pops


class population (object):
    
    def __init__(self, current_genotype_pops, fitnesses, mut_matrix, deathrate = .05, b0 = 2, K = 5e6, genotypes = 16):
        
        self.current_genotype_pops = current_genotype_pops
        self.fitnesses = fitnesses
        self.mut_matrix = mut_matrix
        self.deathrate = deathrate
        self.b0 = b0
        self.K = K
        self.genotypes = genotypes
        
        self.total_pop = np.sum(self.current_genotype_pops)
         
    '''
    Below here I started to define a function that would allow input-control over the drug application. 
    But I don't want to worry about making the simulation extensible between the 16-genotype without temp 
    and the 32-genotype with temp yet so I'm going to just import the array of selection coefficient vectors
    directly and come back to this later once I've validated the simulation on the 16-genotype model. 
    '''
    
    #def makeSvectors(self, start_temp = 20, end_temp = 37, temp_ramp = .02, temp_time_toramp = 10110, start_conc = 1e-2, end_conc = 1e-2, conc_ramp = .002, conc_time_toramp = 10110, times = np.arange(0, 45001)):
        #tt = (end_temp-start_temp)/(1.+np.exp(-temp_ramp*(times - temp_time_toramp)))+start_temp
        #cc = (end_conc-start_conc)/(1.+np.exp(-conc_ramp*(times - conc_time_toramp)))+start_conc
    
    ## Need to figure out how to import sfc and interpolation function
    
    
    def genotype_birthrate(self, sti):
        
        
        birthrate = genotype_birthrate(self.total_pop, self.K, self.b0, sti)
        
        return birthrate
    
    def mutation(self, genotype):
        
        genotype_out = mutation(self.mut_matrix, genotype)
        
        return genotype_out
    
    def sim_step(self, st):
                
        new_genotype_pops = sim_step(self.genotypes, self.current_genotype_pops, self.deathrate, self.total_pop, self.K, self.b0, self.mut_matrix, st)
        
        self.new_genotype_pops = new_genotype_pops
        return self.new_genotype_pops
    
    def mult_steps(self, si):
        
        genotype_pops = []
        for i in tqdm(range(45000)):
            self.sim_step(si[i])
            genotype_pops.append(self.new_genotype_pops)
            self.current_genotype_pops = self.new_genotype_pops
        
        self.genotype_pops = genotype_pops
        
        return self.genotype_pops
    
    def run(self, nsteps, si):
        if (nsteps > si.shape[0]):
            return "ERROR: number of simulation steps exceeds length of selection coefficient matrix"
              
        genotype_pops = []
        
        for  step in tqdm(range(nsteps)):
    
            genotype_pops.append(self.current_genotype_pops)
            
            new_genotype_pops = self.sim_step(si[step])
            
            self.current_genotype_pops = new_genotype_pops
            
            self.total_pop = np.sum(self.current_genotype_pops)
            
        self.genotype_pops = np.array(genotype_pops)
        
        return self.genotype_pops
            
           
    def reset(self):
        current_pops0 = np.zeros(self.genotypes)
        current_pops0[0] = 10000
        self.current_genotype_pops = current_pops0
    
    
pops0_32 = np.zeros(32)
pops0_32[0] = 10000
            
test_pop = population(current_genotype_pops = current_pops0, fitnesses = sCyc[:, 0], mut_matrix = mutdat, deathrate = .05, b0 = 2)
test_pop32 = population(current_genotype_pops = pops0_32, fitnesses = sFlat[:, 0], mut_matrix = mutdat32, deathrate = .05, b0 = 2, genotypes = 32)
