import numpy as np
from tqdm import tqdm
from numba import jit

#sCyc = np.loadtxt("./CDDriving_PythonVersion/sCyc_orig.csv", delimiter = ",")
mutdat = np.loadtxt("./CDDriving_PythonVersion/InputData/mut_matrix_.001.dat", delimiter = ",")
#sCyc_drugRamp = np.loadtxt("./CDDriving_PythonVersion/sCyc_drugRamp_16Allele.csv", delimiter = ",")
s32_flat = np.loadtxt("./Const_st_T27C10e-02.csv", delimiter = ",")
s32_Tramp = np.loadtxt("./st_Tramp_Cconst.csv", delimiter = ",")
s8_Tramp = np.loadtxt("./st_8_from32_Tramp.txt", delimiter = ",")
mutdat32 = np.loadtxt("./mutdat32_fromscript.csv", delimiter = ",")
#s32_cRamp = np.loadtxt("st_Tconst27_Cramp.csv", delimiter = ",")
current_pops0 = (.04/16)*100000*np.ones(16)
current_pops0[14] = .96*100000

mutdat4 = np.loadtxt("mutdat4.txt", delimiter = ",")
Xeq4 = np.loadtxt("4Gen_Xeq.txt", delimiter = ",")
st_4_flat = np.loadtxt("st_4_from32_flat_8820.txt", delimiter = ",")

mutdat8 = np.loadtxt("mutdat8.txt", delimiter = ",")
Xeq8 = np.loadtxt("Xeq_8_orig.txt", delimiter = ",")
st_8_flat = np.loadtxt("st_8_from32_11424.txt", delimiter = ",")

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

    T = total_pop/K 
    if (T > 1):
        m = 0
    else:
        m = b0*(1-T)
    
    for i in range(genotypes):
        rnd_vec_death_i = np.random.rand(int(current_genotype_pops[i]))
        
        for j in rnd_vec_death_i:
            if (j >= deathrate):
                new_genotype_pops[i] += 1

        genotype_birthrate = m*(1+st[i])
        rnd_vec_birth_i = np.random.rand(int(current_genotype_pops[i]))
        
        mut_probs = np.random.rand(int(current_genotype_pops[i]))

        possible_muts = []
        for j in range(len(mut_matrix[i, :])):
            if (0 < mut_matrix[i, j] < .99):
                possible_muts.append(j)

        possible_muts = np.array(possible_muts)

        birth_test = genotype_birthrate - rnd_vec_birth_i
    
        for j in range(len(birth_test)):
            if (birth_test[j] >= 0) and(mut_probs[j] < mut_matrix[i, i]):
                    new_genotype_pops[i] += 1
            if (birth_test[j] >= 0) and (mut_probs[j] > mut_matrix[i, i]):
                    new_genotype_pops[int(np.random.choice(possible_muts))] += 1
            else:
                new_genotype_pops[i] += 0
            
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
            
           
    def reset(self, i):
        current_pops0 = np.zeros(self.genotypes)
        current_pops0[i] = 10000
        self.current_genotype_pops = current_pops0
    
    
            
#test_pop = population(current_genotype_pops = current_pops0, fitnesses = sCyc[:, 0], mut_matrix = mutdat, deathrate = .05, b0 = 2)

current_pops0_32 = np.zeros(32)
current_pops0_32[0] = 10000

test_pop32 = population(current_genotype_pops = current_pops0_32, fitnesses = s32_flat, mut_matrix = mutdat32, deathrate = .05, b0 = 2, K = 5e5, genotypes = 32)

pops_4 = np.round(10000*Xeq4[822, :])

pop4 = population(current_genotype_pops = pops_4, fitnesses = st_4_flat, mut_matrix = mutdat4, deathrate = .05, b0 = 2, K = 5e5, genotypes = 4)

Xeq0_8 = np.loadtxt("Xeq0_8.txt", delimiter = ",")
pops_8 = np.round(10000*Xeq0_8)

pop8 = population(current_genotype_pops = pops_8, fitnesses = st_8_flat, mut_matrix = mutdat8, deathrate = .05, b0 = 2, K = 5e5, genotypes = 8)