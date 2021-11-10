import numpy as np
from tqdm import tqdm
from numba import jit

@jit(nopython=True)
def genotype_birthrate(total_pop, K, b0, sti):
    '''
    This function calculates a birthrate at a given timestep for a single allele. K is the carrying capacity.
    b0 is the maximum birthrate. If total population is greater than the carrying capacity, then the birhtrate is 0.
    Otherwise, the birhtrate is scaled so it is larger when the population is smaller and decreases as population
    approaches the carrying capacity. sti is the selection coefficient of the genotype at the time step and it scales
    the birthrate.
    '''
    T = total_pop/K
    if (T > 1):
        m = 0
    else:
        m = b0*(1-T)

    genotype_birthrate = m*(1+sti)
    return genotype_birthrate

@jit(nopython=True)
def mutation(mut_matrix, genotype):
    '''
    This function defines the possibility of mutation between genotypes, which is occurs at birth and is independent
    of selection. genotype is a number [0, 31] and mut_matrix is a 32x32 matrix describing probability of mutation from
    index i to index j. With a total mutation rate of .05, the diagonal entries will have a value of .95, and for
    row 0, describing mutation from allele 00000(0), entries 00001(1), 00010(2). 00100(4), 01000(8), and 10000(16) will each
    have a value of m = .01. The other entries will be 0.
    '''
    possibles = []

    #Scan the mutation matrix and make a list to hold the alleles the genotype can mutate to, and their probabilities.
    for i in range(len(mut_matrix[genotype, :])):
        if (mut_matrix[genotype, i] > 0):
            possibles.append((i, mut_matrix[genotype, i]))


    possibles = np.array(possibles)
    rnd = np.random.rand()

    #print(possibles)

    #If the random number is on [0, 1-5m), the offspring has the parental genotype. If it's on [1-5m, 1-4m), [1-4m, 1-3m), [1-3m, 1-2m), [1-2m, 1-m), [1-m, 1), it changes to the respective mutated genotype.
    i = 0
    while (rnd >= np.sum(possibles[:i, 1])):
        new_allele = possibles[i, 0]
        i += 1

    return int(new_allele)


@jit(nopython=True)
def sim_step(genotypes, current_genotype_pops, deathrate, total_pop, K, b0, mut_matrix, st):
    '''
    Defines one step of the simulation. genotypes is the number of genotypes in the simulation (32 in our model). Current_genotype pops is a vector of length "genotypes" where each entry is the population of that genotype. deathrate, K, and b0
    are parameters of the simulation. mut_matrix is a matrix of size "genotypes" x "genotypes". describing probability of mutation from index i to index j. With a total mutation rate of .05, the diagonal entries will have a value of .95, and for
    row 0, describing mutation from allele 00000(0), entries 00001(1), 00010(2). 00100(4), 01000(8), and 10000(16) will each have a value of m = .01. The other entries will be 0. st is a vector of length "genotypes" where each entry is the
    selection coefficient for that genotype at a single timestep.
    '''

    #Defube a vector of 0s to hold the populations of each genotype as they get counted up at the timestep.
    new_genotype_pops = np.zeros(genotypes)

    #K is carrying capacity. If population is over carrying capacity, birthrate m = 0. Otherwise the intrinsic birthrate b0 is scaled by how close the population is to K.
    T = total_pop/K
    if (T > 1):
        m = 0
    else:
        m = b0*(1-T)

    for i in range(genotypes):
        #Generate a random number for death probability for each individual of genotype i. Store values in a vector.
        rnd_vec_death_i = np.random.rand(int(current_genotype_pops[i]))

        #Deathrate is a probability, so if the random number > d, the individual gets counted at that time step.
        for j in rnd_vec_death_i:
            if (j >= deathrate):
                new_genotype_pops[i] += 1

        #The birthrate of each genotype is scaled by its selection coefficient at the time step.
        genotype_birthrate = m*(1+st[i])

        #Generate random numbers to realize probabilities of dividing and offspring mutating.
        rnd_vec_birth_i = np.random.rand(int(current_genotype_pops[i]))
        mut_probs = np.random.rand(int(current_genotype_pops[i]))

        #mut_matrix is scanned for entries between 0 and .99 to find the alleles one mutational distance away that the given genotype can mutate into.
        possible_muts = []
        for j in range(len(mut_matrix[i, :])):
            if (0 < mut_matrix[i, j] < .99):
                possible_muts.append(j)


        possible_muts = np.array(possible_muts)

        #birth_test is a vector with length equal to the number of individuals with genotype i. For an entry in birth_test, if it is >0, the individual gives birth. otherwise it does not.
        birth_test = genotype_birthrate - rnd_vec_birth_i

        for j in range(len(birth_test)):
            if (birth_test[j] >= 0) and(mut_probs[j] < mut_matrix[i, i]):
                #The offspring will most likely have the same genotype as the parent.
                    new_genotype_pops[i] += 1
            if (birth_test[j] >= 0) and (mut_probs[j] > mut_matrix[i, i]):
                #Mutation is independent of selection, so probabilities of mutating to any one of the genotypes one mutational distance away are equal.
                    new_genotype_pops[int(np.random.choice(possible_muts))] += 1
            else:
                #No birth from that individual
                new_genotype_pops[i] += 0

    #The if loop cycles through all genotypes before the new_genotype_pops vector is returned.
    return new_genotype_pops

class population (object):

    #For ease of use the simulation is written down in an object-oriented manner.
    def __init__(self, current_genotype_pops, fitnesses, mut_matrix, deathrate = .05, b0 = 2, K = 5e6, genotypes = 32):

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
        #Convenience function to run simulation for 45000 time steps and see what the genotype populations are at the end step.
        for i in tqdm(range(45000)):
            self.sim_step(si[i])
            genotype_pops.append(self.new_genotype_pops)
            self.current_genotype_pops = self.new_genotype_pops

        self.genotype_pops = genotype_pops

        return self.genotype_pops

    def run(self, nsteps, si):
        if (nsteps > si.shape[0]):
            return "ERROR: number of simulation steps exceeds length of selection coefficient matrix" #Selection coefficients need to be provided for every simulation step.

        genotype_pops = []

        for  step in tqdm(range(nsteps)):
            genotype_pops.append(self.current_genotype_pops)
            new_genotype_pops = self.sim_step(si[step])
            self.current_genotype_pops = new_genotype_pops
            self.total_pop = np.sum(self.current_genotype_pops)

        self.genotype_pops = np.array(genotype_pops)

        return self.genotype_pops #self.genotype_pops is an array with a column for every time point containing the frequencies of each genotype as the row entries so you can plot frequencies as function of time.


    def reset(self, i):
        #Unless you want to use the final populations from one simulation as the starting point for the next you will need to reset the current_genotype_pops variable before running the next simulation.
        current_pops0 = np.zeros(self.genotypes)
        current_pops0[i] = 10000
        self.current_genotype_pops = current_pops0
