Hello! If you are reading this you are likely a student in the Weinreich lab who has decided to try their hand at this incomplete project. Good luck! I hope that you're able to succeed where I have failed -- or, at the very least, you learn that the true gift is the friends and knowledge we make along the way ;).

In this ReadMe I'll provide a little background on the science and motivation behind the project and a summary of my contributions and code to get you started. I would love to hear from you if you'd like to talk about the project as well -- for the foreseeable future you'll be able to reach me at my Brown email, jacob_marglous@brown.edu. One final note before we begin -- I studied physics in undergrad only because I was having trouble memorizing things in biology classes, and I only learned my minimal coding skills *ad hoc* for my physics classes and research projects. This is to say that if you find something fishy in the project documents, please ask me before you waste too much time trying to figure it out -- experience teaches that in all likelihood I'm the silly one. 

Without further ado...

## Background 
This project is centered around a subfield in mathematical oncology and population genetics with the buzzword *evolutionary driving*. The fundamental idea is: generally, genotypes correspond with phenotypes. And, we know that over time a biological population (we're thinking about pathogenic populations, so picture a group of cancer cells or a bacterial infection) will respond to environmental pressures by selecting for members of the population with fitter phenotypes. So, if a tumor is treated with chemotherapy or a bacterial infection with antibiotics, the *genotype frequencies* of the genotypes corresponding with fitter phenotypes will increase, and the genotype frequencies of the less fit genotypes wil decrease. When treating disease, this is **bad** -- it means drug resistance, metastasis (cancer), sepsis (infection), etc. 

To a large extent, we can't control this: variation always leads to selection, evolution, etc. What can we control? The environment itself -- the types and sequence of drugs we treat the population with. Here's an example:

	Imagine a population of pathogens with mostly genotype 2, that is susceptible to drug A but resistant to drug B. Then, treating with antibiotic A causes genotype 1 to become nearly fixed (nearly the only genotype) in the population. This happens because when the population is exposed to A, individuals with genotype 1 are able to continue to reproduce the most efficiently. But, genotype 1 is also susceptible to antibiotic B. So if you treat with drug A for a sufficient time for the population to become nearly all genotype 1, then switch to drug B, the population will be nearly or completely wiped out.
	
This is evolutionary driving in a nutshell: you've used evolution to outsmart your pathogen. The population thinks it's winning as it climbs the fitness landscape to a fit genotype, but then you reveal yourself to be the true mastermind as you open a trap door as soon as it reaches the peak! *Checkmate*. You can also imagine how complicated this gets. Here are some good research questions that people have started to answer: which pairs of drugs (for, say, a particular type of lung cancer) induce collateral sensitivity -- i.e., resistance to drug 1 increases sensitivity to drug 2 (Dhawan et al)? Are such pairs even commmon (Nichol et al.)? Could we move beyond pairs and think about chains of 3 or more collaterally sensitive treatments? 

Ultimately, the goal is to use rapid sequencing in the clinic and more complete knowledge of drug sensitivity interactions to tailor treatments to the pathogenic populations in patients in real time. Rather than eventually succumbing to multi-drug-resistant metatstatic cancer, patients could be helped by evolutionarily-informed treatments that keep their tumors susceptible to some drug in the toolbox. Imagine: instead of trying to *cure* cancer, we could domesticate it into a manageable hronic disease.

Setting aside the direct clinical challenges, you might have a couple of concerns. Here are some big ones:
1. You've called this idea *evolutionary driving*, but you're hardly Max Verstappen piloting his Red Bull F1 car here. Maybe you can control the car a little bit by dictating the sequence of the tracks, but after all, there are only so many tracks available to you, and once he's on a given one, he's going to go where he wants to. 
2. How do you even know when the race is over -- that is, how do you know when to switch from drug A to drug B (or when to switch from the Indianapolis Motor Speedway to the Nurburgring (there's a really good documentary about F1 on Netflix called Drive to Survive -- is it obvious that I just finished it?))

This is where this project comes in.

## Motivating counterdiabatic driving
So, two major shortcomings of evolutionary driving in its basic form is that the number of evolutionary trajectories offered to the pathogen is limited, and the time it takes to climb the path is unknown. Getting the timing right is really important in a clinical setting because you want to treat the patient as quickly as possible, but if you make the collateral switch too soon, you'll provide opportunities for resistance to develop rather than susceptibility. A method called counterdiabatic driving  offers a theoretical solution of fine control over the entire treatment, not just the binary switch between drugs. 

Here's the idea: at the start of a treatment plan, the genotypes in the naive population are each going to have some frequency (genotype 1 is .95, genotype 2 is .005, genotype 3 is .01, genotype 4 is .015, etc.). A treatment plan dictates the concentration of drug the population is exposed to at each time on a schedule. At each drug concentration, each of the genotypes has a "fitness", which is a defined metric that's essentially a scaled growth rate between 0 and 1. *Each genotype has a fitness which is a function of drug concentration, and the drug concentration is a function of time.* So over time, the genotype frequencies will change due to the changing fitnesses, which change due to the the changing drug concentration. Imagine a treatment plan where drug concentration increased monotonically from 0, and then stayed at a higher concentration. The genotype frequencie would fluctuate, and then settle on a new equilibrium favoring whichever genotype is fittest at the high concentration (that genotype would not become the *only* genotype at equilibrium, because random mutation always allows for some flux between the genotypes -- but more on the model later). 

You could stop the increasing stage of the drug schedule at any time, and the results would be generally the same -- at any drug concentration, there will be a set of equilibrium genotype frequencies. You can even calculate what this equilibrium frequency will be! (See paper and thesis ) But it takes time for genotype frequencies to change, and that time is related to the lifespan of the organisms -- it takes time for generations to grow up, and old ones to die out. (Just ask Pierre Gasly, who's been waiting for a seat in one of the two Red Bull cars to open up for years now.) So, at any moment *while* the drug concentration is changing, the true population frequencies will lag behind the equiibrium frequency at that time (remember: the equilibrium frequencies are determined by the drug concentration, but the drug concentration is determined by the time.)

So, a question: 
1. Given an original drug treatment plan, can we modify it so the population is at equilibrium frequencies (or as close to them as possible) at every moment on the schedule?

And three more related ones:
		a. How would you do this? 
		b. How could you test it?
		c. And why is it useful?

So, to the first question: if you think back to thermodynamics classes you may rememer a lot of problems about compressing a gas with a piston. To compress the piston from point a to point b instantaneously, you can do work $W = P\Delta V$ to overcome the pressure of the gas pushing against the piston, where $P$ is the final pressure of the gas. But performing this compression instantaneously requires the maximum possible amount of work, because as the gas becomes more and more compressed, its pressure increases. The gas will also be knocked out of equilibrium during the compression because the external pressure you apply will be greater than the pressure of the gas. If you do the compression in $n$$ small steps of $W_n = P_n\delta V$, you'll reduce the amount of work needed because you won't be applying excess pressure at the beginning of the compression. The gas will also re-equilibrate between each mini compression of $\delta V$. The way to minimize the work required (you may guess where this is going) is to perform an infinite number of small steps  of size $dV$. Because there's an infinite number of steps, they also require infinite time to perform, and the gas will always be at equilibrium. This is called *adiabatic* compression or expansion. 

The same idea applies to particles/waves in a box in quantum mechanics. The allowed equilibrium wavelengths of a particle in a box is determined by the Schrodinger equation, which depends on the length of the box. If you change the width of the box, you'll knock your system out of equilibrium. But if you change the width of the box *infinitely slowly*, the quantum system will remain in equilibrium throughout the whole process. 

If we were both (1) interested in changing systems while maintaining them at equilibrium, and (2) lived for an infinite time, this would be useful. Unfortunately, only the first thing is true, and so there's been a lot of effort devoted to finding *shortcuts to adiabaticity*, also called *counterdiabatic driving* methods, that change the energy of systems as quickly as possible while keeping them at equilibrium. 

	One demonstration explained to me is a waiter carrying a tray of full glasses of water. To avoid spilling, she has to walk super slowly from the kitchen to the table. But if she tips the tray forward a little bit to balance her momentum, she can walk more quickly. Tipping the tray forward is a **shortcut to adiabaticity**.
		
How does this relate to our system? It turns out that our evolving system, which is described by a *Wright-Fisher model*, is mathematically equivalent to the waveform governed by the Schrodinger equation. Maybe you can see a big clinical problem -- if we attempted to implement collateral sensitivity and wanted to wait for the population to equilibrate to drug A before dosing with drug B, we'd have to wait forever. But, we can also cross-apply the shortcuts to adiabaticity the physicists have drawn up to our system, so that we can implement our ideal treatment plans before the patient dies. 

Here's a sketch of how this mechanism would work. You'll
1. pick a drug treatment plan
2. Calculate the equilibrium genotype frequencies at every timepoint of the plan. 
	- You'll do this by writing down a function for the flux between each of the genotypes and finding the genotype frequencies at each time point for which that function's derivative is 0. 
3. Use an agent-based model (ABM) to see how the population will actually evolve under that treatment plan. 
	- You should observe the actual realized genotype frquencies **lag** behind the equilibrium frequencies throughout the treatment plan. 
4. Apply counterdiabatic driving to calculate a new treatment plan that matches the population to the original equilibrium frequencies at every time point. 
5. Apply that treatment plan to the ABM and observe that the realized genotype frequencies are now matched to the equilibrium frequencies, as closely as possible. 

This is the outline of the Iram et al paper which inspired this project. Please see that paper, and to a lesser extent, my thesis (both linked below) for more of the mathematical details!

## My contribution
Iram et al took fitness data for 16 bacterial genotypes in response to an antibiotic and did some pretty math to show they can control the rate of change of the genotype distribution as the population evolves in response to that drug. It's really exciting, but also a bit unsatisfying. In the end, the fittest genotype at high drug concentrations will always emerge, so even if we can control the timing of the evolution, we still don't have great control over the endpoint. The Weinreich lab had data on antibiotic fitness for 32 *E. coli* mutants at 6 different temperatures, and we hoped that with the addition of temperature to drug concentration as a second control parameter, we could gain additional control over the evolving population. 

The 32 alleles also represented a combinatorially complete set of alleles. That is, there are 5 loci with a wildtype (0) or mutant (1) allele possible at each locus. So, 00000 is the wildtype, and 11111 is the full mutant. Each allele is 1 mutational distance away from 5 others. 00000 is 1 mutational distance away from 00001, 00010, 00100, 01000, and 10000, and so on for each allele. 

This is described in more detail in my thesis, but the data that was available for the 32 mutants was IC95s -- the drug concentration for each allele, at each of 6 temperatures, where 5% of individuals of that allele were able to survive. As you may imagine, all of the math sketched here requires fitness values for each allele at *every* possible drug concentration. So, I converted IC95, a single fitness value, to a growthrate curve, in order to build growthrate functions dependent on temperature and drug concentration for each of the 32 alleles. 

Then, I wrote an ABM to observe how the population might respond to various drug and temperature treatments. (Iram et al. had built an ABM for their paper, but it was implemented in C and written for their 16 genotypes rather than my 32 -- so, I wrote my own.) The idea is as follows: at every timepoint, each individual

- Can divide (produce a single offspring, plus itself), with probability $b$, determined by the maximal birthrate, population size, and selection coefficient of the individual's genotype at that time step. 
	- If it divides, the offspring can be mutated to any of the mutationally adjacent alleles, each with equal probability $m$. 
- Can die with probability $d$ (independent of giving birth)
Then, you step forward in time and repeat this process for a new set of selection coefficients. 

## What went wrong
Unfortunately, I had trouble calculating equilibrium concentrations for the 32 alleles. Generally, the change in equilibrium frequencies between successive time steps is small, so Iram et al. wrote a function that solves for equilibrium frequencies computationally at the first step, then uses that solution as the starting guess for the equilibrium solution at the next step. For reasons still unknown, for 32 genotypes the solver frequently finds impossible equilibrium frequencies (i.e. <0 or >1). I've tried a variety of ideas (you can see them in the attached scripts) to guide the solver to more reasonable solutions, but I haven't hacked it yet. This is where you come in!

## What needs to be done
The overarching goal of this project is to replicate the Iram et al paper with two control parameters instead of one. 
- To that end, the first major goal moving forward is to find reliable equilibium frequencies for any hypothetical drug-and-temperature treatment plan. 
- With that complete, you'll be able to use the ABM to compare actual behavior of the population to those equilibrium frequencies. 
- Then, you can calculate the CD modification to the original protocol, update the ABM, and observe that the lag from equilibrium is reduced!

### Possible problems
Because I was stuck on equilibrium frequencies, I wasn't able to think about many of the CD driving aspects of the project. Calculating the CD modification to the original protocol is done in the space of selection coefficients. That part shouldn't be too bad, but because the selection coefficients depend on two control parameters rather than one as in Iram et al, I'm not sure how to relay the desired selection coefficients at each time for the CD protocol into a new drug-and-temperature protocol. This is a crucial part of the project and one I hope you'll be able to figure out. 


## Guide to the repository
- Evolution-Counterdiabatic-Driving-1.0 is my local copy of the Github repository containing the work for the Iram et al. paper (linked below). The repository is worth looking through!
	- One very useful thing in this repository: to generate mutation matrices, navigate to `Evolution-Counterdiabatic-Driving-1.0/landscapes/scripts` in the terminal and type
	```
	python make_trans_prob.py [number of genotypes] [total mutation rate of out of genotype] > ./[desiredfilename].csv
	
	```
	to generate a mutation matrix in the correct form for the scripts in this project. 
	
- Empirical is a library that Evolution-Counterdiabatic-Driving-1.0 depends on.
- ScottLabDataProcessing folder contains some plots of the raw IC95 data for the 32 genotypes generated by myself and Nikhil Krishnan, a former student in the Scott lab. It also has some Jupyter notebooks he used (and I adapted) to make those plots from the data.
- ShamreenNotebooks contains the code to calculate and plot equilibrium frequencies, plot results from the ABM, and calculate CD modifications from the original paper, in Jupyter notebook form. 
- weinreich_temp_modifiedlabels.csv is the raw data the project is based on. Each column is a temperature (20, 25, 30, 35, 37, 41C), and each row is MICs for the given allele at those temperatures.
- mutdat32_fromscript.csv is a mutation matrix for 32 alleles with a total mutation rate from each genotype as .001.
- tempConcAnalysis.py contains functions to interpolate the original IC95 data into usable fitness landscapes. 
	- - Heatmaps.ipynb has some plotting functions to visualize the landscapes you make. There's some overlap with tempConcAnalysis.py but I was fighting with seaborn (a plotting library) to make the graphs and I'm scared to delete anything now!
- Running seascapeGeneration.py in the terminal will guide you through setting a treatment protocol and calculating selection coefficients at each time step for each genotype. This script should make tempConcAnalysis.py redundant, but again -- scared to delete it!
- TestingNumericalSolns.py contains functions for calculating equilibrium frequencies based on the selection coefficients that seascapeGeneration makes. These functions are based on the original functions for 16 genotypes in Shamreen's notebooks. They are also the part of the project i've spent a lot of time on because for some treatment plans they can provide unreliable solutions. 
	- DeterministicPrediction.py is an attempt to guide the solver in the right direction when it is being recalcitrant. The thought is that the numerical solver used to calculate equlibrium typically uses solutions from the last time step to guess the solution for the current one. But in the 32-dimensional space it is solving in, the solver is very sensitive to small perturbations away from the desired solution. The frequency updating equation which is set to 0 at equilibrium has a stochastic part and a deterministic part dependent on the genotype selection coefficients (see Iram et al.). My reasoning here was that if I updated the equilibrium genotypes from the last step with the deterministic portion of the updating formula, and repeated this many times, I woud approach the equilibrium frequencies. I could plug this updated "deterministic prediction" into the solver and hopefully with this added push it would find the desired solution. It unfortunately only works for some cases. 
- plotEquilibriumFreqs.py is a script of convenience to plot the output from the solving functions as genotype frequencies vs. time. 
- The ABM folder contains:
	- Simulation_draft2.py, the ABM as written in Python.
	- RunSims.py, a script that made it easier for me to run the simulation from the command line. Be sure to check the simulation parameters in this one before you use it too heavily to make sure they match your needs. 
	
-The Results folder contains:
	- Tables of selection coefficients and calculated equilibrium frequencies, and plots of treatment protocols and equilibrium frequencies for a number of sample protocols. These were generated with the code in ...

	

## Resources
### Papers to get you started (look in their references for more)
- [Iram et al., "Controlling the speed and trajectory of evolution"](http://www.nature.com/articles/s41567-020-0989-3)
	- This is a foundational paper behind this project. There are also good resources cited in this one if you're looking for more information about shortcuts to adiabaticity or evolutionary driving. The SI is long and challenging but useful.
- [Dolson et al, "An introduction to using counterdiabatic driving"](https://www.mitpressjournals.org/doi/abs/10.1162/isal_a_00344)
	- This is a conference abstract by the same authors that may be helpful. 
- [Dhawan et al, "Collateral sensitivity networks"](10.1038/s41598-017-00791-8)
	- A good look at collateral sensitivity. 
- [Nichol et al, "Antibiotic collateral sensitivity is contingent on the repeatability of evolution"](http://www.nature.com/articles/s41467-018-08098-6)
	- A good critical look at the limits of collateral sensitivity, and a motivator for the importance of counterdiabatic driving, 

### Other Resources
- A [Github repository](https://github.com/hincz-lab/Evolution-Counterdiabatic-Driving) with the code from the Iram et al. paper.
- A Twitter thread from Dr. Jacob Scott, one of the coauthors of the Iram et al paper, can be found [here](https://twitter.com/CancerConnector/status/1298270321422692354?s=20). It has lots of good background and links to papers motivating this project. 
-  A similar, more physics-focused threat from Professor Mike Hinczewski, another coauthor and Shamreen Iram's graduate advisor, can be found [here](https://twitter.com/hinczewskilab/status/1301338147251523584).
-  A nice brief conference video from Professor Emily Dolson, another coauthor, is [here] (https://www.youtube.com/watch?v=uJOgtbW57Mw).
-  Here is a recorded meeting from last year where a lot of the contributors to the project discuss it. I found it to be [really helpful background](https://www.youtube.com/watch?v=qonnyPlsyQw).
-  My [thesis](https://www.brown.edu/academics/physics/sites/physics/files/Copy%20of%20Marglous%2C%20Jacob%20POST.pdf) from last year. I made an efford to explain some of the math at a simpler level than some of the resources, so hopefully it's helpful!