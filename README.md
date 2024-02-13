# Neurocompensation_Aging Project 
# Contributions of sub-community based on short and long-range white matter tracts in personalized age-associated neurocompensatory 
# mechanism

## Background
This study investigates how the aging brain adapts its biological parameters, such as coupling between brain areas and time delays, in response to age-related deterioration of white matter tracts. Using fMRI data and computational modeling on large cohorts, the study reveals that short-range connections primarily modulate global coupling strength to compensate for structural loss, while long-range connections contribute to conduction delay. The aging brain calibrates these parameters to sustain functional integrity and achieve optimal metastability. 

#### young
The folder 'young' has young subjects data. First eigth charecters are the id of the subject and last two numeric are the age of the subject. 

#### old
The folder 'old' contains data for old subjects. The first eight characters represent the subject's ID, and the last two digits represent the subject's age.


#### Empirical

The 'empirical_SC_analysis.m' file contains the empirical structural connectivity (SC) analysis.

####  Simulation
The sub-folder 'simulation/' contains all the simulation files and required functions files. 
'BOLD.m' simulation of BOLD signal.

'Joana_Kuramoto_GlobalCoupling.m' is the Kuramoto model simulation.
Ref:  Cabral J, Hugues E, Sporns O, Deco G.
    Role of local network oscillations in resting-state functional connectivity.
   Neuroimage. 2011 Jul 1;57(1):130-9. Epub 2011 Apr 12.

'run_simulation.m'  Runs the simulation for individual subjects and determines the optimal values for k and tau.

'sim_meta_map.m'  Function file for simulation and saving metastability and synchronization maps for individual subjects.

#### Example
In the sub-folder 'examples/', we provide examples on how to run the code. Some simulation results are provided in the 'output' folder.
