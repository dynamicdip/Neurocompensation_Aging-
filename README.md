# Dynamic Neurocompensation associated with the Healthy Aging Project 
Contributions of short- and long-range white matter tracts in dynamic compensation with aging
Priyanka Chakraborty, Suman Saha, Gustavo Deco, Arpan Banerjee, Dipanjan Roy
https://academic.oup.com/cercor/article/35/2/bhae496/7945604


## Background
Optimal brain function is shaped by a combination of global information integration, facilitated by long-range connections, and local processing, which relies on short-range connections and underlying biological factors. With aging, anatomical connectivity undergoes significant deterioration, which affects the brain’s overall function. Despite the structural loss, previous research has shown that normative patterns of functions remain intact across the lifespan, defined as the compensatory mechanism of the aging brain. However, the crucial components in guiding the compensatory preservation of the dynamical complexity and the underlying mechanisms remain uncovered. Moreover, it remains largely unknown how the brain readjusts its biological parameters to maintain optimal brain dynamics with age. In this work, we provide a parsimonious mechanism using a whole-brain generative model to uncover the role of sub-communities, comprised of short-range and long-range connectivity, in driving the dynamic compensation process in the aging brain
We utilized two neuroimaging datasets to demonstrate how short- and long-range white matter tracts affect compensatory mechanisms. We unveil their modulation of intrinsic global scaling parameters, such as global coupling strength and conduction delay, via a personalized large-scale brain model. Our key finding suggests that short-range tracts predominantly amplify global coupling strength with age, potentially representing an epiphenomenon of the compensatory mechanism. This mechanism explains the significance of short-range connections in compensating for the major loss of long-range connections that occurs during aging. This insight could help identify alternative avenues to address aging-related diseases where long-range connections are significantly deteriorated.

Keywords: aging, white matter fiber tracts, subcommunity, metastability, compensation

Overview 

<img width="1562" height="775" alt="Screenshot 2026-03-04 at 8 26 39 AM" src="https://github.com/user-attachments/assets/4ed4d247-8265-4cc0-b5f7-546bb6166e99" />

#### young
The folder 'young' has the young subjects' data. The first eight characters are the id of the subject, and the last two numeric characters are the age of the subject. 

#### old
The folder 'old' contains data for old subjects. The first eight characters represent the subject's ID, and the last two digits represent the subject's age.


#### Empirical

The 'empirical_SC_analysis.m' file contains the empirical structural connectivity (SC) analysis.

####  Simulation
The sub-folder 'simulation/' contains all the simulation files and required function files. 
'BOLD.m' simulation of BOLD signal.

'Joana_Kuramoto_GlobalCoupling.m' is the Kuramoto model simulation.
Ref:  Cabral J, Hugues E, Sporns O, Deco G.
    Role of local network oscillations in resting-state functional connectivity.
   Neuroimage. 2011 Jul 1;57(1):130-9. Epub 2011 Apr 12.

'run_simulation.m'  Simulates individual subjects and determines the optimal values for k and tau.

'sim_meta_map.m'  Function file for simulation and saving metastability and synchronization maps for individual subjects.

#### Example
In the sub-folder 'examples/', we provide examples on how to run the code. Some simulation results are provided in the 'output' folder.
