ReefMod-PALAU_v1.1


Liam Lachs, Newcastle University (liam.lachs@newcastle.ac.uk)

ReefMod-GBR_v6.8 by Yves-Marie Bozec, The University of Queensland (y.bozec@uq.edu.au)

READ ME

========================
General information
========================

Author: Liam Lachs; Yves-Marie Bozec
Contact: liam.lachs@newcastle.ac.uk; y.bozec@uq.edu.au
DOI: TBC
License: CC-BY 4.0
Last updated: 24/09/2024
Related article: Natural selection could determine whether Acropora corals persist under 
                 expected climate change (Science)

========================
Introductory information
========================
This repository contains the scripts of ReefMod-PALAU_v1.1, based on ReefMod-GBR_v6.8. It 
is a coral individual-based model that simulates 95 coral populations across the reef 
system of Palau. This version was used to explore the adaptive potential of corymbose 
Acropora corals under ocean warming, with a revised (compared to ReefMod-GBR_v6.8) 
implementation of coral bleaching. This model tracks coral heat tolerance at the individual 
level. The model is used to explore future adaptive trajectories across three different 
future emissions scenarios (Paris Agreement limiting warming to 2C: SSP1-2.6; 
Middle-of-the-road: SSP2-4.5; Worst-case: SSP5-8.5) and for selective pressure forcings 
(i.e., future marine heatwave stress) from 16 different global climate models.

Files included in the data deposit are the inputs to ReefMod-PALAU v1.1 and include: 

1) Future heat stress forcings
2) Coral parametrisations
3) Reef connectivity matrix

========================== 
Methodological information
==========================
1) System requirements
- All modelling was conducted and tested on MATLAB R2024a
- No non-standard hardware is required to run the main analysis.

2) Installation Instructions
- Install MATLAB (with a license) following https://uk.mathworks.com/help/install/ug/install-products-with-internet-connection.html.
- This should take less than 1 hour, if R is not yet installed.

3) Demonstration and Instructions for use
- Code functionality and descriptions are given as comments in the code.
- This will reproduce the output files underpinning presented in the paper.
