# Tracking bacterial inter- and intraspecific competition in pairwise assays using a transient fluorescent dye

This repository contains data and code for the manuscript titled: "Tracking bacterial inter- and intraspecific competition in pairwise assays using a transient fluorescent dye".

This manuscript presents a new methodology to track relative abundances of two competing bacterial types based on staining the population of one competitor with a transient fluorescent dye. Competitor identity during the competition is then determined based on single cell fluorescence intensity. Classification of stained and unstained cells is based on fluorescence thresholds determined by Receiver-Opterator Curves (ROC). Julius Hoffmann (JH) and Lutz Becks (LB) designed this methodology and the experiments and Ole Johannsen (OJ) provided tools and analyses to determine fluorescence thresholds for classification.  

For the manuscript presented three experiments were performed: a staining test, a competition assay with two genotype of Pseudomonas fluorescens and a competition assay between Pseudomonas fluorescens and Escherichia coli. Data for each experiment is stored in a separate folder of this repository. For the staining test there is data for the classification of stained and unstained cells and the measured fluorescence values. For both competition experiment there is data for the classification of stained and unstained cells, data for the evaluation of competitions on growth agar and optical density measurements from the competition cultures. 

The folder R scripts contains programming code that processes the raw data in the three experiment folders.
The main folder contains three R markdown script that each contain the analyses for one of the three experiement. The scripts can be compiled as html reports. Output from theses script is saved to the "Results" folder. 
