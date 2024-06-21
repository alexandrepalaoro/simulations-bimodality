# Sigmoid allometries generate male dimorphism in secondary sexual traits: a comment on Packard (2023).

Authors: Bruno Buzatto, Glauco Machado, Alexandre V. Palaoro <br>
Journal: Evolutionary Ecology, DOI: 10.1007/s10682-024-10303-6 (<a href="https://link.springer.com/article/10.1007/s10682-024-10303-6" target="_blank">link</a>) <br>
Contact about code, and analyses: alexandre.palaoro@gmail.com

---

### This readme has been divided in three parts. First, we will talk about file structure, then the code, and then the dataset.

##### File structure:

We uploaded four different folders: "code", "data", "figures", and "simul_data". Each folder contain the files we used to run all analyses and the output files that came from the codes (the figures and simulation data). 

##### Code:

The folder "code" contains all code used in this paper. The code is divided in files used to run the analyses of the real data, and files used for the simulations. Codes that start with "Function" are required to run the analyses used in the "Simulating sigmoid allometries and detecting dimorphisms.R". These are just some functions to make it easier to analyze everything. The two code files that start with "simulation_" contain the raw code used for the simulations themselves. If you have the folder "simul_data", you do not need to use these codes and can use the resulting simulations directly with the main code: "Simulating sigmoid allometries and detecting dimorphisms.R"

##### Data:

Our data is divided in the folder "data" and the folder "simul_data".
In the "data" folder, we have the real harvestmen data used in the original papers. The name of the file relates to the species, and inside the folder you will have two columns, one related to the size of the body and another to the size of the weapon. 

In the "simul_data" you will find .RDS files that contains the simulations. You can run your own simulations, but these are the files we used to generate the tables and figures for the paper. 

The code was run with RStudio (v2023.06.2) in R software (v4.3.2). 

#### Figures:

The folder "figures" contain the figures that were taken directly from R. In some of them, we made some minor changes in Powerpoint.

##### Acknowledgments:
We thank Gary Packard for triggering this interesting debate on the topic of male dimorphism and static allometry, and also for revising an early version of the manuscript. We also thank an anonymous reviewer for comments and Matthew Symonds and Emma Sherratt for considering this submission for Evolutionary Ecology.
