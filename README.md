# Proto-Metabolic Reaction Network Simulator (ProMetRNs)

Project title: Simulation of pre-biotic formation of metabolic networks with the appearance of heritable catalysts of random specificity

Author: AC Du Plessis, Prospective MSc Bioinformatics and Computational Biology

## Project Background
Although many _in silico_ simulations to model the Origin of Life (OoL) has been attempted, there is no simulation that encapsulates all aspect. As this is only a master's project, we will did not set out to solve the OoL, but rather to gain some insight into the catalyst vs replicator debate. With this simulation we are working from the assumption that life could have started in an alkaline hydrothermal vent system like the Lost City hydrothermal vent system. Using the minimal metabolism of JCVI syn3A as a guide, we simulate the prebiotic chemical environment of a hydrothermal vent system by adding "catalysts" to _in silico_ pores in a vent and verify whether a chemical reaction network can form based on the chemistry of this JCVI syn3A.  

## Running the scripts

### Preprocessing

#### Running instructions
The preprocessing script is not included in the main program to keep the script that has to be run on the HPC as small as possible. This also allows

#### Input file
As input we used the full_jcvi_rxn_list.csv file which was compiled from the supplementary materials from the Breuer _et al_. (2019) Essential metabolism of a minimal cell. 
Format:
> Reaction nr,Reaction equation
>
> S001,D-glucose 6-phosphate <==> D-fructose 6-phosphate
> 
> S002,D-fructose 6-phosphate + ATP --> ADP + D-fructose 1,6-biphosphate + H+

#### Expected output
Two files will be created: 
1. chem_dict.txt
2. rxn numbers.csv

The chemical dictionary contains all of the metabolites extracted from the metabolic reaction list provided and its assigned value (I had to assign values because when working with pandas dfs later everything got mixed up when string was involved).
Format: 
> D-GLUCOSE 6-PHOSPHATE: 1
>
> D-FRUCTOSE 6-PHOSPHATE: 2 

The reaction number file contains the metabolic reaction information extracted from the provided metabolic reaction list. As with the previous file, every reaction is assigned a number and all metabolites are represented using their metabolite number to avoid running into string handling [problems. 
Format:
> reactions,reactants,products
>
> 1,[1],[2]
>
> 2,"[2, 3]","[11, 4, 16]"

### Part 1 - The reverse simulation
The reverse simulation (starting with the full JCVI reaction network and deleting reactions with every iteration to see how quickly the network disappears) was only run locally. It is small enough to not be run on the HPC. 
#### Dependencies

#### Running instructions

#### Expected output

### Part 2- The proto-metabolic reaction networks

#### Dependencies

#### Running instructions

#### Expected output

### Part 3- Visualizing the output files


## Usage

## Data and dependencies

### Data sources

### Environment
Requirements files

### HPC configurations
My OBS directives and a sample script, I suppose?

## Licence

## Contact and Contributions
