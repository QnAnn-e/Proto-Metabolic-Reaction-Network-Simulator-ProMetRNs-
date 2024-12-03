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

```
Reaction nr,Reaction equation

S001,D-glucose 6-phosphate <==> D-fructose 6-phosphate

S002,D-fructose 6-phosphate + ATP --> ADP + D-fructose 1,6-biphosphate + H+
```

#### Expected output
Two files will be created: 
1. chem_dict.txt
2. rxn numbers.csv

The chemical dictionary contains all of the metabolites extracted from the metabolic reaction list provided and its assigned value (I had to assign values because when working with pandas dfs later everything got mixed up when string was involved).
Format: 

```
D-GLUCOSE 6-PHOSPHATE: 1

D-FRUCTOSE 6-PHOSPHATE: 2 
```

The reaction number file contains the metabolic reaction information extracted from the provided metabolic reaction list. As with the previous file, every reaction is assigned a number and all metabolites are represented using their metabolite number to avoid running into string handling [problems. 
Format:

```
reactions,reactants,products

1,[1],[2]

2,"[2, 3]","[11, 4, 16]"
```


### Part 1 - The reverse simulation
The reverse simulation (starting with the full JCVI reaction network and deleting reactions with every iteration to see how quickly the network disappears) was only run locally. It is small enough to be run locally. 

#### Dependencies

Makes use of the following python modules: 
* matplotlib
* networkx
* numpy
* pandas
* csv
* random
* concurrent.futures

#### Running instructions
1. pip install all of the above modules
2. perform the data preprocessing on the [full_jcvi_rxn_list.csv](Data%20files/full_jcvi_rxn_list.csv) by making use of [data_preprocessing.py](data_preprocessing.py)
3. Add the output files from the data preprocessing, together with [updated smiles.csv](Data%20files/updated_smiles.csv) to the same directory and run the [reverse simulation](Reverse%20simulation/reverse_sim2.py) script


#### Expected output

The script should produce: 
* A "systemfitness.csv" file that will contain the fitness values for all systems after every 5 iterations
* A "full.svg" image file that contains the full JCVI chemical reaction network overview before any of the deletions started happening
* A "reverse_sim_network.gml" file that contains a key and value list that describes the full reaction network.  



### Part 2- The proto-metabolic reaction networks

#### Dependencies

All dependancies for this program can be found in the [requirements](HPC%20setup%20files/requirements.txt) file. The files in main_15000_4000 is used to run a simulation for 15 000 systems over 4 000 iterations and the files in main_60000_8000 to run a simulation for 60 000 systems over 8 000 iterations. 

#### Running instructions

This simulation is very computationally expensive and requires the use of a high perfomance computing (HPC) unit/cluster.

1. Use the same "chem_dict.txt" and "rxn numbers.csv" files from the reverse simulation produced by the [data preprocessing](data_preprocessing.py) script
2. Create a new directory and add the files in [HPC setup files](HPC%20setup%20files) and [Main simulation](Main%simulation)
3. Change the PBS directives in the [simulation](Main%20simulation/simulation2.pbs) file to fit the capacity of the HPC you are using and then submit the script


#### Expected output

The script will generate the following:
1. 4 cythonized files - 2 for the more_cyhtonized.py file and 2 for the fastest_cython.py file
2. A "fitness.csv" file containing the fitnesses for all systems over all iterations
3. A "system.csv" file containing the system information of all systems
4. An output, error and log file to inform if anything out of the ordinary happened during the run


### Part 3- Visualizing the output files
The network visualization can be done by running the "vis1.pbs" script in [Data visualization and analysis](Data%20visualization%20and%20analysis/). The fitness and standard deviation graphs can be generated by running the "vis1.pbs" script in [Fitness Calculation](Fitness%20calculation/). 

## Data and dependencies

### Data sources 

The data used to compile the [JCVI reaction list](Data%20files/full_jcvi_rxn_list.csv) was obtained from:
Marian Breuer, Tyler M. Earnest, Chuck Merryman, Kim S. Wise, Lijie Sun, Michaela R. Lynott, Clyde A. Hutchison, Hamilton O. Smith, John D. Lapek, David J. Gonzalez, Valérie de Crécy-Lagard, Drago Haas, Andrew D. Hanson, Piyush Labhsetwar, John I. Glass, Zaida Luthey-Schulten (2019) Essential metabolism for a minimal cell eLife 8:e36842. 


### Environment
The virtual environment is automatically setup when running the [simulation](Main%20simulation/simulation2.pbs) script. 
It makes use of the [requirements](HPC%20setup%20files/requirements.txt) file together with the "source sim2env/bin/activate" command. 


### HPC configurations

Sample PBS directives:
```
#!/bin/bash
#PBS -N simulation
#PBS -l select=1:ncpus=80:mem=80GB:vnode=n12.hpc
#PBS -l walltime=160:00:00
#PBS -o sim.out
#PBS -e sim.err
#PBS -m abe
#PBS -M youremailaddress
```

## Licence
MIT LICENSE

## Contact

Whether you have questions, ideas, or just want to grab a coffee and chat about the Origin of Life, I’d love to hear from you! You can reach me at:

Email: 26703173@sun.ac.za
