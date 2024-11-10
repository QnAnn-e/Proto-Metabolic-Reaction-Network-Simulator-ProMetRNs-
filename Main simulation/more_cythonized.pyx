# cython: language_level=3

"""
The simulation without any additional functions that will be run on the hpc as a baseline simulation
steps:
1. data preprocessing: extracting the information from the reaction list and turning it into a dictionary and df
2. setting up an empty system
3. generating random catalysts
4. verifying whether they are actually in jcvi and adding them if they are
5. repeating this process x times before removing all catalysts that weren't integrated
6. calculating the fitness """

# importing built-in modules
import networkx as nx
import itertools
from numpy import genfromtxt
import pandas as pd
import csv
from csv import writer
import copy
import concurrent.futures
import re
import random
import time
import logging
from logging.handlers import RotatingFileHandler
from scipy import stats
import sys
import os
from networkx.algorithms import bipartite, community
import matplotlib.pyplot as plt
import pickle
import Cython
from fastest_cython import fitness_calc



class CatalystGenerator():
    """it creates a random number that represents a hypothetical catalyst that could be present in a system
    this is it's own class cause the catalysts will undergo some more intense hacking later on"""

    def __init__(self):
        return

    def creating_catalysts(self, amount):

        """ input: amount = number of catalysts to be chosen
            output: rand cat list = list with specified number of random catalysts"""

        # setting the seed outside of the for loop
        rand_cat_list = []
        while amount != 0:
            rand_cat = random.randrange(1, 99999)
            rand_cat_list.append(rand_cat)
            amount -= 1
        return rand_cat_list


class My_Chemical_Rules():

    def __init__(self):
        pass

    def map_to_jcvi(self, system):

        """ input: system = list of randomly picked reactions
            output: jcvi in sys = list of reactions that is in both the random list and in jcvi"""

        # just keep the rnd numbers associated with jcvi rxns
        jcvi_in_sys = []
        for rxn_nr in system:
            if rxn_nr <= 336: ### change here when you work with full reaction set
                jcvi_in_sys.append(rxn_nr)
        return jcvi_in_sys

    def selected_reactions(self, rxn_df, jcvi_in_sys):

        """ input: rxn df = the original reaction dataframe with the reactants and products of all reactions
                   jcvi in sys = list of reactions in system
            output: filtered df = smaller dataframe with only reaction information for reactions currently in system"""

        # I thought I wouldn't need to create a new df, but alas
        filtered_df = rxn_df.loc[jcvi_in_sys]
        return filtered_df

    def verifying_catalysts(self, filtered_df, jcvi_in_sys):

        """ input: filtered df = smaller df with information for only the reactions currently in the system
                   jcvi in sys = list of reactions in system
            output: connections = list of list grouping the growing reaction chains"""

        # Iterate through each reaction in jcvi_in_sys
        for reaction in jcvi_in_sys:
            # Get the products for the current reaction
            current_products = filtered_df.loc[reaction, 'products']

            # Flag to track whether a connection has been added for this reaction
            connection_added = False

            # Iterate through each row in filtered_df to check for connections
            for index, row in filtered_df.iterrows():
                # Skip the current reaction
                if index == reaction:
                    continue

                # Check if any product of the current reaction is in the reactants of other reactions
                if any(product in row['reactants'] for product in current_products):
                    # If a connection has not been added yet, add the connection
                    if not connection_added:
                        # Append a sub-list with the reaction numbers of both reactions involved in the connection
                        connections.append([reaction, index])
                        # Set the flag to True to indicate that a connection has been added
                        connection_added = True
        return connections


class PreFitness():

    def __init__(self):
        pass

    def getting_those_smiles(self):

        """output: lists containing metabolite names and isomeric smiles ids"""

        all_smiles = open("updated smiles.csv", 'r')
        smiles_file = csv.DictReader(all_smiles)
        met_names = []
        ismiles = []
        for col in smiles_file:
            met_names.append(col['Name'])
            ismiles.append(col['Isomeric Smiles'])
        return met_names, ismiles

    def smiles_to_sie(self, ismiles):

        """ input: isomeric smiles id list
            output shannon information entropy values for all isomeric smiles list"""

        sie = []
        for smiles in ismiles:
            if smiles == 0 or smiles == 'NaN':
                sie.append(0)
            else:
                total_char = len(smiles)
                structures = list(smiles)
                char_ft = []
                unique_char = set(structures)
                for char in unique_char:
                    current_char_ft = structures.count(char)
                    char_ft.append([char, current_char_ft])
                prob_list = []
                for item in char_ft:
                    prob = item[1] / total_char
                    prob_list.append(prob)
                temp_sie = stats.entropy(prob_list, base=2)
                sie.append(temp_sie)
        return sie

    def post_processing(self, sie, met_names):

        """input: sie list and metabolite names list
           output: dictionary with the sie values"""
        # getting the met id's for the metabolites
        met_id = []
        for metabolite in met_names:
            if metabolite in chem_dictionary:
                current_met_id = chem_dictionary[metabolite]
                met_id.append(current_met_id)
            else:
                print("oh no!!", metabolite)
        # now just creating a small new dict
        sie_dictionary = {met_id[i]: sie[i] for i in range(len(met_id))}
        return sie_dictionary


class VisualizeMyDats():

    def __init__(self):
        pass

    @Cython.cdivision(True)
    @Cython.boundscheck(False)
    @Cython.wraparound(False)
    @Cython.nonecheck(False)
    def getting_the_nodes(self, new_system):

        """ input: new system = list of lists showing reaction chains currently in system
            output: reaction list = list of reactions in system (node type =0)
                    metabolite list = list of metabolites in system (node type =1)"""
        # node list 1
        rxn_list = []
        # node list 2
        met_list = []
        # iterate through lists of reactions in system
        for sub_sys in new_system:
            if sub_sys != []:
                # iterate through reactions in sub list
                for rxn in sub_sys:
                    # only add if not in list yet
                    if rxn not in rxn_list:
                        rxn_list.append(rxn)
                    # now adding the mets associated with a reaction to the met list
                    # getting the reactants at specific loc
                    r_met = rxn_df.loc[rxn, 'reactants']
                    # getting the products
                    p_met = rxn_df.loc[rxn, 'products']
                    # combine to make sure we don't have doubles
                    metabolites = r_met + p_met
                    for met in metabolites:
                        # print(met)
                        if met not in met_list:
                            met_list.append(met)
        return rxn_list, met_list

    @Cython.cdivision(True)
    @Cython.boundscheck(False)
    @Cython.wraparound(False)
    @Cython.nonecheck(False)
    def creating_edges(self, r_nodes, m_nodes):

        """ input: r nodes = list of reaction nodes
                   m nodes = list of metabolite nodes
            output: s edges = list of lists containing the edges """
        # creating empty list to add system edges to
        s_edges = []
        # iterating through reaction nodes to specify which edges should be added
        for reaction in r_nodes:
            # getting the reactants
            temp_r = rxn_df.loc[reaction, 'reactants']
            # iterating through list if there is more than 1 reactant
            for rxn in temp_r:
                temp_r_edges = ((f'r' + str(reaction)), rxn)
                # adding the new edges to our edge list
                s_edges.append(temp_r_edges)
            # getting the products
            temp_p = rxn_df.loc[reaction, 'products']
            # iterating through list if there is more than 1 reactant
            for rxn in temp_p:
                # this is important cause this is how i tell the program that the arrow should face away from the thing
                temp_p_edges = (rxn, (f'r'+ str(reaction)))
                # adding the new edges to our edge list
                s_edges.append(temp_p_edges)
        return s_edges


class Model():

    def __init__(self):
        pass

    @Cython.cdivision(True)
    @Cython.boundscheck(False)
    @Cython.wraparound(False)
    @Cython.nonecheck(False)
    def setting_up_system(self, system=None):
        """ this creates an empty system and calls the function that can add catalysts
        input: system = list of lists for reactions
        output: list of the reaction chains formed in this vent """

        # step 1: an empty space for the magic to happen
        if system is None:
            system = []
        # step 2: adding rxn by calling function from random class to add randomly chosen catalysts in range (1, 99 999)
        # function format: class.function(amount)
        first_batch = main_cat.creating_catalysts(10)
        # sift out what doesn't map to jcvi
        jcvi_in_sys = chemical_rules.map_to_jcvi(first_batch)
        # filtering duplicates because we don't quantify reactions at this point
        jcvi_in_sys = list(set(jcvi_in_sys))
        if jcvi_in_sys:
            # then we add the reaction numbers to our system
            system.append(jcvi_in_sys)
        # step 2b: now we create a smaller df to iterate through quicker
        # function format class.function(filtered_df, system_content)
        filtered_df = chemical_rules.selected_reactions(rxn_df, jcvi_in_sys)
        # step 3: function format class.function(filtered_df, system_content)
        reaction_chains = chemical_rules.verifying_catalysts(filtered_df, jcvi_in_sys)
        return reaction_chains

    @Cython.cdivision(True)
    @Cython.boundscheck(False)
    @Cython.wraparound(False)
    @Cython.nonecheck(False)
    def death_of_the_weak(self, current_system):
        remaining_reactions = [reaction_chain for reaction_chain in current_system if len(reaction_chain) > 1]
        return remaining_reactions

### MAIN ###

# preprocessing
chem_dictionary = {}
with open("chem_dict.txt", 'r') as file:
    for line in file:
        line = line.strip()
        if not line:
            continue
        key, value = line.split(':', 1)
        chem_dictionary[key.strip()] = int(value.strip())  # Convert value to integer
rxn_df = pd.read_csv('rxn numbers.csv')
rxn_df['reactants'] = rxn_df['reactants'].apply(lambda x: list(map(int, x.strip('[]').split(', '))))
rxn_df['products'] = rxn_df['products'].apply(lambda x: list(map(int, x.strip('[]').split(', '))))

# rnd rxns
main_cat = CatalystGenerator()

# sys rules
chemical_rules = My_Chemical_Rules()
connections = []

# fitness
sys_fitness = PreFitness()
met_id, ismiles = sys_fitness.getting_those_smiles()
sie = sys_fitness.smiles_to_sie(ismiles)
sie_dict = sys_fitness.post_processing(sie, met_id)

# graphs
instantiate = VisualizeMyDats()

# Instantiate class
bio_system = Model()

# Setting up rotating file handler
handler = RotatingFileHandler('simprog.log', maxBytes=10*1024*1024, backupCount=5)  # 10 MB, keep 5 backups
formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
handler.setFormatter(formatter)

# Setting up logging
logging.basicConfig(level=logging.INFO, handlers=[handler])

@Cython.cdivision(True)
@Cython.boundscheck(False)
@Cython.wraparound(False)
@Cython.nonecheck(False)
def process_pore(pore, lifetime_range):
    new_system = []
    temp_sys_fit1 = []
    logging.info(f'Starting system {pore}')
    for lifetime in range(lifetime_range):
        quotient, remainder = divmod(lifetime, 10)

        # Add catalysts to the system
        new_system = bio_system.setting_up_system(new_system)

        # check whether there are some disconnects after ten iterations
        if remainder == 0:
            new_system = bio_system.death_of_the_weak(new_system)

        # Get nodes and edges
        r_nodes, m_nodes = instantiate.getting_the_nodes(new_system)
        s_edges = instantiate.creating_edges(r_nodes, m_nodes)

        # Calculate fitness if there are metabolites
        if m_nodes:
            total_sys_fitness = fitness_calc(m_nodes, sie_dict)
        else:
            total_sys_fitness = 0

        temp_sys_fit1.append(total_sys_fitness)
    # Return both the system and its fitness

    return copy.deepcopy(new_system), copy.deepcopy(temp_sys_fit1)

@Cython.cdivision(True)
@Cython.boundscheck(False)
@Cython.wraparound(False)
@Cython.nonecheck(False)
def create_systems_and_fitness_concurrent(pore_range, lifetime_range):
    random.seed()
    a_small_vent = []
    fit_small_vent = []

    with concurrent.futures.ProcessPoolExecutor() as executor:
        results = executor.map(process_pore, range(pore_range), [lifetime_range] * pore_range)

        # Collect results from all processes
        for system, fitness in results:
            a_small_vent.append(system)
            with open("systems.csv", "a") as f:
                writer_object = writer(f)
                writer_object.writerow(system)
            fit_small_vent.append(fitness)
            with open("fitness.csv", "a") as f:
                writer_object = writer(f)
                writer_object.writerow(fitness)

    return a_small_vent, fit_small_vent


if __name__ == '__main__':
    logging.info('script started')
    #met_id, sie_dict = bio_system.sie()

    # Run the concurrent processing
    a_small_vent, fit_small_vent = create_systems_and_fitness_concurrent(pore_range=100, lifetime_range=100)

