# Let's try something random
import random
from deap import creator, base, tools, algorithms
import pubchempy as pcp
import re
import pandas as pd
import csv
import numpy as np
from scipy import stats
from data_preprocessing import chem_dictionary

# will technically get this from preprocessing but i'm in a rush and it takes too long to wait for the api
# filepath1 = "central_met_SIE.csv"
# columns = ["Met_id", "Metabolite", "SIE"]
# all_data = pd.read_csv(filepath1, usecols=columns)
# met_id = all_data["Met_id"].values.tolist()
# met_name = all_data["Metabolite"].values.tolist()
# met_SIE = all_data["SIE"].values.tolist()


class Fitness():


    def __init__(self):
        pass

    def fitness_complexity(self,system,amt_cats_added):
        chains = 0
        nodes = 0
        sum_chains = 0
        amt_cats_in_sys = 0
        for chain in system:
            amt_cats_in_sys += 1
            sum_chains += len(chain)
            if len(chain) >= 2:
                chains += 15
        # then maybe normalizing it by dividing through the amount of reactions
        amt_rxn_chains = chains/amt_cats_added
        avg_length_rxn_chain = sum_chains/amt_cats_in_sys
        return amt_rxn_chains, avg_length_rxn_chain

    def shannons_info_ent(self, met_smiles):
        # using the shannon's information entropy to calculate complexity
        # NOW to calculate the Shannon's information entropy just using the characters
        sie = []
        # this is for a single smile
        for smiles in met_smiles:
            if smiles == 0 or smiles == 'NaN':
                sie.append(0)
            else:
                # get the value to divide by
                total_char = len(smiles)
                # turn the smiles str into a smiles list
                structures = list(smiles)
                # count the unique ones in this -single- smile
                char_ft = []
                unique_char = set(structures)
                for char in unique_char:
                    current_char_ft = structures.count(char)
                    char_ft.append([char, current_char_ft])
                # calculate the info entropy
                prob_list = []
                for item in char_ft:
                    # getting the prob
                    prob = item[1] / total_char
                    prob_list.append(prob)
                # calc shannons info ent
                temp_sie = stats.entropy(prob_list, base=2)
                sie.append(temp_sie)
        return sie

    # leaving this to show i thought about this, but it isn't a good indication
    def complexity_smiles_len(self, met_smiles):
        met_lens = []
        for met in met_smiles:
            met_lens.append(len(met))
        return met_lens

    # the idea is to calculate (entropy of products) - (entropy of reactants)
    # and from there whether the fitness increased from the previous state
    def relative_entropy(self, met_smiles, system, rxn_df):
        # calculate the change in entropy for every rxn chain
        #for rxn in system:
        pass


    # def preliminary_fitness_function(self, new_system):
    #     for item in pore:



# main code to test if this works

evaluate = Fitness()

# reading CSV file
all_smiles = open("updated smiles.csv", 'r')
# creating dictreader object
file = csv.DictReader(all_smiles)
met_names = []
ismiles = []

for col in file:
    met_names.append(col['Name'])
    ismiles.append(col['Isomeric Smiles'])

most_sie = evaluate.shannons_info_ent(ismiles)

# for i in range(0, len(ismiles)):
#     print(ismiles[i])
#     print(most_sie[i])
