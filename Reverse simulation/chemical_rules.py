import random
import numpy as np
import pandas as pd

# this is the class with all the functions needed to check whether two reactions can form a chain
class My_Chemical_Rules():

    def __init__(self):
        pass

    def map_to_jcvi(self, system):

        """ input: system = list of randomly picked reactions
            output: jcvi in sys = list of reactions that is in both the random list and in jcvi"""

        # just keep the rnd numbers associated with jcvi rxns
        jcvi_in_sys = []
        for rxn_nr in system:
            if rxn_nr <= 337: ### change here when you work with full reaction set
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
        # just deleting all the empty []'s
        #filtered_connections = [sub_list for sub_list in connections if any(isinstance(item, int) for item in sub_list)]
        return connections


chemical_rules = My_Chemical_Rules()
connections = []

# class 2.0 of this chemical rules need to make room for entropy to do its work and following Enrique's suggestion
# should allow for chemical compounds in the system to form chains with chemical compounds seemingly against the rules
# of this system following this logic: ' Similarly, the molecules that were available on early Earth had many different
# ways to bind together to produce a variety of chemical reactions; the chance of nature generating the right molecular
# structure to enable self-replication is slim."
