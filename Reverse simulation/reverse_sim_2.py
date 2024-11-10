"""reverse simulation for prof
this is where we add all the reactions to the system and take them away at random to see how long it survives....
what is the conditions for survival? probably just checking to see how long before all the chains disappear? """

# existing python modules
import matplotlib.pyplot as plt
import networkx as nx
import itertools
from numpy import genfromtxt
import pandas as pd
import csv
import copy
import random
# import multiprocessing as mp can't run multiple executions of the same process, only different processes in parallel
import time
import concurrent.futures

# my own classes
# this gets me the jcvi reactions in a workable format
from data_preprocessing import chem_dictionary, rxn_df
# this is the method used to keep everything random
from rnd_catalyst_gen import main_cat
# check whether a reaction chain can exist
from chemical_rules import chemical_rules
# to draw the graph
from graphs import instantiate
# the sie function
from fitness import evaluate


random.seed()

# reading CSV file
all_smiles = open("updated smiles.csv", 'r')
# creating dictreader object
file = csv.DictReader(all_smiles)
met_names = []
ismiles = []

for col in file:
    met_names.append(col['Name'])
    ismiles.append(col['Isomeric Smiles'])

# and calculating that SIE
all_sie = evaluate.shannons_info_ent(ismiles)

# getting the met id's for the metabolites
met_id = []
for metabolite in met_names:
    if metabolite in chem_dictionary:
        current_met_id = chem_dictionary[metabolite]
        met_id.append(current_met_id)

# getting the met id's for the metabolites
met_id = []
for metabolite in met_names:
    if metabolite in chem_dictionary:
        current_met_id = chem_dictionary[metabolite]
        met_id.append(current_met_id)


# now just creating a small new dict
sie_dict = {met_id[i]: all_sie[i] for i in range(len(met_id))}


class ReverseSim():

    def __init__(self):
        pass

    def network_stats(self, graph_f):
        """Calculate weakly connected components --You treat the graph as undirected to check how nodes are connected
        regardless of edge direction.
        if > 1 graph has multiple disjoint parts when treated as an undirected graph
        Each weakly connected component represents a separate "island" of nodes that can be navigated between,
        but there is no path between nodes in different weakly connected components."""
        weakly_connected_subgraphs = list(nx.weakly_connected_components(graph_f))
        num_weakly_connected_subgraphs = len(weakly_connected_subgraphs)

        """Calculate strongly connected components-check if all nodes within each component are reachable from each other
        via directed paths, which is more restrictive and specific to directed graphs.
        if > 1 not all nodes are reachable from each other when considering edge direction
        each strongly connected component is a subgraph where you can navigate between all nodes, but you cannot navigate
        between nodes in different strongly connected components when following the direction of the edges"""
        strongly_connected_subgraphs = list(nx.strongly_connected_components(graph_f))
        num_strongly_connected_subgraphs = len(strongly_connected_subgraphs)
        total_dc, average_degree, bipartite_clustering, density, assortativity = instantiate.graph_metrics(graph_f)
        return

    def full_system(self):
        # setting up the full system
        first_batch = list(range(1, 337))
        jcvi_in_sys = chemical_rules.map_to_jcvi(first_batch)
        filtered_df = chemical_rules.selected_reactions(rxn_df, jcvi_in_sys)
        reaction_chains = chemical_rules.verifying_catalysts(filtered_df, jcvi_in_sys)
        return reaction_chains

    def fitness_calc(self, m_nodes, r_nodes, s_edges):
        # approach 1 - using the SIE as an indicator of fitness
        # part 1 -  get the smiles for all the metabolites in the system
        central_sie = []
        for item in m_nodes:
            if item in met_id:
                # get the metabolite SIE
                central_sie.append(sie_dict[item])
            else:
                pass

        # part 2 - summing the SIE
        system_sie = sum(central_sie)

        # part 3 - subtracting the energy
        system_fitness = system_sie - len(m_nodes)
        return system_fitness

    def network_info(self, reaction_chains):
        # stats at iteration 0
        r_nodes, m_nodes = instantiate.getting_the_nodes(reaction_chains)
        s_edges = instantiate.creating_edges(r_nodes, m_nodes)
        graph_f = instantiate.creating_a_graph(r_nodes, m_nodes, s_edges, "full.svg")
        number_of_edges = graph_f.number_of_edges()
        sys_fit = reverse_sim.fitness_calc(m_nodes, r_nodes, s_edges)
        return graph_f, sys_fit


reverse_sim = ReverseSim()
# gml_file = 'reverse_sim_network.gml'
# nx.write_gml(graph_r, gml_file)
#plt.show()

def process_system(j):
    system_fitness = []
    system_fitness.append(j)
    # print(f"system {j}")
    connections = reverse_sim.full_system()
    graph_r = reverse_sim.network_info(connections)
    # reverse_sim.network_stats(graph_r)

    for i in range(0, (len(connections))):
        new_connections = connections
        # random item is [1,2] literally a connection between two nodes being deleted
        random_item = random.choice(new_connections)
        new_connections.remove(random_item)
        if i % 5 == 0:
            graph_r, fitness = reverse_sim.network_info(new_connections)
            system_fitness.append(fitness)
            # reverse_sim.network_stats(graph_r)
    return system_fitness

if __name__ == '__main__':
    with concurrent.futures.ProcessPoolExecutor() as executor:
        system_fitness = executor.map(process_system, range(100))

    # write the system fitness to a file
    with open('system_fitnesses.csv', 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerows(system_fitness)

# instantiate.different_way(graph_r)