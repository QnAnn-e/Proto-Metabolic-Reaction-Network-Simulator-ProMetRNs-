import csv
from itertools import islice
from collections import deque
import numpy as np
import pandas as pd
import networkx as nx
from networkx.algorithms import bipartite, community
import matplotlib.pyplot as plt
import ast
from collections import Counter
from data_preprocessing import rxn_df


def getting_the_nodes(new_system):
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
                r_met = (rxn_df.loc[rxn, 'reactants'])
                # getting the products
                p_met = (rxn_df.loc[rxn, 'products'])
                # combine to make sure we don't have doubles
                metabolites = r_met + p_met
                for met in metabolites:
                    if met not in met_list:
                        met_list.append(met)
    return rxn_list, met_list


def creating_edges(r_nodes, m_nodes):
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
            temp_p_edges = (rxn, (f'r' + str(reaction)))
            # adding the new edges to our edge list
            s_edges.append(temp_p_edges)
    return s_edges


def creating_a_graph(r_nodes, m_nodes, s_edges, graph_name):

    # just one last thing before we go
    new_rxn = []
    for rxn in r_nodes:
        new_rxn.append(f'r' + str(rxn))

    # Creating the empty bipartite canvas
    G = nx.MultiDiGraph()
    plt.figure(figsize=(35, 20))

    # add edges
    G.add_edges_from(s_edges)
    # Adding the system node data
    G.add_nodes_from(new_rxn, bipartite=0)  # R
    G.add_nodes_from(m_nodes, bipartite=1)  # M

    # two node sets U and V
    R = set(new_rxn)
    # labels = {node: node for node in R}
    M = set(m_nodes)
    # labels.update({node: "" for node in M})

    # Draw the graph with different colors for sets U and V
    node_size = [300 if node in R else 200 for node in G.nodes()]
    node_color = [
        'skyblue' if G.nodes[node]['bipartite'] == 0 else '#FFFF007F'
        for node in G.nodes
    ]

    pos_no_overlap = nx.spring_layout(G, k=15, iterations=650, seed=42)
    # pos_no_overlap = nx.bipartite_layout(G, R, align='horizontal', scale=150)
    # pos_no_overlap = nx.shell_layout(G, scale=1500)
    # pos_no_overlap = nx.circular_layout(G, scale=50)
    #pos_rescaled = nx.rescale_layout_dict(pos_no_overlap, scale=500)
    edge_weights = [data['weight'] for _, _, data in G.edges(data=True)]
    nx.draw_networkx_edges(G, pos_no_overlap, width=[w / (len(edge_weights)) for w in edge_weights],
                           edge_color="black", alpha=0.25)
    nx.draw_networkx_nodes(G, pos_no_overlap, node_color=node_color, node_size=node_size)
    nx.draw_networkx_labels(G, pos_no_overlap, font_size=14, font_color="black")

    plt.savefig(graph_name)
    plt.clf()
    plt.close()
    # plt.show()
    return G


file_path = "systems.csv"
chunk_size = 100000  # Number of lines per chunk
parsed_systems = []

with open(file_path, "r") as f:
    chunk = []
    for line_number, line in enumerate(f, start=1):
        try:
            chunk.append(ast.literal_eval(line.strip()))
        except Exception as e:
           # print(f"Skipping line {line_number} due to error: {e}")
            continue
        # Process the chunk when full
        if len(chunk) == chunk_size:
            parsed_systems.extend(chunk)
            chunk = []  # Clear the chunk

    # Process the remaining lines
    if chunk:
        parsed_systems.extend(chunk)

#print(parsed_systems)
new_parsed = []
for system in parsed_systems:
    system_parsed = []
    for item in system:
        if len(item) > 1:
            clean_item = ast.literal_eval(item.strip())
            system_parsed.append(clean_item)
        else:
            system_parsed.append(item)
    new_parsed.append(system_parsed)
print(len(new_parsed))
# df = pd.read_csv('systems.csv', header=None, names=['system_data'])
# df['parsed_system'] = df['system_data'].apply(ast.literal_eval)

total_rxn = []
total_met = []
total_edges = []
for system in new_parsed:
    rxn_list, met_list = getting_the_nodes(system)
    total_rxn.append(rxn_list)
    total_met.append(met_list)
    edge_list = creating_edges(rxn_list, met_list)
    total_edges.append(edge_list)
# getting all the r nodes across all systems
avg_rxn_list = list(set([item for rxns in total_rxn for item in rxns]))
avg_met_list = list(set([item for mets in total_met for item in mets]))
# getting data on the sub systems
number_of_mets_list = [item for sublist in total_met for item in sublist]
met_counts = Counter(number_of_mets_list)
sorted_met_counts = met_counts.most_common(10)
print(sorted_met_counts)
# getting data on the sub systems
number_of_rxn_list = [item for sublist in total_rxn for item in sublist]
rxn_counts = Counter(number_of_rxn_list)
sorted_rxn_counts = rxn_counts.most_common(10)
print(sorted_rxn_counts)
# creating a weighted edge list
total_edges = [item for ed in total_edges for item in ed]
edges_as_tuples = [tuple(edge) for edge in total_edges]
counts = Counter(edges_as_tuples)
avg_edge_list = [(item[0], item[1], {'weight': count}) for item, count in counts.items()]
# creating the graph
# graphs = creating_a_graph(avg_rxn_list, avg_met_list, avg_edge_list, "average")
# gml_file = 'average_network.gml'
# nx.write_gml(graphs, gml_file)