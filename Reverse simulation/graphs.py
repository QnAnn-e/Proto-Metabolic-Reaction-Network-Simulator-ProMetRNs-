import networkx as nx
from networkx.algorithms import bipartite, community
import matplotlib.pyplot as plt
import pandas as pd
from data_preprocessing import rxn_df
import csv
import heapq
import numpy as np
import community as community_louvain


class VisualizeMyDats():

    def __init__(self):
        pass

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
                    r_met = (rxn_df.loc[rxn, 'reactants'])
                    # getting the products
                    p_met = (rxn_df.loc[rxn, 'products'])
                    # combine to make sure we don't have doubles
                    metabolites = r_met + p_met
                    for met in metabolites:
                        if met not in met_list:
                            met_list.append(met)
        return rxn_list, met_list

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

    def creating_a_graph(self, r_nodes, m_nodes, s_edges, graph_name):

        """ input: r nodes = list with reaction nodes
                   m nodes = list with metabolite nodes
                   s edges = list of lists showing all edges
                   graph name = custom graph name to save the figs under
            output: G = the graph"""

        # Creating the empty bipartite canvas
        G = nx.MultiDiGraph()

        # Transforming the data here, otherwise everything crashes
        # just one last thing before we go
        new_rxn = []
        for rxn in r_nodes:
            new_rxn.append(f'r' + str(rxn))

        # Adding the system node data
        G.add_nodes_from(new_rxn, bipartite=0)  # R
        G.add_nodes_from(m_nodes, bipartite=1)  # M

        # add edges
        G.add_edges_from(s_edges)
        # two node sets U and V
        R = set(new_rxn)
        M = set(m_nodes)
        # check that it is bipartite
        bipartite.is_bipartite(G)

        # Extract node attributes for coloring and shaping nodes
        #node_color = [G.nodes[n]['colour'] for n in G.nodes()]
        #node_shape = [G.nodes[n]['node_shape'] for n in G.nodes()]

        # pos = nx.spring_layout(G, scale=1.5*len(m_nodes), iterations=35, k=0.75)
        # pos = nx.spectral_layout(G, scale=1.5 * len(m_nodes))
        # trying the bipartite way -- it's not working
        # pos = nx.bipartite_layout(G, R)
        # nx.draw(G, pos, with_labels=True)

        # show our graph
        # nx.draw_networkx(G, node_color=node_color, pos=pos,  font_size=6, node_size=150, arrowsize=5, width=0.7)
        # Compute the Kamada-Kawai layout
        pos = nx.kamada_kawai_layout(G)

        # Draw the graph with different colors for sets U and V
        node_color = ['skyblue' if G.nodes[node]['bipartite'] == 0 else 'yellow' for node in G.nodes]

        # Plot the graph, ensuring no overlap
        plt.figure(figsize=(14, 7))
        # Avoiding node overlap by adjusting the scale
        pos_no_overlap = nx.rescale_layout_dict(pos, scale=7500)
        nx.draw(G, pos_no_overlap, with_labels=True, node_color=node_color, node_size=60, font_size=4)
        plt.savefig(graph_name)
        plt.clf()
        plt.close()
        # plt.show()
        return G

    def different_way(self, B):
        # Extract nodes of each bipartite set
        top_nodes = {n for n, d in B.nodes(data=True) if d['bipartite'] == 0}
        bottom_nodes = set(B) - top_nodes
        # Bipartite layout
        pos = nx.bipartite_layout(B, top_nodes)
        # Draw the graph with positions specified by bipartite_layout
        nx.draw(B, pos, with_labels=True, node_color=["skyblue" if n in top_nodes else "lightgreen" for n in B])
        # plt.show()
        return

    def graph_metrics(self, G):

        """ input: G = the graph
            output: total_dc = value representing the measure of degree centrality -  how well-connected the network is
                    average_degree = value representing how many connections (edges) each node (vertex) has on average
                    bipartite_clustering = value representing the locally interconnectivity of the nodes
                    density = value showing how close the network is to being fully connected
                    assortativity = value that provides insight into information flow through a system, lower value,
                    faster flow through network"""

        # calculating degree centrality
        total_dc = 0
        degree_centrality = nx.degree_centrality(G)
        for dc in degree_centrality.values():
            total_dc += dc
        # print("total degree centrality:", total_dc)

        # Calculate metrics
        degree_distribution = [d for n, d in G.degree()]
        average_degree = sum(degree_distribution) / len(degree_distribution)
        # print("Average Degree:", sum(degree_distribution) / len(degree_distribution))
        # avg clustering does not work on this type of graph
        bipartite_clustering = nx.bipartite.clustering(G)
        #print("Average Clustering Coefficient:", bipartite_clustering)
        # average path length
        # if nx.is_strongly_connected(G):
        #     avg_path_length = nx.average_shortest_path_length(G)
        #     print("Average Shortest Path Length:", avg_path_length)
        #     diameter = nx.diameter(G)
        #     print("Diameter:", diameter)
        # else:
        #     print("The graph is not connected enough.")

        density = nx.density(G)
        # print("Graph Density:", density)

        assortativity = nx.degree_pearson_correlation_coefficient(G)
        # print("Degree Assortativity:", assortativity)
        return total_dc, average_degree, bipartite_clustering, density, assortativity

    def nmfa(self, G):
        # calculating node properties
        node_properties_U = {node: G.degree(node) for node in G.nodes if G.nodes[node]['bipartite'] == 0}
        node_properties_V = {node: G.degree(node) for node in G.nodes if G.nodes[node]['bipartite'] == 1}

        # print("Node properties (U):", node_properties_U)
        # print("Node properties (V):", node_properties_V)

        # Choose a range of q values
        q_U_values = np.linspace(-10, 10, 50)
        # Scale parameter (can vary depending on the network size)
        epsilon = 0.1
        Z_q_U = []
        # Loop over each q value to compute Z_q_U
        for q_U in q_U_values:
            partition_sum_U = 0
            for value in node_properties_U.values():
                partition_sum_U += value ** q_U
            Z_q_U.append(partition_sum_U / len(node_properties_U))

        # Convert Z_q_U to a numpy array
        Z_q_U = np.array(Z_q_U)
        # Calculate tau_q_U after the loop
        tau_q_U = np.log(Z_q_U) / np.log(1 / epsilon)
        # Calculate D_q_U outside the loop
        D_q_U = tau_q_U / (q_U_values - 1)
        # Calculate alpha_U and f_alpha_U outside the loop
        alpha_U = np.gradient(D_q_U, q_U_values)
        f_alpha_U = q_U_values * alpha_U - D_q_U

        Z_q_V = []
        q_V_values = np.linspace(-10, 10, 50)

        for q_V in q_V_values:
            partition_sum_V = 0
            for value in node_properties_V.values():
                partition_sum_V += value ** q_V
            Z_q_V.append(partition_sum_V / len(node_properties_V))

        Z_q_V = np.array(Z_q_V)
        tau_q_V = np.log(Z_q_V) / np.log(1 / epsilon)
        D_q_V = tau_q_V / (q_V_values - 1)
        alpha_V = np.gradient(D_q_V, q_V_values)
        f_alpha_V = q_V_values * alpha_V - D_q_V

        print("alpha U:", alpha_U)
        plt.plot(alpha_U, f_alpha_U, label='Set U')
        print("alpha V:", alpha_V)
        plt.plot(alpha_V, f_alpha_V, label='Set V')
        plt.xlabel('Singularity Strength (alpha)')
        plt.ylabel('Multifractal Spectrum (f(alpha))')
        plt.title('Multifractal Spectrum of the Bipartite Network')
        plt.legend()
        plt.show()

    def dijkstra(self, G):
        pass


instantiate = VisualizeMyDats()
# r_nodes, m_nodes = instantiate.getting_the_nodes()
# print("rxn nodes:", r_nodes)
# print("met nodes:", m_nodes)
# s_edges = instantiate.creating_edges(r_nodes, m_nodes)
# print("edge list:", s_edges)
# graph = instantiate.creating_a_graph(r_nodes, m_nodes, s_edges)


#with open('edgelist.csv', 'w') as f:
#     writer = csv.writer(f)
#     for item in s_edges:
#         writer.writerow(item)

#
# # Get bipartite layout
# pos = nx.bipartite_layout(graph, graph.nodes())
#
# # Separate nodes and edges
# node_set1 = {n for n, d in graph.nodes(data=True) if d['bipartite'] == 'nodes'}  # Nodes
# node_set2 = set(graph) - node_set1  # Edges
#
# # Draw nodes with different shapes
# nx.draw_networkx_nodes(graph, pos, nodelist=node_set1, node_shape='o', node_color='r')  # Nodes
# nx.draw_networkx_nodes(graph, pos, nodelist=node_set2, node_shape='s', node_color='b')  # Edges
#
# # Draw edges
# nx.draw_networkx_edges(graph, pos)
#
# # Add labels to nodes
# labels = {node: node for node in graph.nodes()}
# nx.draw_networkx_labels(graph, pos, labels=labels)
#
# # Show the plot
# plt.show()

# graph= instantiate.setting_the_scene(new_system)
#
#
# pos = nx.spring_layout(graph)
# nx.draw_networkx(graph, pos)

# Draw nodes of set1 (square)
#nx.draw_networkx_nodes(graph, pos, nodelist=list(rxn), node_shape='s', node_color='lightblue')

# Draw nodes of set2 (circle)
#nx.draw_networkx_nodes(graph, pos, nodelist=list(met), node_shape='o', node_color='orange')


#plt.show()

# def getting_them_nodes_and_edges(self, new_system, filtered_df):
#     # create graphs for some systems
#     graph_edges = []
#     graph_nodes = []
#     for generation in new_system:
#         for reaction in generation:
#             if reaction == None:
#                 pass
#             else:
#                 graph_edges.append('r' + reaction)
#                 for row in rxn_df:
#                     # iterate through df and add reactant and products
#                     if reaction == rxn_df[]:
#
#     return
#
#
# def setting_the_scene(self, new_system):
#     # gentle reminder: new_system format- [[rxn1, rxn2], [rxn4], [rxn12, rxn13], [rxn45], [rxn18]]
#     # Creating the empty bidirectional canvas
#     G = nx.DiGraph()
#     # Adding the data from the specific system in
#     G.add_edges_from(new_system)
#     return G