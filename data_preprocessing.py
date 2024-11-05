import pandas as pd
import re
import csv


class JcviReactionsInput():

    def __init__(self):
        pass

    def unpacking_jcvi(self, filepath1):

        """input: filepath1 = the csv file containing the reaction numbers and reactions
        output: reaction_nr = a list with the reaction numbers
                reaction_details = a list with the reaction information"""

        columns = ["Reaction number", "Reaction equation"]
        all_data = pd.read_csv(filepath1, usecols=columns)
        reaction_nr = all_data["Reaction number"].values.tolist()
        reaction_details = all_data["Reaction equation"].values.tolist()
        return(reaction_nr, reaction_details)


class ReactionsDetails():

    def __init__(self):
        return

    def breakdown(self, reactions):

        """input: reactions = a list of reactions without reaction numbers
            output: caps_individual_reactants = uppercase list of all reactant components in reactions list
                    caps_individual_products = uppercase list of all product components in reactions list"""

        components = []
        more_components = []
        more_pieces = []
        reactant_list = []
        product_list = []
        reaction_dict = {}
        product_dict = {}

        # removing the weird stuff in the equation
        for reaction in reactions:
            if ":" in reaction:
                a, b = reaction.split(":")
                components.append(b)
            else:
                components.append(reaction)

        # dividing into reagents and products
        for item in components:
            if '<==>' in item:
                new = item.replace('<==>', ';')
                more_components.append(new)
            elif '-->' in item:
                new = item.replace('-->', ';')
                more_components.append(new)
            elif '<--' in item:
                new = item.replace('<--', ';')
                more_components.append(new)
            else:
                more_components.append(item)

        # now splitting into ";" separated components
        # producing a list of comma separated strings- before the comma is reactants and after products
        for item in more_components:
            comp = item.split(';')
            more_pieces.append(comp)
        i = 0
        for items in more_pieces:
            # just the normal things
            reactants = more_pieces[i][0]
            reactant_list.append(reactants)

            products = more_pieces[i][1]
            product_list.append(products)
            i += 1

        # Initialize an empty list to store the individual reactants
        individual_reactants = []

        # Iterate through each item in the reactant_list
        for item in reactant_list:
            # Split the item at '+' to get a list of substrates
            substrates = item.split('+')

            # Remove leading and trailing whitespace from each substrate and add to the result list
            individual_reactants.extend([substrate.strip() for substrate in substrates])
            # Adding a character so I can figure out which reactants belong to which reaction
            individual_reactants.extend('*')

        # trying to remove only the values in parentheses
        individual_reactants = [re.sub(r'\(\d+(\.\d+)?\)', '', item).strip() for item in individual_reactants]
        individual_reactants = [re.sub(r'\[.*?]', '', item).strip() for item in individual_reactants]

        # Filter out any empty strings
        individual_reactants = list(filter(None, individual_reactants))
        # make it all caps, i'm done struggling
        caps_individual_reactants = [item.upper() for item in individual_reactants]

        individual_products = []

        # Iterate through each item in the reactant_list
        for item in product_list:

            # Split the item at '+' to get a list of substrates
            resulting_products = item.split('+')

            # Remove leading and trailing whitespace from each substrate and add to the result list
            individual_products.extend([resulting_product.strip() for resulting_product in resulting_products])
            individual_products.extend('*')

        # trying to remove only the values in parentheses
        individual_products = [re.sub(r'\(\d+(\.\d+)?\)', '', item).strip() for item in individual_products]
        individual_products = [re.sub(r'\[.*?]', '', item).strip() for item in individual_products]

        # Filter out any empty strings
        individual_products = list(filter(None, individual_products))
        # make it all caps, i'm done struggling
        caps_individual_products = [item.upper() for item in individual_products]

        return(caps_individual_reactants, caps_individual_products)

    def chemical_dictionaries(self, chemicals):

        """input: chemicals = a list of metabolites
            output: chemical_dict = key and value pair of metabolites and an assigned number"""

        value = range(1, len(chemicals))
        i = 0
        chemical_dict = {}
        for chemical in chemicals:
            if chemical == '*':
                pass
            elif chemical in chemical_dict:
                pass
            else:
                chemical_dict[chemical] = value[i]
                i += 1
        return chemical_dict

    def chem_lists(self, chem_dictionary, reactants, products):

        """ input: chem_dictionary = dictionary with metabolite numbers and names
                   reactants = all caps reactant list
                   products = all caps product list
            output: abstracted_rrxn = list of rxn numbers for the reactants
                    abstracted_prxn = list of rxn numbers for the products"""

        abstracted_rrxn = []
        for i in range(0, len(reactants)):
            compound = reactants[i]
            if compound in chem_dictionary:
                abstracted_rrxn.append(chem_dictionary[compound])
            else:
                abstracted_rrxn.append(compound)
        abstracted_prxn = []
        for i in range(0, len(products)):
            compound = products[i]
            if compound in chem_dictionary:
                abstracted_prxn.append(chem_dictionary[compound])
            else:
                abstracted_prxn.append(compound)

        return abstracted_rrxn, abstracted_prxn

    def reaction_dict(self, num_rxn):
        """ input: num_rxn= reactant or product list
            output: dictionary with the reaction number and the corresponding reactants and products"""
        num_dict = {}
        num_index = 1
        current_reactants = []

        # Iterate through list1 to build the reactants dictionary
        for item in num_rxn:
            if item == '*':
                # End of a reaction, add to dictionary and reset reactants
                num_dict[num_index] = current_reactants
                num_index += 1
                current_reactants = []
            else:
                # Add reactant to the current list
                current_reactants.append(item)

        # Add the last reaction
        if current_reactants:
            num_dict[num_index] = current_reactants
        return num_dict

    # since we are so good with this let's create one more
    def a_sort_of_lookup_table(self, r_dict, p_dict):

        """ input: r_dict = reactant number and reactant details
                  p_dict = product number and product details
            output: lookup_df = pd df with reaction number and corresponding metabolite numbers"""

        lookup_df = pd.DataFrame({'reactants': r_dict.values(), 'products': p_dict.values()}, index=r_dict.keys())
        return lookup_df



# Class ONE
reactions = JcviReactionsInput()
# lists jcvi data
j_reaction_nr, j_reactions = reactions.unpacking_jcvi("C:\\PycharmProjects\\APMES_4.0\\full_jcvi_rxn_list.csv")

# Class TWO
data = ReactionsDetails()
# extracting the relevant info from the csv
reactants, products = data.breakdown(j_reactions)
# compiling a metabolite list
total_chemicals = reactants + products
# key and value pairs for all metabolites
chem_dictionary = data.chemical_dictionaries(total_chemicals)
# optionally writing it to a file if we uncomment this section
with open('chem_dict.txt', 'w') as f:
    for key, value in chem_dictionary.items():
        f.write(f"{key}: {value}\n")
# uniformity otherwise i get string/int handling errors later on
rrxn, prxn = data.chem_lists(chem_dictionary, reactants, products)
# just turning it into dfs for future use
r_dict = data.reaction_dict(rrxn)
p_dict = data.reaction_dict(prxn)
rxn_df = data.a_sort_of_lookup_table(r_dict, p_dict)
rxn_df = rxn_df.rename_axis('reactions')
rxn_df.to_csv('rxn numbers.csv')