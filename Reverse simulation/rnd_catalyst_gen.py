import random
import numpy as np
import pandas as pd
random.seed()


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

    # logic: after doing a statistical analysis on the jcvi central, lipid and nucleotide metabolism data sets, the
    # conclusion was drawn that some reactions were present in multiple subsection
    # following research by NASA on the ‘biased typewriter’ model, it is observed that some molecules and chemical
    # reactions are more likely to occur than others which means to accurately represent the reality of chemistry in
    # HT vents, I would have to introduce a bit of bias in the random reaction function
    # if you have a process that generates these monomers at the right frequency, then you’re going to be able to find
    # the self-replicators much faster
    def but_make_it_biased_towards_selfrep(self, amount):
        return


main_cat = CatalystGenerator()