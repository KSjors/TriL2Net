# -*- coding: utf-8 -*-
"""
Created on Tue Feb 26 14:11:45 2019

@author: Sjors
"""

from utils.help_functions import guid
from datastructures.rooted_level_k_network import *
import itertools


class TrinetSet:
    def __init__(self, trinet_dict: dict, taxa_names: bidict):
        logging.debug("Creating a trinet set.")
        self.trinet_dict = trinet_dict
        self.uid = guid()
        logging.debug("Trinet set has uid {}.".format(self.uid))
        self.taxa_names = taxa_names
        self.number_of_taxa = len(taxa_names)

    @classmethod
    def copy_trinet_set(cls, trinet_set):
        trinet_dict = copy.deepcopy(trinet_set.trinet_dict)
        taxa_names = copy.deepcopy(trinet_set.taxa_names)
        return cls(trinet_dict, taxa_names)

    @classmethod
    def from_triplet_trinet_list(cls, triplet_trinet_list):
        trinet_dict = {}
        taxa_names = bidict()
        number_of_taxa = 0
        for triplet_trinet in triplet_trinet_list:
            triplet = list(triplet_trinet[0])
            triplet.sort()
            for X in triplet:
                if X not in taxa_names.keys():
                    taxa_names[X] = number_of_taxa
                    number_of_taxa += 1
            trinet_dict[str(triplet)] = triplet_trinet[1]
        return cls(trinet_dict, taxa_names)

    def cut_arc_sets_per_triplet(self) -> dict:
        logging.debug("Calculating cut-arc sets per trinet.")
        cut_arc_sets = {}
        for triplet in self.trinet_dict:
            cut_arc_sets[triplet] = self.trinet_dict[triplet].cut_arc_sets()
        return cut_arc_sets

    def suppress_minimal_sink_set(self, mss: list) -> str:
        logging.debug("Suppressing minimal sink set {}.".format(mss))
        mss.sort()
        current_taxa_names = copy.deepcopy(list(self.taxa_names))

        # Remove all binets and trinets that are a subset of minimal sink set mss
        if len(mss) > 2:
            trinet_iterator = itertools.combinations(mss, 3)
            for trinet in trinet_iterator:
                trinet = list(trinet)
                trinet.sort()
                self.trinet_dict.pop(str(trinet))

        binet_iterator = itertools.combinations(mss, 2)
        for binet in binet_iterator:
            third_leaf_iterator = itertools.combinations(set(current_taxa_names) - set(mss), 1)
            for third_leaf in third_leaf_iterator:
                trinet = list(binet) + list(third_leaf)
                trinet.sort()
                self.trinet_dict.pop(str(trinet))

        # Replace name of each leaf in mss with combined name
        new_name = '(' + "".join(mss) + ')'
        for leaf in mss:
            binet_iterator = itertools.combinations(set(current_taxa_names) - set(mss), 2)
            for binet in binet_iterator:
                trinet = list(binet) + list([leaf])
                trinet.sort()
                trinet_network = self.trinet_dict.pop(str(trinet))
                trinet_network.rename_node(leaf, new_name)
                new_trinet = list(binet) + [new_name]
                new_trinet.sort()
                self.trinet_dict[str(new_trinet)] = trinet_network

        for leaf in mss:
            self.remove_leaf_name(leaf)

        self.add_leaf_name(new_name)
        return new_name

    def remove_leaf_name(self, node_name: str) -> int:
        node_number = self.taxa_names.pop(node_name)
        for y in range(node_number + 1, self.number_of_taxa):
            self.taxa_names[self.taxa_names.inverse[y]] -= 1
        self.number_of_taxa -= 1
        return node_number

    def add_leaf_name(self, node_name: str) -> str:
        new_number = self.number_of_taxa
        self.taxa_names.put(node_name, new_number)
        self.number_of_taxa += 1
        return node_name

    def __str__(self) -> str:
        return str(self.taxa_names)
