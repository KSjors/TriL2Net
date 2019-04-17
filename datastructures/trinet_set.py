# -*- coding: utf-8 -*-
"""
Created on Tue Feb 26 14:11:45 2019

@author: Sjors
"""

import logging
from utils.help_functions import guid
from datastructures.rooted_level_k_network import *
import itertools


class TrinetSet:
    def __init__(self, trinet_dict: dict, taxa_names: bidict):
        self.uid = guid()
        self.logger = logging.getLogger("trinet_set.{}".format(self.uid))
        self.trinet_dict = trinet_dict
        self.taxa_names = taxa_names
        self.number_of_taxa = len(taxa_names)

    @classmethod
    def copy_trinet_set(cls, trinet_set):
        logging.debug("Copying trinet_set {}.".format(trinet_set.uid))
        trinet_set.logger.debug("Being copied.")
        trinet_dict = copy.deepcopy(trinet_set.trinet_dict)
        taxa_names = copy.deepcopy(trinet_set.taxa_names)
        trinet_set_copy = cls(trinet_dict, taxa_names)
        trinet_set_copy.logger.debug("Copy from {}.".format(trinet_set.uid))
        return trinet_set_copy

    @classmethod
    def from_trinet_list(cls, trinet_list: list):
        logging.debug("Creating trinet set from triplet trinet list.")
        trinet_dict = {}
        taxa_names = bidict()
        trinet_set = cls(trinet_dict, taxa_names)
        for trinet in trinet_list:
            trinet_set.add_trinet(trinet)

        trinet_set.logger.debug("Created from triplet trinet list.")
        return trinet_set

    def add_trinet(self, trinet: RootedLevelKNetwork):
        """Add trinet to dictionary."""
        triplet = list(trinet.leaf_names)
        self.logger.debug("Adding trinet {} with leaves {}.".format(trinet.uid, triplet))
        assert len(triplet) == 3, "Can not add trinet {} to trinet_set {} as it does not have exactly three leaves.".format(trinet.uid, self.uid)
        triplet.sort()
        assert str(triplet) not in self.trinet_dict, "Can not add trinet {} to trinet_set {} as it already exists.".format(trinet.uid, self.uid)
        self.trinet_dict[str(triplet)] = trinet
        for X in triplet:
            if X not in self.taxa_names.keys():
                self.taxa_names[X] = self.number_of_taxa
                self.number_of_taxa += 1

    def remove_trinet(self, triplet: list):
        """Remove trinet with taxa = triplet from dictionary."""
        self.logger.debug("Removing triplet {}.".format(triplet))
        triplet.sort()
        assert len(triplet) == 3, "Can not remove {} as it is not a triplet.".format(triplet)
        trinet = self.trinet_dict[str(triplet)]
        self.logger.debug("Triplet {} corresponds to trinet {}.".format(triplet, trinet.uid))
        self.trinet_dict.pop(str(triplet))

    def cut_arc_sets_per_triplet(self) -> dict:
        """Calculate cut-arc sets of each trinet."""
        self.logger.debug("Calculating cut-arc sets per trinet.")
        cut_arc_sets = {}
        for triplet in self.trinet_dict:
            cut_arc_sets[triplet] = self.trinet_dict[triplet].cut_arc_sets()
        return cut_arc_sets

    def suppress_minimal_sink_set(self, mss: list) -> str:
        """Suppress minimal sink set."""
        self.logger.debug("Suppressing minimal sink set {}.".format(mss))
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
                trinet_network.to_standard_form()
                new_trinet = list(binet) + [new_name]
                new_trinet.sort()
                self.trinet_dict[str(new_trinet)] = trinet_network

        for leaf in mss:
            self.remove_leaf_name(leaf)

        self.add_leaf_name(new_name)
        return new_name

    def remove_leaf_name(self, node_name: str) -> int:
        """Remove leaf name from taxa_names and lower all numbers of names above by one"""
        self.logger.debug("Removing leaf name {}.".format(node_name))
        node_number = self.taxa_names.pop(node_name)
        for y in range(node_number + 1, self.number_of_taxa):
            self.taxa_names[self.taxa_names.inverse[y]] -= 1
        self.number_of_taxa -= 1
        return node_number

    def add_leaf_name(self, node_name: str) -> str:
        """Add leaf name."""
        self.logger.debug("Adding leaf name {}.".format(node_name))
        new_number = self.number_of_taxa
        self.taxa_names.put(node_name, new_number)
        self.number_of_taxa += 1
        return node_name

    def __str__(self) -> str:
        return str(self.taxa_names)

    def __getstate__(self):
        self.logger = 'trinet_set.{}'.format(self.uid)
        result = copy.deepcopy(self.__dict__)
        self.logger = logging.getLogger(self.logger)
        return result

    def __setstate__(self, d):
        self.__dict__ = d
        self.logger = logging.getLogger(self.logger)
        return self.__dict__

