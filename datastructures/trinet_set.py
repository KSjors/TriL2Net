# -*- coding: utf-8 -*-
"""
Created on Tue Feb 26 14:11:45 2019

@author: Sjors
"""

from datastructures.rooted_level_k_network import *
import itertools


class TrinetSet:
    def __init__(self, trinets):
        self.trinets = {}
        self.taxa_names = bidict()
        self.number_of_taxa = 0
        for trinet in trinets:
            triplet = list(trinet[0])
            triplet.sort()
            for X in triplet:
                if X not in self.taxa_names.keys():
                    self.taxa_names[X] = self.number_of_taxa
                    self.number_of_taxa += 1
            self.trinets[str(triplet)] = trinet[1]

    def cut_arc_sets_per_triplet(self):
        cut_arc_sets = {}
        for triplet in self.trinets:
            cut_arc_sets[triplet] = self.trinets[triplet].cut_arc_sets()
        return cut_arc_sets

    def suppress_minimal_sink_set(self, mss):
        mss.sort()
        current_taxa_names = copy.deepcopy(list(self.taxa_names))

        # Remove all binets and trinets that are a subset of minimal sink set mss
        if len(mss) > 2:
            trinet_iterator = itertools.combinations(mss, 3)
            for trinet in trinet_iterator:
                trinet = list(trinet)
                trinet.sort()
                self.trinets.pop(str(trinet))

        binet_iterator = itertools.combinations(mss, 2)
        for binet in binet_iterator:
            third_leaf_iterator = itertools.combinations(set(current_taxa_names)-set(mss),1)
            for third_leaf in third_leaf_iterator:
                trinet = list(binet) + list(third_leaf)
                trinet.sort()
                self.trinets.pop(str(trinet))

        # Replace name of each leaf in mss with combined name
        new_name = '(' + "".join(mss) + ')'
        for leaf in mss:
            binet_iterator = itertools.combinations(set(current_taxa_names)-set(mss), 2)
            for binet in binet_iterator:
                trinet = list(binet) + list([leaf])
                trinet.sort()
                trinet_network = self.trinets.pop(str(trinet))
                trinet_network._rename_node(leaf, new_name)
                new_trinet = list(binet) + [new_name]
                new_trinet.sort()
                self.trinets[str(new_trinet)] = trinet_network

        for leaf in mss:
            self.remove_leaf_name(leaf)

        self.add_leaf_name(new_name)
        return new_name

    def remove_leaf_name(self, X):
        x = self.taxa_names.pop(X)
        for y in range(x + 1, self.number_of_taxa):
            self.taxa_names[self.taxa_names.inverse[y]] -= 1
        self.number_of_taxa -= 1
        return x

    def add_leaf_name(self, X):
        new_number = self.number_of_taxa
        self.taxa_names.put(X, new_number)
        self.number_of_taxa += 1
        return X

    def __str__(self):
        return str(self.taxa_names)

