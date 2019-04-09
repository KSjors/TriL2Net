import numpy as np
from tarjan import tarjan
from graphviz import Digraph
from utils.help_functions import *


class Omega:
    def __init__(self, trinet_set):
        self.adj_matrix = np.zeros((trinet_set.number_of_taxa,  trinet_set.number_of_taxa))
        self.taxa_names = trinet_set.taxa_names
        cut_arc_sets_per_triplet = trinet_set.cut_arc_sets_per_triplet()

        triplets = itertools.combinations(self.taxa_names.keys(), 3)
        for triplet in triplets:
            triplet = [X for X in triplet]
            triplet.sort()
            cut_arc_sets = cut_arc_sets_per_triplet[str(triplet)]
            for cut_arc_set in cut_arc_sets:
                if len(cut_arc_set) == 2:
                    x = self.taxa_names[cut_arc_set[0]]
                    y = self.taxa_names[cut_arc_set[1]]
                    z = self.taxa_names[[Z for Z in triplet if Z not in cut_arc_set][0]]
                    self.adj_matrix[x, z] += 1
                    self.adj_matrix[y, z] += 1
        self.adj_matrix += np.eye(len(self.taxa_names))*len(self.taxa_names)

    def minimal_sink_sets(self, level=0):
        adj_dict = {}
        for x in range(len(self.adj_matrix)):
            ys = set(np.where(self.adj_matrix[x] <= level)[0])
            ys -= set([x])
            adj_dict[x] = ys

        tj = tarjan(adj_dict)
        return list_of_list_to_name(tj, self.taxa_names.inverse)
#        return adj_dict

    def visualize(self, level=0):
        dot = Digraph()
        for X in self.taxa_names.keys():
            dot.node(str(X), str(X))

        for x in range(len(self.adj_matrix)):
            ys = np.where(self.adj_matrix[x] <= level)[0]
            for y in ys:
                if y == x:
                    continue
                X = self.taxa_names.inverse[x]
                Y = self.taxa_names.inverse[y]
                dot.edge(X, Y)
        dot.render(view=True)
