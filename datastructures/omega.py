import numpy as np
from tarjan import tarjan
from graphviz import Digraph
# from utils.help_functions import *
import logging


class Omega:
    def __init__(self, trinet_info_list):
        self.uid = guid()
        self.logger = logging.getLogger("omega.{}".format(self.uid))
        self.adj_matrix = np.zeros((trinet_set.number_of_taxa,  trinet_set.number_of_taxa))
        self.taxa_names = trinet_set.taxa_names
        cut_arc_sets_per_triplet = trinet_set.cut_arc_sets_per_triplet()

        for triplet, cut_arc_sets in cut_arc_sets_per_triplet.items():
            for cut_arc_set in cut_arc_sets:
                if len(cut_arc_set) == 2:
                    x = self.taxa_names[cut_arc_set[0]]
                    y = self.taxa_names[cut_arc_set[1]]
                    z = self.taxa_names[[Z for Z in triplet if Z not in cut_arc_set][0]]
                    self.adj_matrix[x, z] += 1
                    self.adj_matrix[y, z] += 1
        self.adj_matrix += np.eye(len(self.taxa_names))*len(self.taxa_names)

    def minimal_sink_sets(self, level=0):
        self.logger.debug("Computing minimal sink sets.")
        adj_dict = {}
        for x in range(len(self.adj_matrix)):
            ys = set(np.where(self.adj_matrix[x] <= level)[0])
            ys -= {x}
            adj_dict[x] = ys

        tj = tarjan(adj_dict)
        return list_of_list_to_name(tj, self.taxa_names.inverse)

    def visualize(self, level=0):
        self.logger.debug("Visualizing level {}.".format(level))
        dot = Digraph()
        for node_name in self.taxa_names.keys():
            dot.node(str(node_name), str(node_name))

        for x in range(len(self.adj_matrix)):
            ys = np.where(self.adj_matrix[x] <= level)[0]
            for y in ys:
                if y == x:
                    continue
                node_name_from = self.taxa_names.inverse[x]
                node_name_to = self.taxa_names.inverse[y]
                dot.edge(node_name_from, node_name_to)
        dot.render(view=True)
