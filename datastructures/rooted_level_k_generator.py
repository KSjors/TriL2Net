from datastructures.rooted_level_k_network import RootedLevelKNetwork, TrinetInfoList, TrinetInfo
import numpy as np
from utils.help_functions import guid
import itertools
from bidict import bidict
import logging
import copy


class RootedLevelKGenerator(RootedLevelKNetwork):
    def __init__(self, name, dir_adj_matrix: np.ndarray, symmetrical_nodes: bidict, level: int = 2, dimension: int = 2, check_valid: bool = True):
        logging.debug("Creating generator network.")
        network = RootedLevelKNetwork.from_dir_adj_matrix(dir_adj_matrix=dir_adj_matrix, level=level, dimension=dimension, check_valid=check_valid)
        network.logger.debug("Created for generator network creation.")
        super().__init__(network.adj_matrix, network.node_name_map, leaf_numbers=network.leaf_numbers, level=network.level, dimension=network.dimension)
        self.name = name
        self.logger = logging.getLogger('network.gen.{}'.format(self.uid))
        self.logger.debug("Created through network {}.".format(network.uid))
        self.symmetrical_nodes = symmetrical_nodes

    def build_trinets(self):
        """Build all possible trinets."""
        self.logger.debug("Building all possible trinets.")
        base_net = copy.deepcopy(self)

        reticulations = self.get_leaf_children(set(self.get_reticulations()), 1)

        # --------- Trinets -----------
        # Create iterator of possible combinations of leaves to add
        all_edges = base_net.get_edges(leafless=True)  # level 0 should not use leafless
        number_of_generator_leaves = len(base_net.leaf_numbers)
        edges_iterator = itertools.combinations(all_edges, 3 - number_of_generator_leaves)

        # For each possible combination, create trinet and save it to trinets_gen_sides list
        trinet_info_list = TrinetInfoList()
        for edges in edges_iterator:
            extra_leaf_dict = {}
            # Add extra leaves to base net and save together with underlying generator (self) and added edges
            current_trinet = copy.deepcopy(base_net)
            for edge in edges:
                _, leaf_name = current_trinet.add_leaf_to_edge(edge)
                extra_leaf_dict[leaf_name] = edge
            current_trinet.prune()
            if current_trinet.number_of_internals_leaves_reticulations()[2] != self.level:
                continue
            trinet_info = TrinetInfo(current_trinet, {'generator'        : self, 'generator_name': self.name, 'reticulations': reticulations,
                                                      'extra_leaf_dict'  : copy.deepcopy(extra_leaf_dict), 'level': self.level,
                                                      'symmetrical_nodes': self.symmetrical_nodes})
            trinet_info_list.append(trinet_info)

        # In case generator has only one leaf, also add two leaves to the same edge (only need to do one node of every symmetry pair)
        if number_of_generator_leaves == 1:
            extra_leaf_dict = {}
            edges_iterator = itertools.combinations(all_edges, 1)
            for edge in edges_iterator:
                current_trinet = copy.deepcopy(base_net)
                new_node_name, leaf_name_1 = current_trinet.add_leaf_to_edge(edge[0])
                extra_leaf_dict[leaf_name_1] = edge[0]
                _, leaf_name_2 = current_trinet.add_leaf_to_edge([new_node_name, edge[0][1]])
                extra_leaf_dict[leaf_name_2] = edge[0]

                current_trinet.prune()
                if current_trinet.number_of_internals_leaves_reticulations()[2] != self.level:
                    continue
                trinet_info = TrinetInfo(current_trinet, {'generator'      : self, 'generator_name': self.name, 'reticulations': reticulations,
                                                          'extra_leaf_dict': copy.deepcopy(extra_leaf_dict),
                                                          'level'          : self.level, 'symmetrical_nodes': self.symmetrical_nodes})
                trinet_info_list.append(trinet_info)

        # --------- Binets -----------

        # Create iterator of possible combinations of leaves to add
        all_edges = base_net.get_edges(leafless=True)
        number_of_generator_leaves = len(base_net.leaf_numbers)
        edges_iterator = itertools.combinations(all_edges, 2 - number_of_generator_leaves)

        # For each possible combination, create binet and save it to trinets_gen_sides list
        for edges in edges_iterator:
            extra_leaf_dict = {}
            # Add extra leaves to base net and save together with underlying generator (self) and added edges
            current_trinet = copy.deepcopy(base_net)
            for edge in edges:
                _, leaf_name = current_trinet.add_leaf_to_edge(edge)
                extra_leaf_dict[leaf_name] = edge
            current_trinet.prune()
            if current_trinet.number_of_internals_leaves_reticulations()[2] != self.level:
                continue
            trinet_info = TrinetInfo(current_trinet, {'generator'      : self, 'generator_name': self.name, 'reticulations': reticulations,
                                                      'extra_leaf_dict': extra_leaf_dict, 'level': self.level, 'symmetrical_nodes': self.symmetrical_nodes})
            trinet_info_list.append(trinet_info)
        return trinet_info_list

    def __copy__(self):
        cls = self.__class__
        cp = cls.__new__(cls)
        cp.adj_matrix = copy.copy(self.adj_matrix)
        cp.node_name_map = copy.copy(self.node_name_map)
        cp.leaf_numbers = copy.copy(self.leaf_numbers)
        cp.number_of_nodes = copy.copy(self.number_of_nodes)
        cp.level = copy.copy(self.level)
        cp.dimension = copy.copy(self.dimension)
        cp.uid = guid()
        cp.name = copy.copy(self.name)
        cp.symmetrical_nodes = copy.copy(self.symmetrical_nodes)
        cp.logger = logging.getLogger('network.{}'.format(self.uid))
        cp._biconnected_components = copy.copy(self._biconnected_components)
        cp._partial_ordering = copy.copy(self._partial_ordering)
        cp._cut_arc_matrix = copy.copy(self._cut_arc_matrix)
        cp._cut_arc_sets = copy.copy(self._cut_arc_sets)
        return cp

    def __deepcopy__(self, memo):
        cls = self.__class__
        cp = cls.__new__(cls)
        cp.adj_matrix = copy.deepcopy(self.adj_matrix)
        cp.node_name_map = copy.deepcopy(self.node_name_map)
        cp.leaf_numbers = copy.deepcopy(self.leaf_numbers)
        cp.number_of_nodes = copy.deepcopy(self.number_of_nodes)
        cp.level = copy.deepcopy(self.level)
        cp.dimension = copy.deepcopy(self.dimension)
        cp.uid = guid()
        cp.name = copy.deepcopy(self.name)
        cp.symmetrical_nodes = copy.deepcopy(self.symmetrical_nodes)
        cp.logger = logging.getLogger('network.{}'.format(self.uid))
        cp._biconnected_components = copy.deepcopy(self._biconnected_components)
        cp._partial_ordering = copy.deepcopy(self._partial_ordering)
        cp._cut_arc_matrix = copy.deepcopy(self._cut_arc_matrix)
        cp._cut_arc_sets = copy.deepcopy(self._cut_arc_sets)
        return cp


def create_translation_dict(transformations):
    translation_dict = {}
    for transform in transformations:
        if len(transform) == 2:
            translation_dict[transform[0]] = transform[1]
            translation_dict[transform[1]] = transform[0]
    return translation_dict


def translate(translation_dict, value):
    if type(value) == list:
        result = []
        for v in value:
            result.append(translate(translation_dict, v))
        return result
    else:
        if value in translation_dict:
            return translation_dict[value]
        else:
            return value
