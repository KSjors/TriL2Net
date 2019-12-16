from datastructures.rooted_level_k_network import RootedLevelKNetwork, NetworkSet, NetworkInfo
import numpy as np
from utils.help_functions import guid
import itertools
from bidict import bidict
import logging
import copy


class RootedLevelKGenerator(RootedLevelKNetwork):
    def __init__(self, name, dir_adj_matrix: np.ndarray, symmetrical_nodes: bidict, level: int = 2, dimension: int = 2, check_valid: bool = True,
                 char_type='ALPH'):
        network = RootedLevelKNetwork.from_dir_adj_matrix(dir_adj_matrix=dir_adj_matrix, level=level, dimension=dimension, check_valid=check_valid,
                                                          char_type=char_type)
        super().__init__(network.adj_matrix, network.node_name_map, leaf_numbers=network.leaf_numbers, level=network.level, dimension=network.dimension)
        self.name = name
        self.symmetrical_nodes = symmetrical_nodes

    def build_trinets(self):
        base_net = copy.deepcopy(self)

        reticulations = self.leaves_below_nodes(set(self.reticulations), 1)

        # --------- Trinets -----------
        # Create iterator of possible combinations of leaves to add
        all_edges = base_net.internal_arcs
        number_of_generator_leaves = len(base_net.leaf_numbers)
        edges_iterator = itertools.combinations(all_edges, 3 - number_of_generator_leaves)

        # For each possible combination, create trinet and save it to trinets_gen_sides list
        trinet_info_list = NetworkSet(network_size=3)
        for edges in edges_iterator:
            extra_leaf_dict = {}
            # Add extra leaves to base net and save together with underlying generator (self) and added edges
            current_trinet = copy.deepcopy(base_net)
            for edge in edges:
                _, leaf_name = current_trinet.add_leaf_to_edge(edge, char_type='ALPH')
                extra_leaf_dict[leaf_name] = edge
            current_trinet.prune()
            if current_trinet.number_of_reticulations != self.level:
                continue
            current_trinet.reset_optimization_variables()
            current_trinet.calculate_optimization_variables()
            trinet_info = NetworkInfo(current_trinet, {'generator'        : self, 'generator_name': self.name, 'reticulations': reticulations,
                                                       'extra_leaf_dict'  : copy.deepcopy(extra_leaf_dict), 'strict_level': self.level,
                                                       'symmetrical_nodes': self.symmetrical_nodes})
            trinet_info_list.append(trinet_info)

        # In case generator has only one leaf, also add two leaves to the same edge (only need to do one node of every symmetry pair)
        if number_of_generator_leaves == 1:
            extra_leaf_dict = {}
            edges_iterator = itertools.combinations(all_edges, 1)
            for edge in edges_iterator:
                current_trinet = copy.deepcopy(base_net)
                current_trinet.reset_optimization_variables()
                new_node_name, leaf_name_1 = current_trinet.add_leaf_to_edge(edge[0], char_type='ALPH')
                extra_leaf_dict[leaf_name_1] = edge[0]
                _, leaf_name_2 = current_trinet.add_leaf_to_edge([new_node_name, edge[0][1]], char_type='ALPH')
                extra_leaf_dict[leaf_name_2] = edge[0]

                current_trinet.prune()
                if current_trinet.number_of_reticulations != self.level:
                    continue
                current_trinet.reset_optimization_variables()
                current_trinet.calculate_optimization_variables()
                trinet_info = NetworkInfo(current_trinet, {'generator'      : self, 'generator_name': self.name, 'reticulations': reticulations,
                                                           'extra_leaf_dict': copy.deepcopy(extra_leaf_dict),
                                                           'strict_level'   : self.level, 'symmetrical_nodes': self.symmetrical_nodes})
                trinet_info_list.append(trinet_info)
        return trinet_info_list

    def build_binets(self):
        base_net = copy.deepcopy(self)

        reticulations = self.leaves_below_nodes(set(self.reticulations), 1)

        # --------- Binets -----------
        # Create iterator of possible combinations of leaves to add
        all_edges = base_net.internal_arcs
        number_of_generator_leaves = len(base_net.leaf_numbers)
        edges_iterator = itertools.combinations(all_edges, 2 - number_of_generator_leaves)

        # For each possible combination, create binet and save it to trinets_gen_sides list
        binet_info_list = NetworkSet(network_size=3)
        for edges in edges_iterator:
            extra_leaf_dict = {}
            # Add extra leaves to base net and save together with underlying generator (self) and added edges
            current_binet = copy.deepcopy(base_net)
            current_binet.reset_optimization_variables()
            for edge in edges:
                _, leaf_name = current_binet.add_leaf_to_edge(edge, char_type='ALPH')
                extra_leaf_dict[leaf_name] = edge
            current_binet.prune()
            if current_binet.number_of_reticulations != self.level:
                continue
            current_binet.reset_optimization_variables()
            current_binet.calculate_optimization_variables()
            binet_info = NetworkInfo(current_binet, {'generator'        : self, 'generator_name': self.name, 'reticulations': reticulations,
                                                     'extra_leaf_dict'  : extra_leaf_dict, 'strict_level': self.level,
                                                     'symmetrical_nodes': self.symmetrical_nodes})
            binet_info_list.append(binet_info)
        return binet_info_list

    @property
    def sides(self):
        return self.edge_sides + self.reticulation_sides

    @property
    def edge_sides(self):
        return self.internal_arcs

    @property
    def reticulation_sides(self):
        reticulation_sides = []
        reticulations = self.reticulations
        for reticulation in reticulations:
            reticulation_child = self.nodes_below_nodes({reticulation}, max_depth=1).difference({reticulation}).pop()
            if self.is_leaf_node(reticulation_child):
                reticulation_sides.append(reticulation)
        return reticulation_sides

    @property
    def sets_of_symmetric_edge_sides(self):
        # TODO: make general
        symmetric_sides_sets_per_generator = {
            '1' : {
                ('0', '1'): [('0', '1'), ('0', '1')]
            },
            '2a': {
                ('0', '1'): [('0', '1')],
                ('0', '2'): [('0', '2')],
                ('1', '2'): [('1', '2')],
                ('1', '3'): [('1', '3')],
                ('2', '3'): [('2', '3')]
            },
            '2b': {
                ('0', '1'): [('0', '1')],
                ('0', '2'): [('0', '2')],
                ('1', '3'): [('1', '3')],
                ('1', '4'): [('1', '4')],
                ('3', '2'): [('3', '2')],
                ('3', '4'): [('3', '4')]
            },
            '2c': {
                ('0', '1'): [('0', '1'), ('0', '2')],
                ('1', '3'): [('1', '3'), ('2', '3')],
                ('1', '4'): [('1', '4'), ('2', '4')]
            },
            '2d': {
                ('0', '1'): [('0', '1')],
                ('0', '2'): [('0', '2')],
                ('1', '3'): [('1', '3'), ('1', '3')],
                ('3', '2'): [('3', '2')]
            }
        }
        return symmetric_sides_sets_per_generator[self.name]

    @property
    def sets_of_symmetric_reticulation_sides(self):
        return {key: [key] for key in self.reticulation_sides}

    @property
    def sets_of_symmetric_sides(self):
        return {**self.sets_of_symmetric_edge_sides, **self.sets_of_symmetric_reticulation_sides}

    @property
    def crucial_sets_of_symmetric_sides(self):
        result = []
        for key, sides in self.sets_of_symmetric_sides:
            if len(sides) > 1 and sides[0] == sides[1]:
                result.append(key)
        return result

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
        cp._o_cut_arc_matrix = copy.copy(self._o_cut_arc_matrix)
        cp._o_cut_arc_sets = copy.copy(self._o_cut_arc_sets)
        cp._o_biconnected_components = copy.copy(self._o_biconnected_components)
        cp._o_partial_ordering = copy.copy(self._o_partial_ordering)
        cp._o_leaf_ordering = copy.copy(self._o_leaf_ordering)
        return cp

    def __deepcopy__(self, memo):
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
        cp.symmetrical_nodes = copy.deepcopy(self.symmetrical_nodes)
        cp.logger = logging.getLogger('network.{}'.format(self.uid))
        cp._o_cut_arc_matrix = copy.deepcopy(self._o_cut_arc_matrix)
        cp._o_cut_arc_sets = copy.deepcopy(self._o_cut_arc_sets)
        cp._o_biconnected_components = copy.deepcopy(self._o_biconnected_components)
        cp._o_partial_ordering = copy.deepcopy(self._o_partial_ordering)
        cp._o_leaf_ordering = copy.deepcopy(self._o_leaf_ordering)
        return cp
