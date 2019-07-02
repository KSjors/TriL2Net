from datastructures.rooted_level_k_network import RootedLevelKNetwork
import numpy as np
import itertools
from utils.help_functions import is_symmetry
from bidict import bidict
import logging
import copy


class RootedLevelKGenerator(RootedLevelKNetwork):
    def __init__(self, dir_adj_matrix: np.ndarray, necessary_edges: list, symmetrical_nodes: bidict, level: int = 2, dimension: int = 2):
        logging.debug("Creating generator network.")
        network = RootedLevelKNetwork.from_dir_adj_matrix(dir_adj_matrix=dir_adj_matrix, level=level, dimension=dimension)
        network.logger.debug("Created for generator network creation.")
        super().__init__(network.adj_matrix, network.node_name_map, network.leaf_names, network.level, network.dimension)
        self.logger = logging.getLogger('network.gen.{}'.format(self.uid))
        self.logger.debug("Created through network {}.".format(network.uid))
        self.necessary_edges = necessary_edges
        self.necessary_edges = necessary_edges
        self.symmetrical_nodes = symmetrical_nodes

    def build_trinets(self):
        """Build all possible trinets."""
        self.logger.debug("Building all possible trinets.")
        base_net = RootedLevelKNetwork.copy_network(self)
        # TODO check if extra leaves are added in standard way

        # # Add edges to generator which are necessary (in case of parallel arcs)
        # # Don't start add nodes from the end of these edges, as this will result in symmetries (before and after)
        # dont_end_in_nodes = []
        # necessary_leaf_dict = {}
        # for edge in self.necessary_edges:
        #     internal_name, leaf_name = base_net.add_leaf_to_edge(edge[0], edge[1])
        #     dont_end_in_nodes.append(internal_name)
        #     necessary_leaf_dict[leaf_name] = edge

        # Create iterator of possible combinations of leaves to add
        all_edges = base_net.get_edges(leafless=True)
        number_of_generator_leaves = base_net.number_of_leaves
        edges_iterator = itertools.combinations(all_edges, 3 - number_of_generator_leaves)

        # For each possible combination, create trinet and save it to trinets_gen_sides list
        trinet_info = {}
        for edges in edges_iterator:
            edges = list(edges)
            extra_leaf_dict = {}

            # # Check if any of the extra leaves will be added after a node in dont_start_from_nodes
            # from_nodes = [edge[1] for edge in edges]
            # if not set(from_nodes).isdisjoint(dont_end_in_nodes):
            #     continue

            # Add extra leaves to base net and save together with underlying generator (self) and added edges
            current_trinet = RootedLevelKNetwork.copy_network(base_net)
            for edge in edges:
                _, leaf_name = current_trinet.add_leaf_to_edge(edge[0], edge[1])
                extra_leaf_dict[leaf_name] = edge
            current_trinet.prune()
            if current_trinet.number_of_internals_leaves_reticulations()[2] != self.level:
                continue
            transformations = current_trinet.to_standard_form()
            # extra_leaf_dict.update(necessary_leaf_dict)

            # Ordering needs to be put through in extra_leaf_dict
            translation_dict = create_translation_dict(transformations)
            translated_extra_leaf_dict = {}
            for extra_leaf, on_edge in extra_leaf_dict.items():
                translated_extra_leaf = translate(translation_dict, extra_leaf)
                translated_on_edge = translate(translation_dict, on_edge)
                translated_extra_leaf_dict[translated_extra_leaf] = translated_on_edge

            on_edges = [tuple(on_edge) for leaf, on_edge in sorted(list(translated_extra_leaf_dict.items()))]

            trinet_info[current_trinet] = {'generator': self, 'on_edges': on_edges}

        # In case generator has only one leaf, also add two leaves to the same edge (only need to do one node of every symmetry pair)
        if number_of_generator_leaves == 1:
            extra_leaf_dict = {}
            edges_iterator = itertools.combinations(all_edges, 1)
            for edge in edges_iterator:
                current_trinet = RootedLevelKNetwork.copy_network(base_net)
                new_node_name, leaf_name_1 = current_trinet.add_leaf_to_edge(edge[0][0], edge[0][1])
                extra_leaf_dict[leaf_name_1] = edge[0]
                _, leaf_name_2 = current_trinet.add_leaf_to_edge(new_node_name, edge[0][1])
                extra_leaf_dict[leaf_name_2] = edge[0]
                # extra_leaf_dict.update(necessary_leaf_dict)
                current_trinet.prune()
                if current_trinet.number_of_internals_leaves_reticulations()[2] != self.level:
                    continue
                transformations = current_trinet.to_standard_form()

                # Ordering needs to be put through in extra_leaf_dict
                translation_dict = create_translation_dict(transformations)
                translated_extra_leaf_dict = {}
                for extra_leaf, on_edge in extra_leaf_dict.items():
                    translated_extra_leaf = translate(translation_dict, extra_leaf)
                    translated_on_edge = translate(translation_dict, on_edge)
                    translated_extra_leaf_dict[translated_extra_leaf] = translated_on_edge

                on_edges = [tuple(on_edge) for leaf, on_edge in sorted(list(translated_extra_leaf_dict.items()))]

                trinet_info[current_trinet] = {'generator': self, 'on_edges': on_edges}

        return trinet_info


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
