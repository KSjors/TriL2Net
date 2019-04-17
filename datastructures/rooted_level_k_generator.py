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
        super().__init__(network.adj_matrix, network.node_names, network.leaf_names, network.level, network.dimension)
        self.logger = logging.getLogger('network.gen.{}'.format(self.uid))
        self.logger.debug("Created through network {}.".format(network.uid))
        self.necessary_edges = necessary_edges
        self.symmetrical_nodes = symmetrical_nodes

    def build_trinets(self):
        """Build all possible trinets."""
        self.logger.debug("Building all possible trinets.")
        base_net = RootedLevelKNetwork.copy_network(self)

        # Add edges to generator which are necessary (in case of parallel arcs)
        # Don't start add nodes from the end of these edges, as this will result in symmetries (before and after)
        dont_end_in_nodes = []
        for edge in self.necessary_edges:
            internal_name, leaf_name = base_net.add_leaf_to_edge(edge[0], edge[1])
            dont_end_in_nodes.append(internal_name)

        # Create iterator of possible combinations of leaves to add
        edges = base_net.get_edges(leafless=True)
        number_of_generator_leaves = base_net.number_of_leaves
        extra_leaves_iterator = itertools.combinations(edges, 3 - number_of_generator_leaves)

        # For each possible combination, create trinet and save it to trinets_gen_sides list
        symmetry_check_list = []
        trinets_gen_sides = []
        for extra_leaves in extra_leaves_iterator:
            extra_leaves = list(extra_leaves)
            extra_leaves.sort()

            # Check if any of the extra leaves will be added after a node in dont_start_from_nodes
            from_nodes = [extra_leaf[1] for extra_leaf in extra_leaves]
            if not set(from_nodes).isdisjoint(dont_end_in_nodes):
                continue

            # # Check if symmetry of extra_leaves has already been done
            # if is_symmetry(self.symmetrical_nodes, symmetry_check_list, extra_leaves):
            #     continue
            # else:
            #     symmetry_check_list.append(extra_leaves)

            # Add extra leaves to base net and save together with underlying generator (self) and added edges
            current_trinet = RootedLevelKNetwork.copy_network(base_net)
            for extra_leaf in extra_leaves:
                current_trinet.add_leaf_to_edge(extra_leaf[0], extra_leaf[1])
            current_trinet.to_standard_form()
            added_edges = self.necessary_edges + extra_leaves
            added_edges.sort()
            trinets_gen_sides.append([current_trinet, self, added_edges])

        # In case generator has only one leaf, also add two leaves to the same edge (only need to do one node of every symmetry pair)
        symmetry_check_list = []
        if number_of_generator_leaves == 1:
            extra_leaves_iterator = itertools.combinations(edges, 1)
            for extra_leaves in extra_leaves_iterator:
                # Check if symmetry of extra_leaves has already been done
                if is_symmetry(self.symmetrical_nodes, symmetry_check_list, extra_leaves):
                    continue
                else:
                    symmetry_check_list.append(extra_leaves)
                current_trinet = RootedLevelKNetwork.copy_network(base_net)
                new_node_name, _ = current_trinet.add_leaf_to_edge(extra_leaves[0][0], extra_leaves[0][1])
                current_trinet.add_leaf_to_edge(new_node_name, extra_leaves[0][1])
                current_trinet.to_standard_form()
                trinets_gen_sides.append([current_trinet, self, (extra_leaves[0], copy.deepcopy(extra_leaves[0]))])

        return trinets_gen_sides
