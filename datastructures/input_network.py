import logging
from datastructures.rooted_level_k_network import RootedLevelKNetwork
from datastructures.trinet import Trinet
from bidict import bidict
import itertools
import numpy as np
from tqdm import tqdm



class InputRootedNetwork(RootedLevelKNetwork):
    def __init__(self, dir_adj_matrix: np.ndarray = None, adj_matrix: np.ndarray = None, node_names: bidict = None, leaf_names: set = None, level: int = 2,
                 dimension: int = 2):
        super().__init__(dir_adj_matrix, adj_matrix, node_names, leaf_names, level, dimension)
        self.is_valid()



    def is_valid(self) -> bool:
        """Check if network is valid: has right degrees and number of roots and reticulations"""
        logging.debug("Checking validity.")
        in_degrees = self.get_in_degrees()
        out_degrees = self.get_out_degrees()

        number_of_roots = sum(in_degrees == 0)
        assert number_of_roots == 1, "Network is not valid, it has more than one root ({}).".format(number_of_roots)
        assert sum(in_degrees <= self.dimension-1), "Network is not valid, it has a node with more arcs entering than allowed."
        assert sum(out_degrees <= self.dimension-1), "Network is not valid, it has a node with more arcs exiting than allowed."
        assert sum((in_degrees + out_degrees) <= self.dimension), "Network is not valid, it has a node with more arcs entering and exiting than allowed."
        # TODO: check if leaves in the end of the matrix

        biconnected_components = self.get_biconnected_components()
        for bc in biconnected_components:
            node_numbers_bc = self._get_node_numbers(bc)
            in_sum_bc = in_degrees[node_numbers_bc]
            number_of_reticulations = sum(in_sum_bc > 1)
            assert number_of_reticulations <= self.level, "Biconnected component {} has to many reticulations.".format(bc)
        return True
