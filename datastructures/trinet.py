from datastructures.rooted_level_k_network import RootedLevelKNetwork
import time


class Trinet(RootedLevelKNetwork):
    def __init__(self, network, triplet, level=2, dimension=2):
        super().__init__(adj_matrix=network.adj_matrix, node_names=network.node_names, leaf_names=triplet, level=level, dimension=dimension)
        self.prune(suppress_redundant='none', suppress_parallel=False)
        self.to_standard_form()


