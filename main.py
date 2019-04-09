import data.test_networks as test_networks
from datastructures.rooted_level_k_network import *
from datastructures.trinet_set import *
import data.generators as generator
from data.all_trinets import *
from utils.help_functions import *
from solver import *
import time

import logging

logging.basicConfig(level=logging.DEBUG)

os.environ["PATH"] += os.pathsep + 'C:/Program Files (x86)/Graphviz2.38/bin/'


# regenerate_trinets()
all_generators, all_trinets, all_trinets_gen_sides = get_trinets()


# g = generators.D
# n = RootedLevelKNetwork.from_dir_adj_matrix(g)
# n.visualize()
# print(n.to_df())
# print(n._to_block_form())
# print(n.to_standard_form_gen_2())
# print(n.to_df())
# # graph = test_networks.B
# network = RootedLevelKNetwork.from_connections_dict(graph)
#
# network_trinets = network.get_exhibited_trinets()

