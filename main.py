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


regenerate_trinets()
all_generators, all_trinets, all_trinets_gen_sides = get_trinets()

all_generators[-1].visualize()
# for trinet in all_trinets:
#     trinet.visualize()
#     time.sleep(1)

# graph = test_networks.B
# network = RootedLevelKNetwork.from_connections_dict(graph)
#
# trinet = RootedLevelKNetwork.trinet_from_network(network, {'H', 'F', 'B'})
# print(trinet.to_df(directed=False))
# print(trinet._get_cut_arc_matrix())
# print(trinet.get_biconnected_components())
# trinet.visualize()
#
