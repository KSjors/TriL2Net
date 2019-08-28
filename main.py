import logging

# TODO
#  Isomorphism checking
#  Random level-k network generation

logging_level = logging.INFO
# set up logging to file - see previous section for more details
logging.basicConfig(level=logging_level,
                    format='%(asctime)s %(name)-12s %(levelname)-8s %(message)s',
                    datefmt='%m-%d %H:%M',
                    filename='main.log',
                    filemode='w')
# define a Handler which writes INFO messages or higher to the sys.stderr
console = logging.StreamHandler()
console.setLevel(logging.DEBUG)

fh = logging.FileHandler(filename='file_log.log')
fh.setLevel(logging_level)
# set a format which is simpler for console use
formatter = logging.Formatter('%(name)-12s: %(levelname)-8s %(message)s')
# tell the handler to use this format
console.setFormatter(formatter)
fh.setFormatter(formatter)
# add the handler to the root logger
logging.getLogger('').addHandler(console)
logging.getLogger('network').addHandler(fh)

import os
import data.test_networks as test_networks
from datastructures.rooted_level_k_network import RootedLevelKNetwork
from solver import Solver
import data.generators as generator
from data.all_trinets import get_trinets, regenerate_trinets, pickle_save, pickle_read
from utils.help_functions import *
import timeit


rebuild = {
    'generators': 0
    , 'network' : 0
    , 'trinets' : 0
}

os.environ["PATH"] += os.pathsep + 'C:/Program Files (x86)/Graphviz2.38/bin/'
filename_save_trinets_B = 'data/trinets'

""" Build or load trinets """
if rebuild['generators']:
    regenerate_trinets()
_, standard_trinet_info_list = get_trinets()

""" Build or load network """
if rebuild['network']:
    dct = test_networks.connections_0
    main_network = RootedLevelKNetwork.from_connections_dict(dct)
    dct_1 = test_networks.connections_1
    network_1 = RootedLevelKNetwork.from_connections_dict(dct_1)
    dct_2a = test_networks.connections_2a
    network_2a = RootedLevelKNetwork.from_connections_dict(dct_2a)
    dct_2b = test_networks.connections_2b
    network_2b = RootedLevelKNetwork.from_connections_dict(dct_2b)
    dct_2c = test_networks.connections_2c
    network_2c = RootedLevelKNetwork.from_connections_dict(dct_2c)
    dct_2d = test_networks.connections_2d
    network_2d = RootedLevelKNetwork.from_connections_dict(dct_2d)

    main_network.replace_leaf_with_network('D1', network_1)
    main_network.standardize_internal_node_names()
    # main_network.replace_leaf_with_network('D2', network_2a)
    # main_network.standardize_internal_node_names()
    # main_network.replace_leaf_with_network('D3', network_2b)
    # main_network.standardize_internal_node_names()
    # main_network.replace_leaf_with_network('D4', network_2c)
    # main_network.standardize_internal_node_names()
    # main_network.replace_leaf_with_network('D5', network_2d)
    # main_network.standardize_internal_node_names()

    pickle_save("data/network.pickle", main_network)
main_network = pickle_read("data/network.pickle")
main_network.visualize()

""" Build or load trinets """
if rebuild['trinets'] or rebuild['network']:
    trinet_info_list = main_network.get_exhibited_trinets()
    pickle_save('data/trinets', trinet_info_list)
trinet_info_list = pickle_read('data/trinets')

# """ Solver """
# solver = Solver(trinet_info_list)
#
# solver.solve()


a = '((1, ((2, (3, (4)Y#H1)g)e, (((Y#H1, 5)h, 6)f)X#H2)c)a, ((X#H2, 7)d, 8)b)r'
network = RootedLevelKNetwork.from_enewick(a)
network.visualize()

""" Play around """
# network_solved = solver.solve()
# print(network.stable_ancestors(['C', 'D']))

# trinet.visualize()
# print(network.equal_structure(network_solved))


