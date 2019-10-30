import logging

# TODO
#  -- CORRECT INPUT --
#  LEVEL-1 input
#  -- TEST ALTERED INPUT --
#  Swap leaves in some trinets
#  Swap trinet for other trinet
#  -- PROOF --
#  Pseudocode interesting functions
#  Proof that for correct input output also correct
#  Proof auxilliary graph theorem
#  -- OPTIMIZE --
#  Pruning
#  Isomorphism checking --> separate code for trees?
#  -- RANDOM --
#  Do not add leaf which increases level
#  -- WRITE --


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
from datastructures.rooted_level_k_network import RootedLevelKNetwork, NetworkInfoList
from solver import Solver
import data.generators as generator
from data.all_trinets import get_standard_networks, regenerate_standard_networks, pickle_save, pickle_read
from utils.help_functions import *
import timeit
from tqdm import tqdm
import matplotlib.pyplot as plt
import scipy.stats as stats
from mip import model
from data.ETS_network import ETS_NETWORK_dict

if __name__ == '__main__':
    config = {
        'rebuild_generators': 0
        , 'rebuild_network' : 1
        , 'network_config'  : {
            'build_type'     : 'evolve'
            # dict, enewick, evolve
            , 'base_net_dict': test_networks.ETS_NETWORK
            , 'enewick'      : '(a, b)0'
            , 'evolve_config': {
                'level'               : 2
                , 'times'             : 20
                , 'reticulate_chance' : 0.1
                , 'termination_chance': 0
            }
        }
        , 'rebuild_trinets' : 0
        , 'trinet_config'   : {
            'replace_network_distort' : 0
            , 'switch_leaves_distort' : 0
            , 'move_leaf_distort'     : 0
            , 'remove_network_distort': 0
        }
    }

    os.environ["PATH"] += os.pathsep + 'C:/Program Files (x86)/Graphviz2.38/bin/'

    """ Build or load trinets """
    if config['rebuild_generators']:
        regenerate_standard_networks()
    _, biconnected_trinet_binet_list, all_trinet_list = get_standard_networks()

    """ Build or load network """
    if config['rebuild_network']:
        network_config = config['network_config']
        if network_config['build_type'] == 'dict':
            dct = network_config['base_net_dict']
            network = RootedLevelKNetwork.from_connections_dict(dct)
            network.standardize_internal_node_names()
        elif network_config['build_type'] == 'enewick':
            enewick = network_config['enewick']
            network = RootedLevelKNetwork.from_enewick(enewick)
        elif network_config['build_type'] == 'evolve':
            evolve_config = network_config['evolve_config']
            enewick = '(a, b)0'
            network = RootedLevelKNetwork.from_enewick(enewick)
            network.level = evolve_config['level']
            network.evolve_times(evolve_config['times'], evolve_config['reticulate_chance'])
            network.terminate(evolve_config['termination_chance'])
            network.standardize_node_names()

        pickle_save("data/network.pickle", network)
    network = pickle_read('data/network.pickle')
    # network.visualize()

    """ Build or load trinets """
    # if config['rebuild_network']:
    #     trinet_info_list = network.network_set(max_processes=6, progress_bar=True)
    #     pickle_save('data/trinets', trinet_info_list)
    # trinet_info_list = pickle_read('data/trinets')

    """ Distort trinets """
    # if config['rebuild_trinets'] or config['rebuild_network']:
    #     distorted_trinet_info_list = copy.deepcopy(trinet_info_list)
    #     trinet_config = config['trinet_config']
    #     distorted = False
    #     if trinet_config['replace_network_distort'] != 0:
    #         distorted_trinet_info_list.replace_network_distort(trinet_config['replace_network_distort'], all_trinet_list)
    #         distorted = True
    #     if trinet_config['switch_leaves_distort'] != 0:
    #         distorted_trinet_info_list.switch_leaves_distort(trinet_config['switch_leaves_distort'])
    #         distorted = True
    #     if trinet_config['move_leaf_distort'] != 0:
    #         distorted_trinet_info_list.move_leaf_distort(trinet_config['move_leaf_distort'])
    #         distorted = True
    #     if trinet_config['remove_network_distort'] != 0:
    #         distorted_trinet_info_list.remove_network_distort(trinet_config['remove_network_distort'])
    #         distorted = True
    #     pickle_save('data/distorted_trinets', distorted_trinet_info_list)
    #     if distorted:
    #         print(distorted_trinet_info_list.summary())
    # distorted_trinet_info_list = pickle_read('data/distorted_trinets')

    # lst = network.get_exhibited_clusters(progress_bar=True)
    network.visualize()
    T1 = 0
    T2 = 0
    k = 5
    for _ in range(k):
        t1 = time.time()
        lst1 = NetworkInfoList.exhibited_networks_of_network(network, 3, progress_bar=True, max_processes=4)
        t2 = time.time()
        lst2 = NetworkInfoList.exhibited_networks_of_network_2(network, 3, progress_bar=True, max_processes=4)
        t3 = time.time()
        T1 += t2 - t1
        T2 += t3 - t2

    print(len(lst1), len(lst2))
    print(T1/k, T2/k)

    # ETS_NETWORK = RootedLevelKNetwork.from_connections_dict(ETS_NETWORK_dict)
    # print(ETS_NETWORK.full_network_set_consistency(ETS_NETWORK, max_processes=6, progress_bar=True))
    # # ETS_NETWORK.standardize_node_names()
    # # ETS_NETWORK = RootedLevelKNetwork.get_network_below_node(ETS_NETWORK, '12')
    # # ETS_NETWORK.standardize_node_names()
    # # ETS_NETWORK.visualize(file_path='example_network', internal_node_labels=True, edge_labels=Falsex    , rankdir='LR', format='pdf')
