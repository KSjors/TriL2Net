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
from datastructures.rooted_level_k_network import RootedLevelKNetwork, NetworkSet, Omega
from data.reference_networks import get_standard_binets_trinets, regenerate_standard_binets_trinets, pickle_save, pickle_read
from utils.help_functions import *
from data.ETS_network import ETS_NETWORK_dict
from data.generators import *
from config import settings

if __name__ == '__main__':
    config = {
        'rebuild_generators': 0
        , 'rebuild_network' : 0
        , 'network_config'  : {
            'build_type'     : 'evolve'  # dict, enewick, evolve
            , 'base_net_dict': ETS_NETWORK_dict
            , 'enewick'      : '(a, b)0'
            , 'evolve_config': {
                'level'               : 2
                , 'times'             : 20
                , 'reticulate_chance' : 0.3
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
        regenerate_standard_binets_trinets()
    _, biconnected_trinet_binet_list, all_trinet_list = get_standard_binets_trinets()

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
            network.evolve_times(evolve_config['times'], evolve_config['reticulate_chance'], progress_bar=True)
            network.terminate_percentage_leaves(evolve_config['termination_chance'])
            network.standardize_node_names()

        pickle_save("data/network.pickle", network)
    network = pickle_read('data/network.pickle')
    if config['rebuild_network']:
        network.visualize()

    if config['rebuild_trinets'] or config['rebuild_network']:
        trinet_info_list = NetworkSet.induced_strict_network_set(network, 3, 1, True, method='Iterative')
        pickle_save('data/trinets', trinet_info_list)
    # trinet_info_list = pickle_read('data/trinets')

    """ Distort trinets """
    if config['rebuild_trinets'] or config['rebuild_network']:
        distorted_trinet_info_list = copy.deepcopy(trinet_info_list)
        trinet_config = config['trinet_config']
        distorted = False
        if trinet_config['replace_network_distort'] != 0:
            distorted_trinet_info_list.uniform_distort(trinet_config['replace_network_distort'], all_trinet_list)
            distorted = True
        if trinet_config['switch_leaves_distort'] != 0:
            distorted_trinet_info_list.switch_leaves_distort(trinet_config['switch_leaves_distort'])
            distorted = True
        if trinet_config['move_leaf_distort'] != 0:
            distorted_trinet_info_list.move_leaf_distort(trinet_config['move_leaf_distort'])
            distorted = True
        if trinet_config['remove_network_distort'] != 0:
            distorted_trinet_info_list.deletion_distort(trinet_config['remove_network_distort'])
            distorted = True
        pickle_save('data/distorted_trinets', distorted_trinet_info_list)
        if distorted:
            print(distorted_trinet_info_list.summary())
    # distorted_trinet_info_list = pickle_read('data/distorted_trinets')



    # enewick = '((((c, d)3, b)2, a)1, (((e, f)6, (((h, i)9, (j)10)8, g)7)5, ((10, (k, 14)13)12, ((l)14, m)15)11)4)0'
    # b = RootedLevelKNetwork.from_enewick(enewick)
    # T = NetworkSet.displayed_trees(b)
    # ts = []
    # for t in T.per_network_info():
    #     ts.append(t.network)
    #
    # T1 = RootedLevelKNetwork.restrict(ts[0], list('afhijk'))
    # T2 = RootedLevelKNetwork.restrict(ts[-1], list('cehijl'))
    #
    # T1.visualize(internal_node_labels=False, rankdir='LR', file_path=None, format='pdf')
    # T2.visualize(internal_node_labels=False, rankdir='LR', file_path=None, format='pdf')
    #
    # T1.visualize(internal_node_labels=False, rankdir='LR', file_path='exhibited_cluster_1', format='pdf')
    # T2.visualize(internal_node_labels=False, rankdir='LR', file_path='exhibited_cluster_2', format='pdf')
