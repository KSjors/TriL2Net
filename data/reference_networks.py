import data.generators as generators
import logging
import pickle
import itertools
import time
from tqdm import tqdm
import copy
from datastructures.rooted_level_k_network import NetworkSet, NetworkInfo, RootedLevelKNetwork


def pickle_save(filename, data):
    pickle_out = open(filename, 'wb')
    pickle.dump(data, pickle_out)
    pickle_out.close()


def pickle_read(filename):
    pickle_in = open(filename, 'rb')
    result = pickle.load(pickle_in)
    pickle_in.close()
    return result


def regenerate_standard_binets_trinets() -> None:
    """Regenerate and save all possible trinets."""
    logging.debug("Regenerating all possible trinets and saving them.")

    # All generators per level
    all_generators = {
        0  : [generators.generator_level0]
        , 1: [generators.generator_level1]
        , 2: [generators.generator_A, generators.generator_B, generators.generator_C, generators.generator_D]
    }

    # Get biconnected binets and trinets for each generator
    biconnected_trinet_list = NetworkSet(equal_naming=False)
    biconnected_binet_list = NetworkSet(equal_naming=False)
    for level, generator_list in all_generators.items():
        for generator in generator_list:
            generator_trinet_info_list = generator.build_trinets()
            generator_binet_info_list = generator.build_binets()
            biconnected_trinet_list.extend(generator_trinet_info_list)
            biconnected_binet_list.extend(generator_binet_info_list)

    biconnected_trinet_list.set_multiplicities_to_one()
    biconnected_binet_list.set_multiplicities_to_one()

    # From binets create trinets with two biconnected components
    two_component_trinet_list = NetworkSet()
    two_binet_infos_iterator = itertools.product(biconnected_binet_list.per_network_info(), repeat=2)
    for binet_infos in two_binet_infos_iterator:
        for index, leaf_name in enumerate(binet_infos[0].network.leaf_names):
            two_component_trinet = copy.deepcopy(binet_infos[0].network)
            two_component_trinet.replace_leaf_with_network(leaf_name, binet_infos[1].network, replace_names=True, char_type='ALPH')
            two_component_trinet_list.append(NetworkInfo(two_component_trinet))

    two_component_trinet_list.extend(biconnected_trinet_list)
    biconnected_trinet_list.extend(biconnected_binet_list)
    for network_info in biconnected_binet_list.per_network_info():
        network_info.network.reset_optimization_variables()
        network_info.network.calculate_optimization_variables()
    biconnected_trinet_list.calculate_info()
    two_component_trinet_list.set_multiplicities_to_one()
    for network_info in two_component_trinet_list.per_network_info():
        network_info.network.reset_optimization_variables()
        network_info.network.calculate_optimization_variables()
    two_component_trinet_list.calculate_info()

    pickle_out = open("data/all_networks_save.pickle", 'wb+')
    data = [all_generators, biconnected_trinet_list, two_component_trinet_list]
    pickle.dump(data, pickle_out)
    pickle_out.close()


def get_standard_binets_trinets() -> (list, NetworkSet):
    """Read and retrieve all possible trinets."""
    logging.debug("Reading and retrieving all possible trinets.")
    pickle_in = open("data/all_networks_save.pickle", 'rb')
    result = pickle.load(pickle_in)
    pickle_in.close()
    all_generators, biconnected_trinet_binet_list, all_trinet_list = result[0], result[1], result[2]
    return all_generators, biconnected_trinet_binet_list, all_trinet_list


LEVEL_1_2_GENERATORS, BICONNECTED_BINET_TRINET_LIST, TRINET_LIST = get_standard_binets_trinets()

# for ti in TRINET_LIST.per_network_info():
#     n = ti.network
#     if n.number_of_reticulations == 1 and n.biconnected:
#         ti.network.visualize()
#
# kk
