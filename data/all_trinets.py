import data.generators as generators
import logging
import pickle
import itertools
import time
from tqdm import tqdm
import copy
from datastructures.rooted_level_k_network import NetworkInfoList, NetworkInfo, RootedLevelKNetwork


def pickle_save(filename, data):
    pickle_out = open(filename, 'wb')
    pickle.dump(data, pickle_out)
    pickle_out.close()


def pickle_read(filename):
    pickle_in = open(filename, 'rb')
    result = pickle.load(pickle_in)
    pickle_in.close()
    return result


def regenerate_standard_networks() -> None:
    """Regenerate and save all possible trinets."""
    logging.debug("Regenerating all possible trinets and saving them.")

    # All generators per level
    all_generators = {
        0: [generators.generator_level0]
        , 1: [generators.generator_level1]
        , 2: [generators.generator_A, generators.generator_B, generators.generator_C, generators.generator_D]
    }

    # Get biconnected binets and trinets for each generator
    biconnected_trinet_list = NetworkInfoList(network_size=3)
    biconnected_binet_list = NetworkInfoList(network_size=3)
    for level, generator_list in all_generators.items():
        for generator in generator_list:
            generator_trinet_info_list = generator.build_trinets()
            generator_binet_info_list = generator.build_binets()
            biconnected_trinet_list.extend(generator_trinet_info_list)
            biconnected_binet_list.extend(generator_binet_info_list)

    biconnected_trinet_list.uniquify(count=False)
    biconnected_binet_list.uniquify(count=False)

    # From binets create trinets with two biconnected components
    two_component_trinet_list = NetworkInfoList(network_size=3)
    two_binet_infos_iterator = itertools.combinations(biconnected_binet_list, 2)
    for binet_infos in two_binet_infos_iterator:
        previous_two_component_trinet = None
        for index, leaf_name in enumerate(binet_infos[0].network.leaf_names):
            two_component_trinet = copy.deepcopy(binet_infos[0].network)
            two_component_trinet.replace_leaf_with_network(leaf_name, binet_infos[1].network, replace_names=True, char_type='ALPH')
            if index == 0 or not two_component_trinet.equal_structure(previous_two_component_trinet):
                two_component_trinet_list.append(NetworkInfo(two_component_trinet))
            previous_two_component_trinet = two_component_trinet

    two_component_trinet_list.extend(biconnected_trinet_list)
    biconnected_trinet_list.extend(biconnected_binet_list)
    two_component_trinet_list.uniquify(count=False)
    biconnected_trinet_list.calculate_info()
    two_component_trinet_list.calculate_info()

    pickle_out = open("data/all_networks_save.pickle", 'wb')
    data = [all_generators, biconnected_trinet_list, two_component_trinet_list]
    pickle.dump(data, pickle_out)
    pickle_out.close()


def get_standard_networks() -> (list, NetworkInfoList):
    """Read and retrieve all possible trinets."""
    logging.debug("Reading and retrieving all possible trinets.")
    pickle_in = open("data/all_networks_save.pickle", 'rb')
    result = pickle.load(pickle_in)
    pickle_in.close()
    all_generators, biconnected_trinet_binet_list, all_trinet_list = result[0], result[1], result[2]
    return all_generators, biconnected_trinet_binet_list, all_trinet_list
