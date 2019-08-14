import data.generators as generators
import logging
import pickle
from datastructures.rooted_level_k_network import TrinetInfoList, TrinetInfo


def pickle_save(filename, data):
    pickle_out = open(filename, 'wb')
    pickle.dump(data, pickle_out)
    pickle_out.close()


def pickle_read(filename):
    pickle_in = open(filename, 'rb')
    result = pickle.load(pickle_in)
    pickle_in.close()
    return result


def regenerate_trinets() -> None:
    """Regenerate and save all possible trinets."""
    logging.debug("Regenerating all possible trinets and saving them.")

    all_generators = {
        0: [generators.generator_level0]
        , 1: [generators.generator_level1]
        , 2: [generators.generator_A, generators.generator_B, generators.generator_C, generators.generator_D]
    }

    trinet_info_list = TrinetInfoList()
    for level, generator_list in all_generators.items():
        for generator in generator_list:
            generator_trinet_info_list = generator.build_trinets()
            trinet_info_list += generator_trinet_info_list

    pickle_out = open("data/all_trinets_save.pickle", 'wb')
    data = [all_generators, trinet_info_list]
    pickle.dump(data, pickle_out)
    pickle_out.close()


def get_trinets() -> (list, TrinetInfoList):
    """Read and retrieve all possible trinets."""
    logging.debug("Reading and retrieving all possible trinets.")
    pickle_in = open("data/all_trinets_save.pickle", 'rb')
    result = pickle.load(pickle_in)
    pickle_in.close()
    all_generators, trinet_info_list = result[0], result[1]
    return all_generators, trinet_info_list
