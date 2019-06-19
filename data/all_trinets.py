import data.generators as generators
from datastructures.rooted_level_k_network import *
from bidict import bidict
import pickle


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
        0 : [generators.generator_level0]
        , 1: [generators.generator_level1]
        , 2: [generators.generator_A, generators.generator_B, generators.generator_C, generators.generator_D]
    }
    # 1, 2, 15, 6, 2, 6
    # 0 - 1 - 3 - 18 - 24 - 26 - 32
    trinet_lookup_dict = {}
    for level, generator_list in all_generators.items():
        for generator in generator_list:
            trinet_info = generator.build_trinets()
            for trinet, info in trinet_info.items():
                trinet_lookup_dict[trinet] = {'generator': generator, 'on_edges': info['on_edges']}

    pickle_out = open("data/all_trinets_save.pickle", 'wb')
    data = [all_generators, trinet_lookup_dict]
    pickle.dump(data, pickle_out)
    pickle_out.close()


def get_trinets() -> (list, list, list):
    """Read and retrieve all possible trinets."""
    logging.debug("Reading and retrieving all possible trinets.")
    pickle_in = open("data/all_trinets_save.pickle", 'rb')
    result = pickle.load(pickle_in)
    pickle_in.close()
    all_generators, trinet_lookup_dict = result[0], result[1]
    return all_generators, trinet_lookup_dict
