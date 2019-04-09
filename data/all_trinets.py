import data.generators as generators
from datastructures.rooted_level_k_network import *
import pickle


def regenerate_trinets() -> None:
    """Regenerate and save all possible trinets."""
    logging.debug("Regenerating all possible trinets and saving them.")
    all_generators = [generators.generator_level0, generators.generator_level1, generators.generator_A, generators.generator_B, generators.generator_C, generators.generator_D]
    all_trinets_gen_sides = []
    for generator in all_generators:
        all_trinets_gen_sides += generator.build_trinets()

    pickle_out = open("data/all_trinets_save.pickle", 'wb')
    data = [all_generators, all_trinets_gen_sides]
    pickle.dump(data,  pickle_out)
    pickle_out.close()


def get_trinets() -> (list, list, list):
    """Read and retrieve all possible trinets."""
    logging.debug("Reading and retrieving all possible trinets.")
    pickle_in = open("data/all_trinets_save.pickle", 'rb')
    result = pickle.load(pickle_in)
    pickle_in.close()
    all_generators, all_trinets_gen_sides = result[0], result[1]
    all_trinets = [x[0] for x in all_trinets_gen_sides]
    return all_generators, all_trinets, all_trinets_gen_sides
