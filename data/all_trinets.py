import data.generators as generators
from utils.help_functions import trinets_of_generator
from datastructures.rooted_level_k_network import *
import pickle
import os


def regenerate_trinets():
    all_generators_dir_adj_matrix = [generators.A, generators.B, generators.C, generators.D, generators.cactus]
    all_trinets_gen_sides = []
    all_generators = []
    for current_generator_dir_adj_matrix in all_generators_dir_adj_matrix:
        current_generator = network_from_dir_adj_matrix(current_generator_dir_adj_matrix)
        all_generators.append(current_generator)
        all_trinets_gen_sides += trinets_of_generator(current_generator)

    pickle_out = open("data/all_trinets_save.pickle", 'wb')
    data = [all_generators, all_trinets_gen_sides]
    pickle.dump(data,  pickle_out)
    pickle_out.close()


def get_trinets():
    pickle_in = open("data/all_trinets_save.pickle", 'rb')
    result = pickle.load(pickle_in)
    pickle_in.close()
    all_generators, all_trinets_gen_sides = result[0], result[1]
    all_trinets = [x[0] for x in all_trinets_gen_sides]
    return all_generators, all_trinets, all_trinets_gen_sides
