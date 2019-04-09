import data.test_networks as test_networks
from datastructures.rooted_level_k_network import *
from datastructures.trinet_set import *
import data.generators as generator
from data.all_trinets import *
from utils.help_functions import *
from solver import *
import time
import pickle
import logging


filename_save_trinets_B = "data/trinets_B"

logging.basicConfig(level=logging.DEBUG)

os.environ["PATH"] += os.pathsep + 'C:/Program Files (x86)/Graphviz2.38/bin/'

# regenerate_trinets()
all_generators, all_trinets, all_trinets_gen_sides = get_trinets()

# dct = test_networks.B
# network = RootedLevelKNetwork.from_connections_dict(dct)
# trinets = network.get_exhibited_trinets()
# pickle_save(filename_save_trinets_B, trinets)
trinets = pickle_read(filename_save_trinets_B)
trinet_set = TrinetSet.from_triplet_trinet_list(trinets)
solver = Solver(trinet_set)
solver.state()

