import logging

# set up logging to file - see previous section for more details
logging.basicConfig(level=logging.DEBUG,
                    format='%(asctime)s %(name)-12s %(levelname)-8s %(message)s',
                    datefmt='%m-%d %H:%M',
                    filename='log.log',
                    filemode='w')
# define a Handler which writes INFO messages or higher to the sys.stderr
console = logging.StreamHandler()
console.setLevel(logging.DEBUG)
# set a format which is simpler for console use
formatter = logging.Formatter('%(name)-12s: %(levelname)-8s %(message)s')
# tell the handler to use this format
console.setFormatter(formatter)
# add the handler to the root logger
logging.getLogger('').addHandler(console)

import os
import data.test_networks as test_networks
from datastructures.rooted_level_k_network import *
from datastructures.trinet_set import *
from solver import *
import data.generators as generator
from data.all_trinets import *
from utils.help_functions import *

os.environ["PATH"] += os.pathsep + 'C:/Program Files (x86)/Graphviz2.38/bin/'
filename_save_trinets_B = 'data/trinets_B'

# regenerate_trinets()
# all_generators, all_trinets, all_trinets_gen_sides = get_trinets()

dct = test_networks.B
network = RootedLevelKNetwork.from_connections_dict(dct)
# network.visualize()
# time.sleep(0.2)
# pickle_save("data/network_B.pickle", network)
# network = pickle_read("data/network_B.pickle")
# network.visualize()

# trinets = network.get_exhibited_trinets()
# pickle_save('data/trinets_B', trinets)
trinets = pickle_read(filename_save_trinets_B)

trinet_set = TrinetSet.from_trinet_list(trinets)

solver = Solver(trinet_set)