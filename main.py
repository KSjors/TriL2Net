import logging

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
from datastructures.rooted_level_k_network import *
from datastructures.trinet_set import *
from solver import *
import data.generators as generator
from data.all_trinets import *
from utils.help_functions import *

rebuild = {
    'generators': 0
    , 'network':  0
    , 'trinets':  0
}

os.environ["PATH"] += os.pathsep + 'C:/Program Files (x86)/Graphviz2.38/bin/'
filename_save_trinets_B = 'data/trinets_B'

""" Build or load trinets """
if rebuild['generators']:
    regenerate_trinets()
all_generators, trinet_lookup_dict = get_trinets()

for trinet, dct in trinet_lookup_dict.items():
    if dct['generator'].level == 1:
        trinet.visualize()

""" Build or load network """
if rebuild['network']:
    dct = test_networks.connections_big
    network = RootedLevelKNetwork.from_connections_dict(dct)
    pickle_save("data/network_B.pickle", network)
network = pickle_read("data/network_B.pickle")

""" Build or load trinets """
if rebuild['trinets'] or rebuild['network']:
    trinets = network.get_exhibited_trinets()
    pickle_save('data/trinets_B', trinets)
trinets = pickle_read('data/trinets_B')

""" Build or load trinet set """
trinet_set = TrinetSet.from_trinet_list(trinets)
# trinet_set.remove_trinet(["J", 'H', 'G'])

""" Solver """
solver = Solver(trinet_set)

""" Play around """

network.visualize()

a1 = solver.solve()
