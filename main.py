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

os.environ["PATH"] += os.pathsep + 'C:/Program Files (x86)/Graphviz2.38/bin/'
filename_save_trinets_B = 'data/trinets_B'


""" Build or load trinets """
regenerate_trinets()
all_generators, all_trinets, all_trinets_gen_sides = get_trinets()

""" Build or load network """
# dct = test_networks.connections_big
# network = RootedLevelKNetwork.from_connections_dict(dct)
# pickle_save("data/network_B.pickle", network)
# network = pickle_read("data/network_B.pickle")
# network.visualize()


""" Build or load trinets """
# trinets = network.get_exhibited_trinets()
# pickle_save('data/trinets_B', trinets)
# trinets = pickle_read('data/trinets_B')

""" Build or load trinet set """
# trinet_set = TrinetSet.from_trinet_list(trinets)


""" Play around """
i=0
# 24, 25, 26, 27, 28, 29
# 30, 31
#
for tn in all_trinets:
    print("--------------", i, "-------------")
    i += 1
    print(tn)
    print("")