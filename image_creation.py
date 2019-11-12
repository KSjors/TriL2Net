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
from datastructures.rooted_level_k_network import RootedLevelKNetwork, NetworkInfoList
from solver import Solver
import data.generators as generator
from data.all_trinets import get_standard_networks, regenerate_standard_networks, pickle_save, pickle_read
from utils.help_functions import *
import timeit
from tqdm import tqdm
import matplotlib.pyplot as plt
import scipy.stats as stats
from mip import model
from data.ETS_network import ETS_NETWORK_dict


enewick = "(c, (b, c)1)0"
network = RootedLevelKNetwork.from_enewick(enewick, check_valid=False)
network.visualize(internal_node_labels=False, edge_labels=False, rankdir='LR', format='pdf', file_path='extinction_example_p2')
